#include <iostream>
#include <cmath>
#include <flens/flens.cxx>
#include <sndfile.hh>
#include <sndfile.h>

using namespace std;
using namespace flens;

typedef GeMatrix<FullStorage<double>> MatrixType;
typedef DenseVector<Array<double>> VectorType;
const Underscore<VectorType::IndexType> _;

// 12AX7 parameters
const double mi = 100.0;
const double E_x = 1.4;
const double K_g1 = 1060.0;
const double K_p = 600.0;
const double K_vb = 300.0;
const double g_cf = 0.00001;
const double g_co = -0.2;

// circuit parameters
const double r_1 = 68e3; const double g_1 = 1.0f / r_1;
const double r_g1 = 1e6; const double g_g1 = 1.0f / r_g1;
const double r_k1 = 2.7e3; const double g_k1 = 1.0f / r_k1;
const double r_a1 = 100e3; const double g_a1 = 1.0f / r_a1;
const double r_2 = 470e3; const double g_2 = 1.0f / r_2;
const double r_g2 = 1e6; const double g_g2 = 1.0f / r_g2;
const double r_k2 = 1.8e3; const double g_k2 = 1.0f / r_k2;
const double r_a2 = 100e3; const double g_a2 = 1.0f / r_a2;
const double r_l = 4e6; const double g_l = 1.0f / r_l;
const double c_1 = 1e-6;
const double c_2 = 22e-9;
const double u_n = 400;

// simulation parameters
const int format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
const int channels = 1;
double f_s;
const double t_s = 1.0 / f_s;
const int buffer_seconds = 1;
const int buffer_size = int(f_s) * buffer_seconds;
int process_seconds = 10;
const char* in_file_name;
const char* out_file_name;

// ----------------------------------------------------------------------------
// - TRIODE MODEL                                                             -
// ----------------------------------------------------------------------------
// sign function
template <typename T> 
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// grid current
inline double i_g(double u_gk)
{
    return (u_gk >= g_co) ? g_cf * pow(u_gk - g_co, 1.5) : 0;
}

// triode plate current
inline double i_a(double u_ak, double u_gk)
{
    double E_1 = u_ak / K_p * log(1.0 + exp(K_p * (1.0 / mi + u_gk / sqrt(K_vb + u_ak * u_ak))));
    return pow(E_1, E_x) / K_g1 * (1.0 + sgn(E_1));
}
// ------------------------------------------------------------------------------------------------
// - THE END OF TRIODE MODEL                                                                      -
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// - CIRCUIT MODEL                                                                                -
// ------------------------------------------------------------------------------------------------
// returns the value of a function present in the equation with index i
double f_bias(int triode, int i, VectorType &x)
{
    double u_k1 = x(1),
           u_a1 = x(2),
           i_g1 = i_g(0 - u_k1),
           i_a1 = i_a(u_a1 - u_k1, 0 - u_k1),
           g_k,
           g_a;

    if(triode == 1)
    {
        g_k = g_k1;
        g_a = g_a1;
    }
    if(triode == 2)
    {
        g_k = g_k2;
        g_a = g_a2;
    }

    if(i == 1)
        return u_k1 * g_k - i_a1 - i_g1;
    if(i == 2)
        return (u_n - u_a1) * g_a - i_a1;
}
// ------------------------------------------------------------------------------------------------
// - THE END OF CIRCUIT MODEL                                                                                -
// ------------------------------------------------------------------------------------------------

// calculates a Jacobian matrix
double df_bias(int triode, int equation, int x_i, VectorType &x)
{
    long double h = 0.00001;
    VectorType x_plus = x, x_minus = x;
    x_plus(x_i) += h;
    x_minus(x_i) -= h;
    return double((long double)(f_bias(triode, equation, x_plus) - f_bias(triode, equation, x_minus)) / (2.0f * h));
}
void J_bias(int triode, VectorType &x, MatrixType &j)
{
    for(int row = 1; row <= 2; ++row)
        for(int column = 1; column <= 2; ++column)
        {
            j(row, column) = df_bias(triode, row, column, x);
        }
}

class Block
{
protected:
    DenseVector<Array<int> > piv;
    VectorType F;
    VectorType F_prev;
    MatrixType Jacobi;
    MatrixType Jacobi_prev;
    VectorType x;
    VectorType x_prev;
    VectorType u_cp;
    VectorType x_plus;
    VectorType x_minus;
    float out_divider;
    long long iterated;
    int equation_count;

    // calculate the Jacobian matrix
    void J(VectorType &x, double u_in, VectorType &u_cp)
    {
        double h = 1e-10;
        for(int i = 1; i <= equation_count; ++i)
        {
            x_plus = x_minus = x;
            x_plus(i) += h;
            x_minus(i) -= h;
            Jacobi(_, i) = (f(x_plus, u_in, u_cp) - f(x_minus, u_in, u_cp))
                           / (2.0 * h);
        }
    }

    void resetState()
    {
        VectorType x_bias_1 = get_bias(1),
                   x_bias_2 = get_bias(2);
        x = 0.0, x_bias_1(1), x_bias_1(2), 0.0,
            0.0, x_bias_2(1), x_bias_2(2);
        u_cp = x_bias_1(1), x_bias_1(2);
    }

public:
    Block(float out_div=1.0, int e_q=0) : out_divider(out_div), iterated(0), equation_count(e_q)
    {
        piv = DenseVector<Array<int> >(equation_count);
        x = VectorType(equation_count);
        x_prev = VectorType(equation_count);
        x_plus = VectorType(equation_count);
        x_minus = VectorType(equation_count);
        u_cp = VectorType(2);
        Jacobi = MatrixType(equation_count, equation_count);
        F = VectorType(equation_count);
        resetState();
    }

    virtual VectorType f(VectorType &x, double u_in, VectorType &u_cp) = 0;

    double process(float u_in)
    {
        for(int iteration = 0; ; ++iteration)
        {
            // save the previous state vector value
            x_prev = x;

            // calculate the Jacobi's matrix
            J(x, u_in, u_cp);

            // the x vector becomes the rgiht hand side value of the equation J * x(t+1) = J * x(t) - f(x(t))
            x = Jacobi * x_prev - f(x_prev, u_in, u_cp);

            // solve the linear system to get the new state vector value
            lapack::sv(Jacobi, piv, x);

            if(is_last_iteration(x_prev, x))
                break;
        }

        double u_k1 = x(2),
               u_a1 = x(3),
               u_2 = x(4),
               u_c2 = u_a1 - u_2;

        u_cp = u_k1,
               u_c2;

        return x(4) / out_divider;
    }

    VectorType get_bias(int triode)
    {
        DenseVector<Array<int> > piv(2);
        VectorType F_bias(2);
        VectorType x_bias(2);
        VectorType x_prev_bias(2);
        MatrixType Jacobi_bias(2, 2);

        x_bias = 2, 300;

        for(int iteration = 0; ; ++iteration)
        {
            // get the Jacobi's matrix
            J_bias(triode, x_bias, Jacobi_bias);

            // inverse the Jacobi's matrix using LU decomposition
            lapack::trf(Jacobi_bias, piv);
            lapack::tri(Jacobi_bias, piv);

            // get the F vector
            for(int equation = 1; equation <= 2; ++equation)
            F_bias(equation) = f_bias(triode, equation, x_bias);

            x_prev_bias = x_bias;
            blas::mv(NoTrans, -1.0, Jacobi_bias, F_bias, 1.0, x_bias);

            if(is_last_iteration(x_prev_bias, x_bias))
                break;
        }

        return x_bias;
    }

    bool is_last_iteration(VectorType x1, VectorType x2, double epsilon = 1e-5)
    {
        for(int i = 1; i <= x1.length(); ++i)
            if(abs(x1(i) - x2(i)) > epsilon)
                return false;
        return true;
    }

};

class FirstBlock : public Block {
    VectorType f(VectorType &x, double u_in, VectorType &u_cp)
    {
        double u_g1 = x(1),
            u_k1 = x(2),
            u_a1 = x(3),
            u_2 = x(4),
            u_g2 = x(5),
            u_k2 = x(6),
            u_a2 = x(7),
            u_c1 = u_k1,
            u_c2 = u_a1 - u_2,
            u_c1p = u_cp(1),
            u_c2p = u_cp(2),
            i_g1 = i_g(u_g1 - u_k1),
            i_a1 = i_a(u_a1 - u_k1, u_g1 - u_k1),
            i_g2 = i_g(u_g2 - u_k2),
            i_a2 = i_a(u_a2 - u_k2, u_g2 - u_k2);

        VectorType results(equation_count);
        results =
            (u_in - u_g1) * g_1 - u_g1 * g_g1 - i_g1,
            u_c1p - u_c1 - (u_k1 * g_k1 - i_a1 - i_g1) / (c_1 * f_s),
            (u_n - u_a1) * g_a1 - (u_2 - u_g2) * g_2 - i_a1,
            u_c2p - u_c2 + (u_2 - u_g2) * g_2 / (c_2 * f_s),
            (u_2 - u_g2) * g_2 - u_g2 * g_g2 - i_g2,
            u_k2 * g_k2 - i_g2 - i_a2,
            (u_n - u_a2) * g_a2 - u_a2 * g_l - i_a2;

        return results;
    }
public:
    FirstBlock(double out_div=1.0) : Block(out_div, 7)
    {
        //
    }
};

class SecondBlock : public Block {
    VectorType f(VectorType &x, double u_in, VectorType &u_cp)
    {
        double u_g1 = x(1),
            u_k1 = x(2),
            u_a1 = x(3),
            u_2 = x(4),
            u_g2 = x(5),
            u_k2 = x(6),
            u_a2 = x(7),
            u_c2 = u_a1 - u_2,
            u_c2p = u_cp(2),
            i_g1 = i_g(u_g1 - u_k1),
            i_a1 = i_a(u_a1 - u_k1, u_g1 - u_k1),
            i_g2 = i_g(u_g2 - u_k2),
            i_a2 = i_a(u_a2 - u_k2, u_g2 - u_k2),
            g_1 = g_2;

        VectorType results(equation_count);
        results =
            (u_in - u_g1) * g_1 - u_g1 * g_g1 - i_g1,
            u_k1 * g_k1 - i_a1 - i_g1,
            (u_n - u_a1) * g_a1 - (u_2 - u_g2) * g_2 - i_a1,
            u_c2p - u_c2 + (u_2 - u_g2) * g_2 / (c_2 * f_s),
            (u_2 - u_g2) * g_2 - u_g2 * g_g2 - i_g2,
            u_k2 * g_k2 - i_g2 - i_a2,
            (u_n - u_a2) * g_a2 - u_a2 * g_l - i_a2;

        return results;
    }
public:
    SecondBlock(double out_div=1.0) : Block(out_div, 7)
    {
        //
    }
};

int main(int argc, char **argv)
{

    f_s = atoi(argv[3]);
    int start_seconds = atoi(argv[4]);
    process_seconds = atoi(argv[5]);

    SndfileHandle in_file(argv[1], SFM_READ, format, channels, int(f_s));
    SndfileHandle out_file(argv[2], SFM_WRITE, format, channels, int(f_s));

    if(!in_file || !out_file)
        return -1;

    float u_in;
    double u_processed;
    float u_out;

    FirstBlock block1 = FirstBlock();
    SecondBlock block2 = SecondBlock(125.0f);

    in_file.seek(start_seconds * f_s, SEEK_SET);

    for(int sample = 0; sample < f_s * process_seconds; ++sample)
    {
        in_file.readf(&u_in, 1);

        u_processed = block1.process(u_in);
        u_out = block2.process(u_processed);

        if(sample % 10000 == 0)
          cout << sample << ',' << 10 * u_in << ',' << 125.0f * u_out << endl;

        if(sample % 10000 == 0)
            cout << 100.0f * (1.0f * sample) / (1.0f * f_s * process_seconds) << "%..." << endl;

        out_file.write(&u_out, 1);
    }

    return 0;
}
