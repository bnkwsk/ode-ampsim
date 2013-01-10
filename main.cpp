#include <iostream>
#include <cmath>
#include <flens/flens.cxx>
#include <sndfile.hh>
#include <sndfile.h>

using namespace std;
using namespace flens;

typedef GeMatrix<FullStorage<double> > MatrixType;
typedef DenseVector<Array<double> > VectorType;

// 12AX7 parameters
const double mi = 100.0f;
const double E_x = 1.4f;
const double K_g1 = 1060.0f;
const double K_p = 600.0f;
const double K_vb = 300.0f;
const double g_cf = 0.00001;
const double g_co = -0.2;

// circuit parameters
const double r_1 = 68000; const double g_1 = 1.0f / r_1;
const double r_g1 = 1000000; const double g_g1 = 1.0f / r_g1;
const double r_k1 = 2700; const double g_k1 = 1.0f / r_k1;
const double r_a1 = 100000; const double g_a1 = 1.0f / r_a1;
const double r_2 = 470000; const double g_2 = 1.0f / r_2;
const double r_g2 = 1000000; const double g_g2 = 1.0f / r_g2;
const double r_k2 = 1800; const double g_k2 = 1.0f / r_k2;
const double r_a2 = 100000; const double g_a2 = 1.0f / r_a2;
const double r_l = 4000000; const double g_l = 1.0f / r_l;
const double c_1 = 0.000001f;
const double c_2 = 0.000000022f;
const double u_n = 400;
const int equation_count = 7;

// simulation parameters
const int format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
const int channels = 1;
const double f_s = 44100;
const double t_s = 1.0 / f_s;
const int buffer_seconds = 1;
const int buffer_size = int(f_s) * buffer_seconds;
const int process_seconds = 30;
const char* in_file_name = "blues_in.wav";
const char* out_file_name = "out.wav";

// ------------------------------------------------------------------------------------------------
// - TRIODE MODEL                                                                                 -
// ------------------------------------------------------------------------------------------------
// sign function
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// grid current
double i_g(double u_gk)
{
    if(u_gk >= g_co)
        return g_cf * pow(u_gk - g_co, 1.5f);
    return 0;
}

// triode plate current
double i_a(double u_ak, double u_gk)
{
    double E_1 =  u_ak / K_p * log(1.0f + exp(K_p * (1.0f / mi + u_gk / sqrt(K_vb + pow(u_ak, 2)))));
    return pow(E_1, E_x) / K_g1 * (1.0f + sgn(E_1));
}
// ------------------------------------------------------------------------------------------------
// - THE END OF TRIODE MODEL                                                                      -
// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// - CIRCUIT MODEL                                                                                -
// ------------------------------------------------------------------------------------------------
// returns the value of a function present in the equation with index i
double f(int i, VectorType &x, double u_in, VectorType &u_cp)
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

    if(i == 1)
        return (u_in - u_g1) * g_1 - u_g1 * g_g1 - i_g1;
    if(i == 2)
        return u_c1p - u_c1 - (u_k1 * g_k1 - i_a1 - i_g1) / (c_1 * f_s);
    if(i == 3)
        return (u_n - u_a1) * g_a1 - (u_2 - u_g2) * g_2 - i_a1;
    if(i == 4)
        return u_c2p - u_c2 + (u_2 - u_g2) * g_2 / (c_2 * f_s);
    if(i == 5)
        return (u_2 - u_g2) * g_2 - u_g2 * g_g2 - i_g2;
    if(i == 6)
        return u_k2 * g_k2 - i_g2 - i_a2;
    if(i == 7)
        return (u_n - u_a2) * g_a2 - u_a2 * g_l - i_a2;
}

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

// returns a derivative by x_i of the a function present in the equation with index i
double df(int equation, int x_i, VectorType &x, double u_in, VectorType &u_cp)
{
    long double h = 0.0000001f;
    VectorType x_h = x;
    x_h(x_i) += h;
    return double((long double)(f(equation, x_h, u_in, u_cp) - f(equation, x, u_in, u_cp)) / h);
}
double df_bias(int triode, int equation, int x_i, VectorType &x)
{
    long double h = 0.0000001f;
    VectorType x_h = x;
    x_h(x_i) += h;
    return double((long double)(f_bias(triode, equation, x_h) - f_bias(triode, equation, x)) / h);
}

// calculates a Jacobian matrix
void J(VectorType &x, double u_in, VectorType u_cp, MatrixType &j)
{
    for(int row = 1; row <= equation_count; ++row)
        for(int column = 1; column <= equation_count; ++column)
        {
            j(row, column) = df(row, column, x, u_in, u_cp);
        }
}
void J_bias(int triode, VectorType &x, MatrixType &j)
{
    for(int row = 1; row <= 2; ++row)
        for(int column = 1; column <= 2; ++column)
        {
            j(row, column) = df_bias(triode, row, column, x);
        }
}

bool is_last_iteration(VectorType x1, VectorType x2, double epsilon = 0.00001f)
{
    for(int i = 1; i <= x1.length(); ++i)
        if(abs(x1(i)- x2(i)) > epsilon)
            return false;
    return true;
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
        blas::mv(NoTrans, -1.0f, Jacobi_bias, F_bias, 1.0f, x_bias);

        if(is_last_iteration(x_prev_bias, x_bias))
            break;
    }

    return x_bias;
}

int main()
{
    DenseVector<Array<int> > piv(equation_count);
    VectorType F(equation_count);
    MatrixType Jacobi(equation_count, equation_count);
    VectorType x(equation_count);
    VectorType x_0(equation_count);
    VectorType u_cp(2);
    VectorType x_prev(equation_count);

    SndfileHandle in_file(in_file_name, SFM_READ, format, channels, int(f_s));
    SndfileHandle out_file(out_file_name, SFM_WRITE, format, channels, int(f_s));

    if(!in_file || !out_file)
        return -1;

    VectorType x_bias_1 = get_bias(1),
               x_bias_2 = get_bias(2);

    x_0 = 0.0f, x_bias_1(1), x_bias_1(2), 0.0f,
          0.0f, x_bias_2(1), x_bias_2(2);

    u_cp = x_bias_1(1), x_bias_1(2);

    float u_in;
    float u_out;

    x = x_0;

    for(int sample = 0; sample < f_s * process_seconds; ++sample)
    {
        in_file.readf(&u_in, 1);
        for(int iteration = 0; ; ++iteration)
        {
            // get the Jacobi's matrix
            J(x, u_in, u_cp, Jacobi);

            // inverse the Jacobi's matrix using LU decomposition
            lapack::trf(Jacobi, piv);
            lapack::tri(Jacobi, piv);

            // get the F vector
            for(int equation = 1; equation <= equation_count; ++equation)
                F(equation) = f(equation, x, u_in, u_cp);

            x_prev = x;
            blas::mv(NoTrans, -1.0f, Jacobi, F, 1.0f, x);

            if(is_last_iteration(x_prev, x))
                break;
        }

        u_out = x(4) / 25.0f;

        double u_k1 = x(2),
               u_a1 = x(3),
               u_2 = x(4),
               u_c1 = u_k1,
               u_c2 = u_a1 - u_2;

        u_cp(1) = u_k1;
        u_cp(2) = u_c2;

        //cout << sample << ',' << 10 * u_in << ',' << u_out << endl;

        if(sample % 3000 == 0)
            cout << 100.0f * (1.0f * sample) / (1.0f * f_s * process_seconds) << "%..." << endl;
        out_file.write(&u_out, 1);
    }
}
