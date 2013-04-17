#include <cmath>
#include <thread>
#include <vector>
#include <iostream>

#include <flens/flens.cxx>

#include "ProcessorElement.hpp"
#include "Sample.hpp"

typedef flens::GeMatrix<flens::FullStorage<double>> MatrixType;
typedef flens::DenseVector<flens::Array<double>> VectorType;
typedef flens::DenseVector<flens::Array<int>> IntegerVectorType;

class PreampBlock;

class NoConvergenceException
{
    PreampBlock &block;
public:
    NoConvergenceException(PreampBlock &b) : block(b)
    {
        //
    }

    PreampBlock &getBlock()
    {
        return block;
    }
};

class PreampBlock : public ProcessorElement
{
protected:
    IntegerVectorType piv;
    MatrixType Jacobi;
    VectorType x;
    VectorType x_prev;
    VectorType u_cp;
    VectorType i_cp;
    VectorType x_plus;
    VectorType x_minus;
    float out_divider;
    int equation_count;
    double f_s;

    const flens::Underscore<VectorType::IndexType> _;

    // calculate the Jacobi's matrix
    void J(VectorType &x, double u_in, VectorType &u_cp, VectorType &i_cp)
    {
        double h = 1e-10;
        for(int i = 1; i <= equation_count; ++i)
        {
            x_plus = x_minus = x;
            x_plus(i) += h;
            x_minus(i) -= h;
            Jacobi(_, i) =
                (f(x_plus, u_in, u_cp, i_cp) - f(x_minus, u_in, u_cp, i_cp)) /
                (2.0 * h);
        }
    }

    // calculate the Jacobi's matrix for getting the bias
    void J_bias(VectorType &x, MatrixType &J_b)
    {
        double h = 1e-10;
        VectorType x_plus(4);
        VectorType x_minus(4);
        for(int i = 1; i <= 4; ++i)
        {
            x_plus = x_minus = x;
            x_plus(i) += h;
            x_minus(i) -= h;
            J_b(_, i) = (f_bias(x_plus) - f_bias(x_minus)) / (2.0 * h);
        }
    }

    virtual VectorType f(VectorType &x, double u_in, VectorType &u_cp, VectorType &i_cp) {};
    virtual VectorType f_bias(VectorType &x) {};
    virtual void updateState() {};

public:
    void resetState()
    {
        VectorType x_bias = get_bias();
        x = 0.0, x_bias(1), x_bias(2), 0.0,
            0.0, x_bias(3), x_bias(4);
        u_cp = x_bias(1), x_bias(2);
        i_cp = 0.0, 0.0;
        x_prev = x;
    }

    PreampBlock(double _f_s, float out_div=1.0, int e_q=0) : ProcessorElement(), out_divider(out_div), equation_count(e_q), f_s(_f_s)
    {
        piv = IntegerVectorType(equation_count);
        x = VectorType(equation_count);
        x_prev = VectorType(equation_count);
        x_plus = VectorType(equation_count);
        x_minus = VectorType(equation_count);
        u_cp = VectorType(2);
        i_cp = VectorType(2);
        Jacobi = MatrixType(equation_count, equation_count);
    }

    void setUCp(VectorType u)
    {
        u_cp = u;
        x(2) = u(1);
        x(3) = u(2);
    }

    VectorType getState()
    {
        return x;
    }

    void run()
    {
        bool last = false;
        Sample in, out;
        double u_in;
        while(!last)
        {
            in = input->pop();
            u_in = in.getValue();
            last = in.isLast();
            output->push(Sample(process(u_in), last));
        }
    }

    double process(double u_in)
    {
        // predict the state vector value using simple extrapolation
        //x_predict = 2 * x_prev - x_prev_prev;
        //x_prev_prev = x_prev;
        //x = x_predict;

        for(int iteration = 0; ; ++iteration)
        {
            // save the previous state vector
            x_prev = x;

            // calculate the Jacobi's matrix
            J(x, u_in, u_cp, i_cp);

            // the x vector becomes the rgiht hand side value of the equation J * x(t+1) = J * x(t) - f(x(t))
            x = Jacobi * x_prev - f(x_prev, u_in, u_cp, i_cp);

            // solve the linear system to get the new state vector value
            flens::lapack::sv(Jacobi, piv, x);

            if(is_last_iteration(x_prev, x))
                break;
        }

        for(int i = 1; i <= equation_count; ++i)
            if(x(i) != x(i))
            {
                throw NoConvergenceException(*this);
            }

        updateState();
        
        return x(4) / out_divider;
    }

    VectorType get_bias()
    {
        IntegerVectorType piv(4);
        VectorType x_bias(4);
        VectorType x_prev_bias(4);
        MatrixType Jacobi_bias(4, 4);
        x_bias = 2, 300, 2, 300;

        for(int iteration = 0; ; ++iteration)
        {
            x_prev_bias = x_bias;

            // get the Jacobi's matrix
            J_bias(x_bias, Jacobi_bias);
            //J_bias(triode, x_bias, Jacobi_bias);

            //flens::lapack::trf(Jacobi_bias, piv);
            //flens::lapack::tri(Jacobi_bias, piv);
            x_bias = Jacobi_bias * x_prev_bias - f_bias(x_prev_bias);

            flens::lapack::sv(Jacobi_bias, piv, x_bias);

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