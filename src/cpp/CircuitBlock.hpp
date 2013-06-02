#pragma once

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

template<typename Circuit>
class CircuitBlock;

class NoConvergenceException
{
    ProcessorElement &block;
public:
    NoConvergenceException(ProcessorElement &b) : block(b)
    {
        //
    }

    ProcessorElement &getBlock()
    {
        return block;
    }
};

template<typename Circuit>
class CircuitBlock : public ProcessorElement
{
private:
    IntegerVectorType piv;
    MatrixType Jacobi;
    VectorType x;
    VectorType xBackup;
    VectorType xPrev;
    VectorType xPrevBackup;
    VectorType xPredict;
    VectorType xPrevPrev;
    VectorType uCP;
    VectorType iCP;
    VectorType xPlus;
    VectorType xMinus;
    Circuit circuit;
    int capacitorCount;
    double fS;
    int equationCount;
    const flens::Underscore<VectorType::IndexType> _;

    // calculate the Jacobi's matrix
    void J(VectorType &x, double uIn, VectorType &uCP, VectorType &iCP)
    {
        double h = 1e-10;
        for(int i = 1; i <= equationCount; ++i)
        {
            xPlus = xMinus = x;
            xPlus(i) += h;
            xMinus(i) -= h;
            Jacobi(_, i) =
                (circuit.f(xPlus, uIn, uCP, iCP, fS) - circuit.f(xMinus, uIn, uCP, iCP, fS)) /
                (2.0 * h);
        }
    }

    // calculate the Jacobi's matrix for getting the bias
    void J_bias(VectorType &x, MatrixType &J_b)
    {
        double h = 1e-10;
        VectorType xPlus(equationCount);
        VectorType xMinus(equationCount);
        for(int i = 1; i <= equationCount; ++i)
        {
            xPlus = xMinus = x;
            xPlus(i) += h;
            xMinus(i) -= h;
            J_b(_, i) = (f_bias(xPlus) - f_bias(xMinus)) / (2.0 * h);
        }
    }

    VectorType f_bias(VectorType &x) {
        VectorType iCP(capacitorCount);
        VectorType uCP(capacitorCount);
        for(unsigned int i = 1; i <= capacitorCount; ++i)
            iCP(i) = uCP(i) = 0.0;
        return circuit.f(x, 0.0, uCP, iCP, fS);
    };

    bool inline isLastIteration(VectorType x1, VectorType x2, double epsilon=1e-4)
    {
        for(int i = 1; i <= x1.length(); ++i)
        {
            if(x1(i) - x2(i) > epsilon || x1(i) - x2(i) < -epsilon)
                return false;
        }
        return true;
    }

    double process(double uIn)
    {
        do
        {
            // save the previous state vector
            xPrev = x;

            // calculate the Jacobi's matrix
            J(x, uIn, uCP, iCP);

            // the x vector becomes the rgiht hand side value of the equation J * x(t+1) = J * x(t) - f(x(t))
            x = Jacobi * xPrev - circuit.f(xPrev, uIn, uCP, iCP, fS);

            // solve the linear system to get the new state vector value
            flens::lapack::sv(Jacobi, piv, x);
        } while(!isLastIteration(xPrev, x));
     
        // test the convergence
        for(int i = 1; i <= equationCount; ++i)
            if(x(i) != x(i))
                throw NoConvergenceException(*this);

        circuit.updateState(x, uCP, iCP, fS);
        
        return circuit.getOutput(x);
    }

    VectorType getBias(VectorType xBiasStart)
    {
        IntegerVectorType piv(equationCount);
        VectorType x_bias(equationCount);
        VectorType x_prev_bias(equationCount);
        VectorType uCP(capacitorCount);
        MatrixType Jacobi_bias(equationCount, equationCount);
        double _fS = fS;
        fS = 0.0;
        x_bias = xBiasStart;

        do
        {
            x_prev_bias = x_bias;

            J_bias(x_bias, Jacobi_bias);
            x_bias = Jacobi_bias * x_prev_bias - f_bias(x_prev_bias);
            flens::lapack::sv(Jacobi_bias, piv, x_bias);
        } while(!isLastIteration(x_prev_bias, x_bias));
        
        fS = _fS;
        return x_bias;
    }
public:
    CircuitBlock(double _fS) : ProcessorElement(), fS(_fS)
    {
        equationCount = circuit.getEquationCount();
        capacitorCount = circuit.getCapacitorCount();
        piv = IntegerVectorType(equationCount);
        x = VectorType(equationCount);
        xBackup = VectorType(equationCount);
        xPrev = VectorType(equationCount);
        xPrevBackup = VectorType(equationCount);
        xPredict = VectorType(equationCount);
        xPrevPrev = VectorType(equationCount);
        xPlus = VectorType(equationCount);
        xMinus = VectorType(equationCount);
        uCP = VectorType(capacitorCount);
        iCP = VectorType(capacitorCount);
        Jacobi = MatrixType(equationCount, equationCount);

        x = xPrev = getBias(circuit.getXBiasStart());
        circuit.resetState(x, uCP, iCP);
    }

    void run()
    {
        bool last = false;
        Sample in;
        double uIn;
        double u_in_prev; // as the left hand side of the interpolation
        double u_in_inter; // as the left hand side of the interpolation
        double u_out;
        double originalFS = fS;
        double maxScalingFactor = 512;

        double scalingFactor; // splitted quantity
        while(!last)
        {
            in = input->pop();
            u_in_prev = uIn;
            uIn = in.getValue();
            last = in.isLast();
            xBackup = x; // before anything goes wrong
            xPrevBackup = xPrev; // before anything goes wrong
            fS = originalFS;
            try
            {
                u_out = process(uIn);
                output->push(Sample(u_out, last));
            }
            catch(NoConvergenceException e1)
            {
                for(scalingFactor = 2.0; scalingFactor <= maxScalingFactor; scalingFactor *= 2.0)
                {
                    fS = scalingFactor * originalFS;
                    try {
                        x = xBackup;
                        xPrev = xPrevBackup;
                        for(double i = 1; i <= scalingFactor; i += 1.0)
                        {
                            u_in_inter = i / scalingFactor * uIn + (1.0 - i / scalingFactor) * u_in_prev;
                            u_out = process(u_in_inter);
                        }
                        output->push(Sample(u_out, last));
                        break;
                    }
                    catch(NoConvergenceException e2)
                    {
                        // if the upsampling wasn't successfull return the previous value
                        if(scalingFactor == maxScalingFactor)
                            output->push(Sample(circuit.getOutput(xBackup), last));
                    }
                }
            }
        }
    }
};