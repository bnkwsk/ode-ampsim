#pragma once

#include <cmath>

#include "Sample.hpp"
#include "ProcessorElement.hpp"
#include "CircuitBlock.hpp"
#include "PipeliningProcessor.hpp"
#include "CabinetSimulator.hpp"
#include "BilinearInterpolator.hpp"

// 12AX7 triode parameters
const double mi = 100.0;
const double E_x = 1.4;
const double K_g1 = 1060.0;
const double K_p = 600.0;
const double K_vb = 300.0;
const double g_cf = 0.00001;
const double g_co = -0.2;

class GuitarPreampBase
{
protected:
    // circuit parameters
    const double r_1 = 68e3; const double g_1 = 1.0f / r_1;
    const double r_g1 = 1e6; const double g_g1 = 1.0f / r_g1;
    const double r_k1 = 2700; const double g_k1 = 1.0f / r_k1;
    const double r_a1 = 100e3; const double g_a1 = 1.0f / r_a1;
    const double r_2 = 470e3; const double g_2 = 1.0f / r_2; // 1e6
    const double r_g2 = 1e6; const double g_g2 = 1.0f / r_g2; // 470e3
    const double r_k2 = 1800; const double g_k2 = 1.0f / r_k2;
    const double r_a2 = 100e3; const double g_a2 = 1.0f / r_a2;
    const double r_3 = 470e3; const double g_3 = 1.0f / r_3;
    const double r_g3 = 470e3; const double g_g3 = 1.0f / r_g3;
    const double r_k3 = 1800; const double g_k3 = 1.0f / r_k3;
    const double r_a3 = 100e3; const double g_a3 = 1.0f / r_a3;
    const double r_4 = 470e3; const double g_4 = 1.0f / r_4;
    const double r_g4 = 470e3; const double g_g4 = 1.0f / r_g4;
    const double r_k4 = 1800; const double g_k4 = 1.0f / r_k4;
    const double r_a4 = 100e3; const double g_a4 = 1.0f / r_a4;
    const double r_l = 4e6; const double g_l = 1.0f / r_l;
    const double c_1 = 1e-6;
    const double c_2 = 22e-9;
    const double c_3 = 22e-9;
    const double c_4 = 22e-9;
    const double c_5 = 22e-9;
    const double u_n = 400;

    double outDivider;
    VectorType xBiasStart;
    int equationCount, capacitorCount;

    // grid current
    inline double i_g(double u_gk)
    {
        return (u_gk >= g_co) ? g_cf * pow(u_gk - g_co, 1.5) : 0;
    }

    // triode plate current
    static BilinearInterpolator i_a;

public:
    GuitarPreampBase(int _equationCount, int _capacitorCount, double _outDivider=1.0) : equationCount(_equationCount), capacitorCount(_capacitorCount), outDivider(_outDivider)
    {
        xBiasStart = VectorType(equationCount);
    }

    VectorType getXBiasStart()
    {
        return xBiasStart;
    }

    int getEquationCount()
    {
        return equationCount;
    }
    int getCapacitorCount()
    {
        return capacitorCount;
    }
};

// sign function    
template <typename T> 
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
// anode current function
double _i_a(double u_ak, double u_gk)
{
    double E_1 = u_ak / K_p * log(1.0 + exp(K_p * (1.0 / mi + u_gk / sqrt(K_vb + u_ak * u_ak))));
    return pow(E_1, E_x) / K_g1 * (1.0 + sgn(E_1));
}
BilinearInterpolator GuitarPreampBase::i_a = BilinearInterpolator(_i_a, -2000.0, 2000.0, -200.0, 200.0, 5.0, 0.05);

class GuitarPreampFirstBlockCircuit : public GuitarPreampBase
{
public:
    void resetState(VectorType x_bias, VectorType &uCP, VectorType &iCP)
    {
        uCP = x_bias(2), x_bias(3) - x_bias(4), x_bias(7) - x_bias(8);
        iCP = 0.0, 0.0, 0.0;
    }

    VectorType f(VectorType &x, double u_in, VectorType &uCP, VectorType &iCP, double fS) //__attribute__((always_inline))
    {
        double u_g1 = x(1),
            u_k1 = x(2),
            u_a1 = x(3),
            u_2 = x(4),
            u_g2 = x(5),
            u_k2 = x(6),
            u_a2 = x(7),
            u_3 = x(8),
            u_c1 = u_k1,
            u_c2 = u_a1 - u_2,
            u_c3 = u_a2 - u_3,
            u_c1p = uCP(1),
            u_c2p = uCP(2),
            u_c3p = uCP(3),
            i_c1p = iCP(1),
            i_c2p = iCP(2),
            i_c3p = iCP(3),
            i_g1 = i_g(u_g1 - u_k1),
            i_a1 = i_a(u_a1 - u_k1, u_g1 - u_k1),
            i_g2 = i_g(u_g2 - u_k2),
            i_a2 = i_a(u_a2 - u_k2, u_g2 - u_k2);

        VectorType results(equationCount);
        results =
#include "generated/FirstBlock.cpp.part"
        ;
        return results;
    }

    inline void updateState(VectorType x, VectorType &uCP, VectorType &iCP, double fS)
    {
        double u_k1 = x(2),
               u_a1 = x(3),
               u_2 = x(4),
               u_a2 = x(7),
               u_3 = x(8),
               u_c2 = u_a1 - u_2,
               u_c3 = u_a2 - u_3;

        iCP = c_1 * 2 * (u_k1 - uCP(1)) * fS - iCP(1),
              c_2 * 2 * (u_c2 - uCP(2)) * fS - iCP(2),
              c_3 * 2 * (u_c3 - uCP(3)) * fS - iCP(3);

        uCP = u_k1,
              u_c2,
              u_c3;
    }

    inline double getOutput(VectorType &x)
    {
        return x(4) / outDivider;
    }

    GuitarPreampFirstBlockCircuit() : GuitarPreampBase(8, 3)
    {
        xBiasStart = 0.0, 2.0, 300.0, 0.0, 0.0, 2.0, 300.0, 0.0;
    }
};

class GuitarPreampSecondBlockCircuit : public GuitarPreampBase
{
public:
    void resetState(VectorType x_bias, VectorType &uCP, VectorType &iCP)
    {
        uCP = x_bias(3) - x_bias(4), x_bias(7) - x_bias(8);
        iCP = 0.0, 0.0;
    }

    VectorType f(VectorType &x, double u_in, VectorType &uCP, VectorType &iCP, double fS) //__attribute__((always_inline))
    {
        double u_g2 = x(1),
            u_k2 = x(2),
            u_a2 = x(3),
            u_3 = x(4),
            u_g3 = x(5),
            u_k3 = x(6),
            u_a3 = x(7),
            u_4 = x(8),
            u_c3 = u_a2 - u_3,
            u_c4 = u_a3 - u_4,
            u_c3p = uCP(1),
            u_c4p = uCP(2),
            i_c3p = iCP(1),
            i_c4p = iCP(2),
            i_g2 = i_g(u_g2 - u_k2),
            i_a2 = i_a(u_a2 - u_k2, u_g2 - u_k2),
            i_g3 = i_g(u_g3 - u_k3),
            i_a3 = i_a(u_a3 - u_k3, u_g3 - u_k3);

        VectorType results(equationCount);
        results =
#include "generated/SecondBlock.cpp.part"
        ;
        return results;
    }

    inline void updateState(VectorType x, VectorType &uCP, VectorType &iCP, double fS)
    {
        double u_k2 = x(2),
               u_a2 = x(3),
               u_3 = x(4),
               u_a3 = x(7),
               u_4 = x(8),
               u_c3 = u_a2 - u_3,
               u_c4 = u_a3 - u_4;

        iCP = c_3 * 2 * (u_c3 - uCP(1)) * fS - iCP(1),
               c_4 * 2 * (u_c4 - uCP(2)) * fS - iCP(2);

        uCP = u_c3, u_c4;
    }

    inline double getOutput(VectorType &x)
    {
        return x(4) / outDivider;
    }
    GuitarPreampSecondBlockCircuit() : GuitarPreampBase(8, 2)
    {
        xBiasStart = 0.0, 2.0, 300.0, 0.0, 0.0, 2.0, 300.0, 0.0;
    }
};


class GuitarPreampThirdBlockCircuit : public GuitarPreampBase
{
public:
    void resetState(VectorType x_bias, VectorType &uCP, VectorType &iCP)
    {
        uCP = x_bias(3) - x_bias(4), x_bias(7) - x_bias(8);
        iCP = 0.0, 0.0;
    }

    VectorType f(VectorType &x, double u_in, VectorType &uCP, VectorType &iCP, double fS) //__attribute__((always_inline))
    {
        double u_g3 = x(1),
            u_k3 = x(2),
            u_a3 = x(3),
            u_4 = x(4),
            u_g4 = x(5),
            u_k4 = x(6),
            u_a4 = x(7),
            u_5 = x(8),
            u_c4 = u_a3 - u_4,
            u_c5 = u_a4 - u_5,
            u_c4p = uCP(1),
            u_c5p = uCP(2),
            i_c4p = iCP(1),
            i_c5p = iCP(2),
            i_g3 = i_g(u_g3 - u_k3),
            i_a3 = i_a(u_a3 - u_k3, u_g3 - u_k3),
            i_g4 = i_g(u_g4 - u_k4),
            i_a4 = i_a(u_a4 - u_k4, u_g4 - u_k4);

        VectorType results(equationCount);
        results =
#include "generated/ThirdBlock.cpp.part"
        ;
        return results;
    }

    inline void updateState(VectorType x, VectorType &uCP, VectorType &iCP, double fS)
    {
        double u_k3 = x(2),
               u_a3 = x(3),
               u_4 = x(4),
               u_a4 = x(7),
               u_5 = x(8),
               u_c4 = u_a3 - u_4,
               u_c5 = u_a4 - u_5;

        iCP = c_4 * 2 * (u_c4 - uCP(1)) * fS - iCP(1),
               c_5 * 2 * (u_c5 - uCP(2)) * fS - iCP(2);

        uCP = u_c4, u_c5;
    }

    inline double getOutput(VectorType &x)
    {
        return x(8) / outDivider;
    }

    GuitarPreampThirdBlockCircuit() : GuitarPreampBase(8, 2, 300.0)
    {
        xBiasStart = 0.0, 2.0, 300.0, 0.0, 0.0, 2.0, 300.0, 0.0;
    }
};

typedef PipeliningProcessor<Sample, ProcessorElement> ProcessorType;

class GuitarAmp : public ProcessorType
{
    double fS;
    ProcessorElement *firstBlock;
    ProcessorElement *secondBlock;
    ProcessorElement *thirdBlock;
    CabinetSimulator *cabinet;

public:
    GuitarAmp(double _fS, const char *impulseResponsePath) : PipeliningProcessor(), fS(_fS)
    {
        firstBlock = new CircuitBlock<GuitarPreampFirstBlockCircuit>(fS);
        secondBlock = new CircuitBlock<GuitarPreampSecondBlockCircuit>(fS);
        thirdBlock = new CircuitBlock<GuitarPreampThirdBlockCircuit>(fS);
        cabinet = new CabinetSimulator(impulseResponsePath);
        elements.push_back(firstBlock);
        elements.push_back(secondBlock);
        elements.push_back(thirdBlock);
        elements.push_back(cabinet);
        setup();
    }

    double getFS()
    {
        return fS;
    }
};