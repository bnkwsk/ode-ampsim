#include <cmath>

#include "Sample.hpp"
#include "ProcessorElement.hpp"
#include "PreampBlock.hpp"
#include "PipeliningProcessor.hpp"
#include "CabinetSimulator.hpp"


class GuitarPreampBlock : public PreampBlock
{
protected:
    // 12AX7 triode parameters
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
    const double r_l = 940e3; const double g_l = 1.0f / r_l;
    const double c_1 = 1e-6;
    const double c_2 = 22e-9;
    const double u_n = 400;

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

    virtual VectorType f_bias(VectorType &x)
    {
        double u_k1 = x(1),
               u_a1 = x(2),
               u_k2 = x(3),
               u_a2 = x(4),
               i_g1 = i_g(0 - u_k1),
               i_a1 = i_a(u_a1 - u_k1, 0 - u_k1),
               i_g2 = i_g(0 - u_k2),
               i_a2 = i_a(u_a2 - u_k2, 0 - u_k2);
        
        VectorType results(4);

        results = u_k1 * g_k1 - i_a1 - i_g1,
                  (u_n - u_a1) * g_a1 - i_a1,
                  u_k2 * g_k2 - i_a2 - i_g2,
                  (u_n - u_a2) * g_a2 - i_a2;

        return results;
    }

    void updateState()
    {
        double u_k1 = x(2),
               u_a1 = x(3),
               u_2 = x(4),
               u_c2 = u_a1 - u_2;

        i_cp = c_1 * 2 * (u_k1 - u_cp(1)) * f_s - i_cp(1),
               c_2 * 2 * (u_c2 - u_cp(2)) * f_s - i_cp(2);

        u_cp = u_k1,
               u_c2;
    }
public:
    GuitarPreampBlock(double f_s, float out_div, int e_q=0) : PreampBlock(f_s, out_div, e_q)
    {
        resetState();
    }
};


class GuitarPreampFirstBlock : public GuitarPreampBlock
{
    VectorType f(VectorType &x, double u_in, VectorType &u_cp, VectorType &i_cp)
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
            i_c1p = i_cp(1),
            i_c2p = i_cp(2),
            i_g1 = i_g(u_g1 - u_k1),
            i_a1 = i_a(u_a1 - u_k1, u_g1 - u_k1),
            i_g2 = i_g(u_g2 - u_k2),
            i_a2 = i_a(u_a2 - u_k2, u_g2 - u_k2);

        VectorType results(equation_count);
        results =
#include "generated/FirstBlock.cpp.part"
        ;
        return results;
    }
public:
    GuitarPreampFirstBlock(double f_s, double out_div=1.0) : GuitarPreampBlock(f_s, out_div, 7)
    {
        //
    }
};

class GuitarPreampNextBlock : public GuitarPreampBlock
{
    VectorType f(VectorType &x, double u_in, VectorType &u_cp, VectorType &i_cp)
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
            i_c2p = i_cp(2),
            i_g1 = i_g(u_g1 - u_k1),
            i_a1 = i_a(u_a1 - u_k1, u_g1 - u_k1),
            i_g2 = i_g(u_g2 - u_k2),
            i_a2 = i_a(u_a2 - u_k2, u_g2 - u_k2),
            g_1 = g_2;

        VectorType results(equation_count);
        results =
#include "generated/SecondBlock.cpp.part"
        ;
        return results;
    }
public:
    GuitarPreampNextBlock(double f_s, double out_div=1.0) : GuitarPreampBlock(f_s, out_div, 7)
    {
        //
    }
};

class GuitarPreampThirdBlock : public GuitarPreampNextBlock
{
    const double r_l = 4e6; const double g_l = 1.0f / r_l;
    public:
        GuitarPreampThirdBlock(double f_s, double out_div=1.0) : GuitarPreampNextBlock(f_s, out_div)
        {
            //
        }
};

typedef PipeliningProcessor<Sample, ProcessorElement> ProcessorType;


class GuitarAmp : public ProcessorType
{
    double f_s;
    GuitarPreampFirstBlock *firstBlock;
    GuitarPreampNextBlock *secondBlock;
    GuitarPreampThirdBlock *thirdBlock;
    GuitarPreampNextBlock *fourthBlock;
    CabinetSimulator *cabinet;

public:
    GuitarAmp(double _f_s, const char *impulseResponsePath) : PipeliningProcessor(), f_s(_f_s)
    {
        firstBlock = new GuitarPreampFirstBlock(f_s, 10.0);
        secondBlock = new GuitarPreampNextBlock(f_s, 10.0);
        thirdBlock = new GuitarPreampThirdBlock(f_s, 5.0);
        fourthBlock = new GuitarPreampNextBlock(f_s, 150.0);
        cabinet = new CabinetSimulator(impulseResponsePath);
        elements.push_back(firstBlock);
        elements.push_back(secondBlock);
        elements.push_back(thirdBlock);
        elements.push_back(fourthBlock);
        elements.push_back(cabinet);
        setup();
    }

    double getFS()
    {
        return f_s;
    }
};