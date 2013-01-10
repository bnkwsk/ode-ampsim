#include <cmath>
#include <string>
#include <iostream>

using namespace std;

#include "spline.hpp"

// 12AX7 parameters
const double mi = 100.0f;
const double E_x = 1.4f;
const double K_g1 = 1060.0f;
const double K_p = 600.0f;
const double K_vb = 300.0f;
const double g_cf = 0.00001;
const double g_co = -0.2;

template <typename T> inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class PlateCurrent
{
public:
    PlateCurrent()
    {
        double u_ak;
        for(int u_ak_i = 0; u_ak_i < u_ak_n; ++u_ak_i)
        {
            u_ak = u_ak_min + (double)u_ak_i * (u_ak_max - u_ak_min) / (double)u_ak_n;
            for(int i = 0; i < u_gk_n; ++i)
            {
                const_u_ak_u_gk[u_ak_i][i] = u_gk_min + (double)i * (u_gk_max - u_gk_min) / (double)u_gk_n;
                const_u_ak_i_a[u_ak_i][i] = i_a(u_ak, const_u_ak_u_gk[u_ak_i][i]);
            }
	    const_u_ak_i_a_second_derivatives[u_ak_i] = spline_cubic_set(u_gk_n, const_u_ak_u_gk[u_ak_i], const_u_ak_i_a[u_ak_i], 0, 0, 0, 0);
        }
    }
    double get(double u_ak, double u_gk)
    {
        last_u_ak_left_n = floor((u_ak - u_ak_min) / ((u_ak_max - u_ak_min) / u_ak_n));
        last_left_i_a = spline_cubic_val(u_gk_n, const_u_ak_u_gk[last_u_ak_left_n], const_u_ak_i_a[last_u_ak_left_n], const_u_ak_i_a_second_derivatives[last_u_ak_left_n], u_gk, &last_i_a_derivative, &last_i_a_second_derivative);
        last_right_i_a = spline_cubic_val(u_gk_n, const_u_ak_u_gk[last_u_ak_left_n + 1], const_u_ak_i_a[last_u_ak_left_n + 1], const_u_ak_i_a_second_derivatives[last_u_ak_left_n + 1], u_gk, &last_i_a_derivative, &last_i_a_second_derivative);
        return last_left_i_a;
    }
private:
    const static int u_gk_n = 500;
    const static double u_gk_min = -200.0f;
    const static double u_gk_max = 200.0f;
    const static int u_ak_n = 500;
    const static double u_ak_min = 0.0f;
    const static double u_ak_max = 400.0f;
    double const_u_ak_u_gk[u_ak_n][u_gk_n];
    double const_u_ak_i_a[u_ak_n][u_gk_n];
    double *const_u_ak_i_a_second_derivatives[u_ak_n];
    double last_i_a_derivative;
    double last_i_a_second_derivative;
    int last_u_ak_left_n;
    double last_left_i_a;
    double last_right_i_a;
public:
    static double i_a(double u_ak, double u_gk)
    {
        double E_1 =  u_ak / K_p * log(1.0f + exp(K_p * (1.0f / mi + u_gk / sqrt(K_vb + pow(u_ak, 2)))));
        return pow(E_1, E_x) / K_g1 * (1.0f + sgn(E_1));
    }
};


int main(int argc, char* argv[])
{
    PlateCurrent i_a;
    
    for(int i = 0; i < 10000; ++i)
        i_a.i_a(350, 1.0f + (double)i * 0.01f);
    
    return 0;
}
