#ifndef INTEGRATE_HPP_INCLUDED
#define INTEGRATE_HPP_INCLUDED

#include <vector>
#include <assert.h>

#include "legendre_weights.h"
#include "legendre_nodes.h"

namespace integrate
{
    ///https://pomax.github.io/bezierinfo/legendre-gauss.html
    ///https://cbeentjes.github.io/files/Ramblings/QuadratureSphere.pdf
    ///http://homepage.divms.uiowa.edu/~atkinson/papers/SphereQuad1982.pdf
    template<typename T, typename U>
    inline
    auto integrate_1d_raw(const T& func, int n, const U& upper, const U& lower)
    {
        std::vector<float> weights = get_legendre_weights(n);
        std::vector<float> nodes = get_legendre_nodes(n);

        using variable_type = decltype(func(0.f));

        variable_type sum = 0;

        for(int j=0; j < n; j++)
        {
            float w = weights[j];
            float xj = nodes[j];

            U value = ((upper - lower)/2.f) * xj + (upper + lower) / 2.f;

            auto func_eval = w * func(value);

            sum = sum + func_eval;
        }

        return ((upper - lower) / 2.f) * sum;
    }

    template<typename T, typename U>
    inline
    auto integrate_1d(const T& func, int n, const U& upper, const U& lower)
    {
        using variable_type = decltype(func(0.f));
        variable_type sum =  0;

        int pieces = 1;
        U step = (upper - lower) / pieces;

        for(int i=0; i < pieces; i++)
        {
            sum += integrate_1d_raw(func, n, (i + 1) * step + lower, i * step + lower);
        }

        return sum;
    }

    template<typename T>
    inline
    auto spherical_integrate(const T& f_theta_phi, int n, float radius = 1)
    {
        float iupper = 2 * M_PI;
        float ilower = 0;

        float jupper = M_PI;
        float jlower = 0;

        ///https://cbeentjes.github.io/files/Ramblings/QuadratureSphere.pdf7 7
        ///0 -> 2pi, phi
        auto outer_integral = [&](float phi)
        {
            auto inner_integral = [&](float theta){return sin(theta) * f_theta_phi(theta, phi);};

            return integrate_1d(inner_integral, n, jupper, jlower);
        };

        ///expand sphere area over unit sphere
        return integrate_1d(outer_integral, n, iupper, ilower) * radius * radius;
    }

    ///cartesian domain is [-radius, +radius]
    template<typename T>
    inline
    auto cartesian_integrate(const T& f_pos, int n, float cartesian_radius = 1.f, float sphere_radius = 1.f)
    {
        auto cartesian_function = [&](float theta, float phi)
        {
            vec3f pos = {cartesian_radius * cos(phi) * sin(theta), cartesian_radius * sin(phi) * sin(theta), cartesian_radius * cos(theta)};

            return f_pos(pos);
        };

        return spherical_integrate(cartesian_function, n, sphere_radius);
    }
}

#endif // INTEGRATE_HPP_INCLUDED
