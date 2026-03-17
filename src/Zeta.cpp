/**
 * @file Zeta.cpp
 * @brief Hurwitz zeta function implementation via Euler-Maclaurin formula.
 *
 * Adapted from https://github.com/andrewfowlie/thermal_funcs (Fowlie, Andrew)
 * and http://fredrikj.net/math/hurwitz_zeta.pdf (Johansson, Fredrik).
 *
 * Modified to remove the GNU Scientific Library (GSL) dependency: the
 * Pochhammer rising-factorial (gsl_sf_poch) is replaced by an inline
 * pure-C++ computation so that the R package builds without a GSL system
 * requirement.
 */

#include <complex>
#include <cmath>
#include "Zeta.h"
using namespace std;

/* -----------------------------------------------------------------------
 * Bernoulli numbers divided by their factorial  B_{2n} / (2n)!
 * (first 51 entries, index 0 = B_0/0! = 1).
 * --------------------------------------------------------------------- */
static constexpr int B_2n_fact_size = 51;
static constexpr double B_2n_fact[B_2n_fact_size] = {
    1.0,
    0.08333333333333333,
    -0.001388888888888889,
     3.306878306878307e-5,
    -8.267195767195768e-7,
     2.08767569878681e-8,
    -5.284190138687493e-10,
     1.3382536530684679e-11,
    -3.3896802963225827e-13,
     8.586062056277845e-15,
    -2.174868698558062e-16,
     5.50900282836023e-18,
    -1.3954464685812525e-19,
     3.534707039629467e-21,
    -8.953517427037546e-23,
     2.267952452337683e-24,
    -5.744790668872202e-26,
     1.455172475614865e-27,
    -3.6859949406653103e-29,
     9.336734257095045e-31,
    -2.36502241570063e-32,
     5.990671762482135e-34,
    -1.51745488446829e-35,
     3.843758125454189e-37,
    -9.73635307264669e-39,
     2.466247044200681e-40,
    -6.247076741820743e-42,
     1.5824030244644914e-43,
    -4.008273685948936e-45,
     1.0153075855569557e-46,
    -2.5718041582418717e-48,
     6.514456035233815e-50,
    -1.6501309906896525e-51,
     4.179830628539476e-53,
    -1.0587634667702908e-54,
     2.6818791912607708e-56,
    -6.793279351107421e-58,
     1.7207577616681404e-59,
    -4.3587303293488934e-61,
     1.1040792903684668e-62,
    -2.7966655133781345e-64,
     7.084036501679471e-66,
    -1.794407408289224e-67,
     4.545287063611096e-69,
    -1.1513346631982053e-70,
     2.9163647710923614e-72,
    -7.387238263497337e-74,
     1.871209311763795e-75,
    -4.739828557761799e-77,
     1.2006125993354507e-78,
    -3.0411872415142924e-80
};

/**
 * Rising factorial (Pochhammer symbol): (s)_n = s*(s+1)*...*(s+n-1).
 * Implemented without GSL.
 */
static inline double poch(double s, int n)
{
    double result = 1.0;
    for (int i = 0; i < n; ++i)
        result *= (s + static_cast<double>(i));
    return result;
}

/* Partial sum: S(s, a, N) = sum_{k=0}^{N-1} (a+k)^{-s} */
static complex<double> S(double s, complex<double> a, int N)
{
    complex<double> sum = 0.0;
    for (int k = 0; k < N; ++k)
        sum += pow(a + static_cast<complex<double>>(k), -s);
    return sum;
}

/* Integral approximation: I(s, a, N) = (a+N)^{1-s} / (s-1) */
static complex<double> I(double s, complex<double> a, int N)
{
    return pow(a + static_cast<complex<double>>(N), 1.0 - s) / (s - 1.0);
}

/* Asymptotic correction via Euler-Maclaurin with Bernoulli numbers */
static complex<double> T(double s, complex<double> a, int N, int M)
{
    const complex<double> d     = a + static_cast<complex<double>>(N);
    const complex<double> factor = pow(d, -s);

    if (M > B_2n_fact_size)
        M = B_2n_fact_size;

    complex<double> sum = 0.0;
    for (int k = 1; k <= M; ++k)
    {
        /* Pochhammer: (s)_{2k-1} */
        const double rising = poch(s, 2 * k - 1);
        sum += B_2n_fact[k] * rising / pow(d, 2 * k - 1);
    }

    return factor * (0.5 + sum);
}

static complex<double> hurwitz_zeta(double s, complex<double> a, int N)
{
    if (N > B_2n_fact_size)
        N = B_2n_fact_size;
    return S(s, a, N) + I(s, a, N) + T(s, a, N, N);
}

double real_hurwitz_zeta(double s, double a, int N)
{
    return hurwitz_zeta(s, static_cast<complex<double>>(a), N).real();
}
