#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "helper.h"
#include <cmath>
#include <vector>
#include <array>

struct EllipticIntegralCoefficient {
    const std::array<double, 5> a = {{
        1.38629436112,
        0.09666344259,
        0.03459092383,
        0.03742563713,
        0.01451196212
    }};
    
    const std::array<double, 5> b = {{
        0.50000000000,
        0.12498593597,
        0.06880248576,
        0.03328355346,
        0.00441787012
    }};
};
EllipticIntegralCoefficient EllipticCoefficient;

inline double calcRegularK(double m) {
    double m1 = 1 - m;
    double FirstTerm = 0.0, SecondTerm = 0.0;
    for (int i = EllipticCoefficient.a.size() - 1; i >= 0; --i)
        FirstTerm = FirstTerm * m1 + EllipticCoefficient.a[i];

    for (int i = EllipticCoefficient.b.size() - 1; i > 0; --i)
        SecondTerm = SecondTerm * m1 + EllipticCoefficient.b[i];
    return FirstTerm + SecondTerm * std::log(1/m1);
}

inline double calcLogK(double m) {
    double m1 = 1 - m;
    return EllipticCoefficient.b[0] * std::log(1/m1);
}

inline double calcGammaCoefficient(double nr, double a, double b) {
    return 2.0 * nr / std::sqrt(a + b);
}

inline double calcDeltaCofficient(double nr, double nz, double r, double r_target, double zd, double a, double b) {
    return r * (nr * (r - r_target) + nz * zd) / ((a - b) * std::sqrt(a + b));
}

void calcK1K2Integrand(double z_target, double r_target, double z, double r, double nz, double nr, double &DevGreen, double &Green) {
    double zd = z - z_target;
    double a = r_target * r_target + r * r + zd * zd, b = 2.0 * r * r_target;
    double m = 4.0 * r * r_target / ((r + r_target) * (r + r_target) + zd * zd);
    double gamma = calcGammaCoefficient(nr, a, b);
    double delta = calcDeltaCofficient(nr, nz, r, r_target, zd, a, b);
    double Km = std::comp_ellint_1(std::sqrt(m));
    DevGreen = (gamma + delta) * std::comp_ellint_2(std::sqrt(m)) - gamma * Km;
    Green = 4.0 * r / std::sqrt(a + b) * Km;
}

void calcRegularK1K2Integrand(double z_target, double r_target, double z, double r, double alpha, double &DevGreen_reg, double &Green_reg) {
    double zd = z - z_target;
    double a = r_target * r_target + r * r + zd * zd, b = 2.0 * r * r_target;
    double m = 4.0 * r * r_target / ((r + r_target) * (r + r_target) + zd * zd);
    double Km_reg = std::comp_ellint_1(std::sqrt(m)) - EllipticCoefficient.b[0] * std::log(1/(1.0 - m));
    DevGreen_reg = - 2.0 * std::sin(alpha) / std::sqrt(a + b) * (std::comp_ellint_2(std::sqrt(m)) - Km_reg);
    Green_reg = 4.0 * r / std::sqrt(a + b) * Km_reg;
}

void calcLogK1K2Integrand(double z_target, double r_target, double z, double r, double alpha, double &DevGreen_log, double &Green_log) {
    double zd = z - z_target;
    double a = r_target * r_target + r * r + zd * zd, b = 2.0 * r * r_target;
    double m = 4.0 * r * r_target / ((r + r_target) * (r + r_target) + zd * zd);
    double Km_log = EllipticCoefficient.b[0] * std::log(1/(1.0 - m));
    DevGreen_log = 2.0 * std::sin(alpha) / std::sqrt(a + b) * Km_log;
    Green_log = 4.0 * r / std::sqrt(a + b) * Km_log;
}

#endif