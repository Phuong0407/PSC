#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "helper.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <utility>
#include <array>
#include <stdexcept>

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

inline double calcRegularG(double m) {
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

inline DoubleArr calcGreenAndDev(const Point2D &Target, const Point2D &Ref, const NormVec &ElemNormVec) {
    double z = Ref.z, r = Ref.r;
    double ri = Target.r, zi = Target.r;
    double zd = z - zi, b = 2 * r * ri, a = ri * ri + r * r + zd * zd;
    double m = 4 * r * ri / ((r + ri) * (r + ri) + zd * zd);
    double nz = ElemNormVec.nz, nr = ElemNormVec.nr;

    if (m < 0 || m >= 1)
        throw std::domain_error("Parameter m outside the scope (0, 1)! Func calcGreenAndDev!");

    double GreenFunc = 4.0 * std::comp_ellint_1(m) / std::sqrt(a+b);

    double FirstTerm = 2 * nr / (std::sqrt(a + b));
    double SecondTerm = r * (nr * (r - ri) + nz * zd) / ((a - b) * std::sqrt(a + b));
    double GreenDev = (FirstTerm + SecondTerm) * std::comp_ellint_2(m) - FirstTerm * std::comp_ellint_1(m);

    DoubleArr Result;
    Result.push_back(GreenFunc);
    Result.push_back(GreenDev);
    return Result;
}

inline DoubleArr calcGreenAndDevRegularPart(const Point2D &Target, const Point2D &Ref, double Alpha) {
    double z = Ref.z, r = Ref.r;
    double ri = Target.r, zi = Target.r;
    double zd = z - zi, b = 2 * r * ri, a = ri * ri + r * r + zd * zd;
    double m = 4 * r * ri / ((r + ri) * (r + ri) + zd * zd);

    double RegularG = calcRegularG(m);
    double GreenFuncRegular = 4.0 * RegularG / std::sqrt(a+b) * r;
    double GreenDevRegular = - 2 * std::sin(Alpha) / (std::sqrt(a + b)) * (std::comp_ellint_2(m) - RegularG);

    DoubleArr Result;
    Result.push_back(GreenFuncRegular);
    Result.push_back(GreenDevRegular);
    return Result;
}

inline DoubleArr calcGreenAndDevLogPart(const Point2D &Target, const Point2D &Ref, double Alpha) {
    double z = Ref.z, r = Ref.r;
    double ri = Target.r, zi = Target.r;
    double zd = z - zi, b = 2 * r * ri, a = ri * ri + r * r + zd * zd;
    double m = 4 * r * ri / ((r + ri) * (r + ri) + zd * zd);

    double LogG = calcLogK(m);
    double GreenFuncLog = 4.0 * LogG / std::sqrt(a+b) * r;
    double GreenDevLog = 2 * std::sin(Alpha) / (std::sqrt(a + b)) * LogG;

    DoubleArr Result;
    Result.push_back(GreenFuncLog);
    Result.push_back(GreenDevLog);
    return Result;
}

#endif