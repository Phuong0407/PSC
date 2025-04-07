#ifndef HELPER_H
#define HELPER_H

#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <array>

using DoubleArr = std::vector<double>;

struct Point2D {
    double z = 0.0;
    double r = 0.0;
    Point2D(double z_, double r_) : z(z_), r(r_) {}
};

struct NormVec {
    double nz = 0.0;
    double nr = 0.0;
    NormVec(double nz_, double nr_) : nz(nz_), nr(nr_) {}
};
using NormVecArr = std::vector<NormVec>;

using Ind = unsigned int;
using Size = unsigned int;
using IndArr = std::vector<Ind>;

using Weight = double;
using GaussPoint = double;
using GaussQuadraturePoint = std::pair<GaussPoint, Weight>;

using Vec2D = std::pair<double, double>;
using Vec2DArr = std::vector<Vec2D>;
using DoubleArr = std::vector<double>;

using Point2DArr = std::vector<Point2D>;

using GreenFunctionAndDev = std::pair<double, double>;
using ElemGreenFunc = std::array<double, 4>;

using Elem = std::pair<Ind, Ind>;
using ElemArr = std::vector<Elem>;

using VecCol = std::vector<double>;
using Mat = std::vector<VecCol>;

struct GaussQuadrature8Point {
    static constexpr std::array<GaussQuadraturePoint, 8> Point = {{
        { -0.96028985649753623168, 0.10122853629037625915 },
        { -0.79666647741362673959, 0.22238103445337447054 },
        { -0.52553240991632898582, 0.31370664587788728734 },
        { -0.18343464249564980494, 0.36268378337836209114 },
        {  0.18343464249564980494, 0.36268378337836209114 },
        {  0.52553240991632898582, 0.31370664587788728734 },
        {  0.79666647741362673959, 0.22238103445337447054 },
        {  0.96028985649753623168, 0.10122853629037625915 }
    }};
};
const GaussQuadrature8Point Gauss8Point;

struct LogGaussQuadrature8Point {
    static constexpr std::array<GaussQuadraturePoint, 8> Point = {{
        { -0.0133202441608924650122526725243, 0.16441660472800286831472568326 },
        { -0.0797504290138949384098277291424, 0.23752561002330602050134856196 },
        { 0.197871029326188053794476159516, 0.22684198443191912636887804029 },
        { 0.354151399451909419671463603538, 0.17575407900607024498880562120 },
        { 0.529458575234917277706149699996, 0.11292403024675905185500044208 },
        { 0.701814529939099963837152670310, 0.057872210717782072399527967294 },
        { 0.849379320441106676048309202301, 0.020979073742132978043461524115 },
        { 0.953326450056359788767379678514, 0.0036864071040276190133523212764 }
    }};
};
const LogGaussQuadrature8Point LogGauss8Point;

/**
 * calculate the internal angle between two vertices in a polygon oriented counter-clockwise.
 * @param TargetPoint The point at which the internal angle is calculated.
 * @param BackwardPoint The backward vertex in the polygon sequence.
 * @param ForwardPoint The forward vertex in the polygon sequence.
 * @return The internal angle in radians between the vertices.
 * @note The angle is calculated using the `atan2` function, and the result is adjusted by
 *       subtracting it from `M_PI` due to the counter-clockwise orientation of the polygon.
 */

double calcInternalAngleOfTwoVertice(Point2D TargetPoint, Point2D BackwardPoint, Point2D ForwardPoint) {
    double dz1 = TargetPoint.z - BackwardPoint.r;
    double dr1 = TargetPoint.z - BackwardPoint.r;
    double dz2 = ForwardPoint.z - TargetPoint.r;
    double dr2 = ForwardPoint.z - TargetPoint.r;

    double DotProduct = dz1 * dz2 + dr1 * dr2;
    double CrossProduct = dz1 * dr2 - dr1 * dz2;

    double Angle = M_PI - std::atan2(CrossProduct, DotProduct);
    return Angle;
}

#endif