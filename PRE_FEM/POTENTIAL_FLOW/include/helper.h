#ifndef HELPER_H
#define HELPER_H

#include <cmath>
#include <array>
#include <vector>
#include <utility>
#include <cstddef>
#include <unordered_set>

constexpr unsigned int NumberOfNodePerElem = 3;

using Ind = unsigned int;
using Size = unsigned int;
using IndArr = std::vector<Ind>;
using IndSet = std::unordered_set<Ind>;

using Point2D = std::pair<double, double>;
using CoordArr = std::vector<Point2D>;

using ElemConn = std::array<Ind, NumberOfNodePerElem>;
using ElemConnArr = std::vector<ElemConn>;

using ColVec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;
using SolVec = ColVec;

using DirichletBoundArr = std::vector<std::pair<Ind, double>>;

using NodeConnSet = std::unordered_set<Ind>;
using NodeConnSetArr = std::vector<NodeConnSet>;
using NodeConnArr = std::vector<std::vector<Ind>>;

using PointDev = std::array<double, 5>;
using DevArr = std::vector<std::array<double, 5>>;

constexpr Ind INVALID_ELEMENT_INDEX = std::numeric_limits<Ind>::max();

void allocateVecCol(ColVec& RawVect, Size Col);
void allocateMat(Mat& RawMat, Size Row, Size Col);

inline double calcTriangularArea(double x1, double y1, double x2, double y2, double x3, double y3) {
    return 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
}

#endif