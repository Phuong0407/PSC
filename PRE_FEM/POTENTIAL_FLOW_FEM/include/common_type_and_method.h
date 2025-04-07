#ifndef COMMON_TYPE_AND_METHOD_H
#define COMMON_TYPE_AND_METHOD_H

#include <cmath>
#include <array>
#include <vector>
#include <utility>
#include <cstddef>
#include <unordered_set>

using Ind = unsigned int;
using Size = unsigned int;
using NodeSol = std::vector<double>;
using IndArr = std::vector<unsigned int>;
using TriConn = std::array<unsigned int, 3>;
using ElemConn = std::vector<std::array<unsigned int, 3>>;
using Point2D = std::pair<double, double>;
using ColVect = std::vector<double>;
using Mat = std::vector<std::vector<double>>;

using CoordArr = std::vector<std::pair<double, double>>;
using PointDev = std::array<double, 5>;
using DevArr = std::vector<std::array<double, 5>>;

using NodeConnSet = std::unordered_set<unsigned int>;
using NodeConnArr = std::vector<std::vector<unsigned int>>;
using NodeConnSetArr = std::vector<std::unordered_set<unsigned int>>;

constexpr int NumNodePerElem = 3;
constexpr Ind INVALID_ELEMENT_INDEX = std::numeric_limits<Ind>::max();

void allocateMat(Mat& RawMat, Size Row, Size Col);
void allocateColVect(ColVect& RawVect, Size Col);

inline double calcTriangularArea(double z1, double r1, double z2, double r2, double z3, double r3) {
    return 0.5 * std::abs(z1 * (r2 - r3) + z2 * (r3 - r1) + z3 * (r1 - r2));
}


#endif