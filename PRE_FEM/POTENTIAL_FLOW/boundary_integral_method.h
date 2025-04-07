#ifndef BOUNDARY_INTEGRAL_METHOD_H
#define BOUNDARY_INTEGRAL_METHOD_H

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include "helper.h"
#include "green_function.h"
#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <array>

DoubleArr calcK1K2(const Point2D &Target, const Point2D &Ref1, const Point2D &Ref2, const NormVec &ElemNormVec);

DoubleArr calcK1K2Singular(const Point2D &Target, const Point2D &Ref1, const Point2D &Ref2, const double SolidAngle);

void buildSystemOfEquation(const Point2DArr &NodeCoord, const ElemArr &ElemConn, const NormVecArr &OutwardNormVec, const DoubleArr &SolidAngle, const DoubleArr &DevPhi, Mat &A, VecCol &B);

inline bool isNodeBelongingToElem(Ind NodeInd, Ind ElemInd, const ElemArr &ElemConn);

#endif
