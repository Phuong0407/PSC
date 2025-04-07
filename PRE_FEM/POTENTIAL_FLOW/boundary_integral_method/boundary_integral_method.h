#ifndef BOUNDARY_INTEGRAL_METHOD_H
#define BOUNDARY_INTEGRAL_METHOD_H

#include "green_function.h"

inline bool isNodeBelongingToElem(Ind NodeInd, Ind ElemInd, const ElemArr &ElemConn);
void calcK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double nz, double nr, double &K11, double &K12, double &K21, double &K22);
void calcRegularK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double alpha, double &K11, double &K12, double &K21, double &K22);
void calcLogK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double alpha, double &K11, double &K12, double &K21, double &K22);
void calcSingularK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double alpha, double &K11, double &K12, double &K21, double &K22);
void buildSystemOfEquation(const Point2DArr &NodeCoord, const ElemArr &ElemConn, const NormVecArr &OutwardNormVec, const DoubleArr &SolidAngle, const DoubleArr &DevPhi, Mat &A, VecCol &B);

/**
 * Computes the integrals of K1 and K2 over a specified region using 8 points Gaussian quadrature.
 * @param z_target The reference value of the z-coordinate for the integrand calculation.
 * @param r_target The reference value of the r-coordinate for the integrand calculation.
 * @param z1 The starting value of the z-coordinate for the integration range.
 * @param r1 The starting value of the r-coordinate for the integration range.
 * @param z2 The ending value of the z-coordinate for the integration range.
 * @param r2 The ending value of the r-coordinate for the integration range.
 * @param nz The z-component of the outward normal vector of the element.
 * @param nr The r-component of the outward normal vector of the element.
 * @param K11 The result of the integral for the first K1 component. This value is modified in the function.
 * @param K12 The result of the integral for the second K1 component. This value is modified in the function.
 * @param K21 The result of the integral for the first K2 component. This value is modified in the function.
 * @param K22 The result of the integral for the second K2 component. This value is modified in the function.
 */
void calcK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double nz, double nr, double &K11, double &K12, double &K21, double &K22) {
    K11 = 0.0, K12 = 0.0, K21 = 0.0, K22 = 0.0;
    double ElemLength = std::sqrt((z1 - z2) * (z1 - z2) + (r1 - r2) * (r1 - r2));
    double Scale = ElemLength / 2.0;

    for (std::size_t i = 0; i < 8; ++i) {
        double xi = Gauss8Point.Point[i].first;
        double z = z1 * (1 - xi) / 2.0 + z2 * (1 + xi) / 2.0;
        double r = r1 * (1 - xi) / 2.0 + r2 * (1 + xi) / 2.0;
        double weight = Gauss8Point.Point[i].second;
        double DevGreen = 0.0, GreenFunc = 0.0;
        calcK1K2Integrand(z_target, r_target, z, r, nz, nr, DevGreen, GreenFunc);
        K11 += (1 - xi) / 2.0 * DevGreen * weight;
        K12 += (1 + xi) / 2.0 * DevGreen * weight;
        K21 += (1 - xi) / 2.0 * GreenFunc * weight;
        K22 += (1 + xi) / 2.0 * GreenFunc * weight;
    }
    K11 *= Scale;
    K12 *= Scale;
    K21 *= Scale;
    K22 *= Scale;
}

/**
 * Computes the integrals of the regular terms K1 and K2 over a specified region
 * using 8 points Gaussian quadrature.
 * @param z_target The reference value of the z-coordinate for the integrand calculation.
 * @param r_target The reference value of the r-coordinate for the integrand calculation.
 * @param z1 The starting value of the z-coordinate for the integration range.
 * @param r1 The starting value of the r-coordinate for the integration range.
 * @param z2 The ending value of the z-coordinate for the integration range.
 * @param r2 The ending value of the r-coordinate for the integration range.
 * @param alpha A constant that affects the integrand.
 * @param K11 The result of the integral for the first K1 component. This value is modified in the function.
 * @param K12 The result of the integral for the second K1 component. This value is modified in the function.
 * @param K21 The result of the integral for the first K2 component. This value is modified in the function.
 * @param K22 The result of the integral for the second K2 component. This value is modified in the function.
 */
void calcRegularK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double alpha, double &K11, double &K12, double &K21, double &K22) {
    K11 = 0.0, K12 = 0.0, K21 = 0.0, K22 = 0.0;
    double ElemLength = std::sqrt((z1 - z2) * (z1 - z2) + (r1 - r2) * (r1 - r2));
    double Scale = ElemLength / 2.0;

    for (std::size_t i = 0; i < 8; ++i) {
        double xi = Gauss8Point.Point[i].first;
        double z = z1 * (1 - xi) / 2.0 + z2 * (1 + xi) / 2.0;
        double r = r1 * (1 - xi) / 2.0 + r2 * (1 + xi) / 2.0;
        double weight = Gauss8Point.Point[i].second;
        double DevGreen = 0.0, GreenFunc = 0.0;
        calcRegularK1K2Integrand(z_target, r_target, z, r, alpha, DevGreen, GreenFunc);
        K11 += (1 - xi) / 2.0 * DevGreen * weight;
        K12 += (1 + xi) / 2.0 * DevGreen * weight;
        K21 += (1 - xi) / 2.0 * GreenFunc * weight;
        K22 += (1 + xi) / 2.0 * GreenFunc * weight;
    }
    K11 *= Scale;
    K12 *= Scale;
    K21 *= Scale;
    K22 *= Scale;
}

/**
 * Computes the integrals of the logarithmic terms K1 and K2 over a specified region
 * using 8 points Gaussian quadrature.
 * @param z_target The reference value of the z-coordinate for the integrand calculation.
 * @param r_target The reference value of the r-coordinate for the integrand calculation.
 * @param z1 The starting value of the z-coordinate for the integration range.
 * @param r1 The starting value of the r-coordinate for the integration range.
 * @param z2 The ending value of the z-coordinate for the integration range.
 * @param r2 The ending value of the r-coordinate for the integration range.
 * @param alpha A constant that affects the integrand.
 * @param K11 The result of the integral for the first K1 component. This value is modified in the function.
 * @param K12 The result of the integral for the second K1 component. This value is modified in the function.
 * @param K21 The result of the integral for the first K2 component. This value is modified in the function.
 * @param K22 The result of the integral for the second K2 component. This value is modified in the function.
 */
void calcLogK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double alpha, double &K11, double &K12, double &K21, double &K22) {
    K11 = 0.0, K12 = 0.0, K21 = 0.0, K22 = 0.0;
    double ElemLength = std::sqrt((z1 - z2) * (z1 - z2) + (r1 - r2) * (r1 - r2));
    double Scale = ElemLength;

    for (std::size_t i = 0; i < 8; ++i) {
        double xi = (1.0 + Gauss8Point.Point[i].first) / 2.0;
        double z = z1 * (1 - xi) + z2 * xi;
        double r = r1 * (1 - xi) + r2 * xi;
        double weight = Gauss8Point.Point[i].second / 2.0;
        double DevGreenLog = 0.0, GreenFuncLog = 0.0;
        calcLogK1K2Integrand(z_target, r_target, z, r, alpha, DevGreenLog, GreenFuncLog);
        K11 += (1 - xi) * DevGreenLog * weight;
        K12 += xi * DevGreenLog * weight;
        K21 += (1 - xi) * GreenFuncLog * weight;
        K22 += xi * GreenFuncLog * weight;
    }
    K11 *= Scale;
    K12 *= Scale;
    K21 *= Scale;
    K22 *= Scale;
}

/**
 * Computes the integrals of the singular terms K1 and K2 over a specified region
 * using 8 points Gaussian quadrature.
 *
 * @param z_target The reference value of the z-coordinate for the integrand calculation.
 * @param r_target The reference value of the r-coordinate for the integrand calculation.
 * @param z1 The starting value of the z-coordinate for the integration range.
 * @param r1 The starting value of the r-coordinate for the integration range.
 * @param z2 The ending value of the z-coordinate for the integration range.
 * @param r2 The ending value of the r-coordinate for the integration range.
 * @param alpha A constant that affects the integrand.
 * @param K11 The result of the integral for the first K1 component. This value is modified in the function.
 * @param K12 The result of the integral for the second K1 component. This value is modified in the function.
 * @param K21 The result of the integral for the first K2 component. This value is modified in the function.
 * @param K22 The result of the integral for the second K2 component. This value is modified in the function.
 */
void calcSingularK1K2Integral(double z_target, double r_target, double z1, double r1, double z2, double r2, double alpha, double &K11, double &K12, double &K21, double &K22) {
    double K11_reg = 0.0, K12_reg = 0.0, K21_reg = 0.0, K22_reg = 0.0;
    double K11_log = 0.0, K12_log = 0.0, K21_log = 0.0, K22_log = 0.0;
    calcRegularK1K2Integral(z_target, r_target, z1, r1, z2, r2, alpha, K11_reg, K12_reg, K21_reg, K22_reg);
    calcLogK1K2Integral(z_target, r_target, z1, r1, z2, r2, alpha, K11_log, K12_log, K21_log, K22_log);
    K11 = K11_reg + K11_log, K12 = K12_reg + K12_log, K21 = K21_reg + K21_log, K22 = K22_reg + K22_log;
}

/**
 * Builds a system of linear equations for the boundary element method (BEM)
 * to solve the potential flow around a 2D body.
 *
 * @param NodeCoord The coordinates of the nodes in the BEM mesh.
 * @param ElemConn The connectivity of the elements in the BEM mesh.
 * @param OutwardNormVec The outward normal vectors of the elements.
 * @param SolidAngle The solid angles of the nodes.
 * @param DevPhi The deviations of the potential at the nodes.
 * @param A The system matrix (modified in the function).
 * @param B The right-hand side vector (modified in the function).
 */
// void buildSystemOfEquation(const Point2DArr &NodeCoord, const ElemArr &ElemConn, const NormVecArr &OutwardNormVec, const DoubleArr &SolidAngle, const DoubleArr &DevPhi, Mat &A, VecCol &B) {
//     A.assign(NodeCoord.size(), VecCol(NodeCoord.size(), 0.0));
//     B.assign(NodeCoord.size(), 0.0);

//     for (Ind NodeInd = 0; NodeInd < NodeCoord.size(); ++NodeInd) {
//         double z_target = NodeCoord[NodeInd].z, r_target = NodeCoord[NodeInd].r;
//         double alpha = SolidAngle[NodeInd];
//         for (Ind ElemInd = 0; ElemInd < ElemConn.size(); ++ElemInd) {
//             const auto &[FirstPointInd, SecondPointInd] = ElemConn[ElemInd];
//             double z1 = NodeCoord[FirstPointInd].z, r1 = NodeCoord[FirstPointInd].r;
//             double z2 = NodeCoord[SecondPointInd].z, r2 = NodeCoord[SecondPointInd].r;

//             double K11 = 0.0, K12 = 0.0, K21 = 0.0, K22 = 0.0;

//             if (isNodeBelongingToElem(NodeInd, ElemInd, ElemConn)) {
//                 calcSingularK1K2Integral(z_target, r_target, z1, r1, z2, r2, alpha, K11, K12, K21, K22);
//             } else {
//                 const auto &[nz, nr] = OutwardNormVec[ElemInd];
//                 calcK1K2Integral(z_target, r_target, z1, r1, z2, r2, nz, nr, K11, K12, K21, K22);
//             }
//             A[NodeInd][FirstPointInd] += K11;
//             A[NodeInd][SecondPointInd] += K12;
//             B[NodeInd] += K21 * DevPhi[FirstPointInd] + K22 * DevPhi[SecondPointInd];
//         }
//         A[NodeInd][NodeInd] = alpha;
//     }
// }

void buildSystemOfEquation(const Point2DArr &NodeCoord, const ElemArr &ElemConn, const NormVecArr &OutwardNormVec, const DoubleArr &SolidAngle, const FluxArr &DevPhi, Mat &A, VecCol &B) {
    std::cout << NodeCoord.size() << std::endl;
    A.assign(NodeCoord.size(), VecCol(NodeCoord.size(), 0.0));
    B.assign(NodeCoord.size(), 0.0);

    for (Ind NodeInd = 0; NodeInd < NodeCoord.size(); ++NodeInd) {
        double z_target = NodeCoord[NodeInd].z, r_target = NodeCoord[NodeInd].r;
        double alpha = SolidAngle[NodeInd];
        for (Ind ElemInd = 0; ElemInd < ElemConn.size(); ++ElemInd) {
            const auto &[FirstPointInd, SecondPointInd] = ElemConn[ElemInd];
            double z1 = NodeCoord[FirstPointInd].z, r1 = NodeCoord[FirstPointInd].r;
            double z2 = NodeCoord[SecondPointInd].z, r2 = NodeCoord[SecondPointInd].r;

            double K11 = 0.0, K12 = 0.0, K21 = 0.0, K22 = 0.0;

            if (isNodeBelongingToElem(NodeInd, ElemInd, ElemConn)) {
                calcSingularK1K2Integral(z_target, r_target, z1, r1, z2, r2, alpha, K11, K12, K21, K22);
            } else {
                const auto &[nz, nr] = OutwardNormVec[ElemInd];
                calcK1K2Integral(z_target, r_target, z1, r1, z2, r2, nz, nr, K11, K12, K21, K22);
            }
            A[NodeInd][FirstPointInd] += K11;
            A[NodeInd][SecondPointInd] += K12;
            B[NodeInd] += K21 * DevPhi[ElemInd][0] + K22 * DevPhi[ElemInd][1];
        }
        A[NodeInd][NodeInd] = alpha;
    }
}

/**
 * Checks if a given node belongs to a specified element in the boundary element method (BEM) mesh.
 *
 * @param NodeInd The index of the node to check.
 * @param ElemInd The index of the element to check.
 * @param ElemConn The connectivity of the elements in the BEM mesh.
 *
 * @return True if the node belongs to the element, false otherwise.
 */
inline bool isNodeBelongingToElem(Ind NodeInd, Ind ElemInd, const ElemArr &ElemConn) {
    const auto& Element = ElemConn[ElemInd];
    return (NodeInd == Element.first || NodeInd == Element.second);
}

#endif