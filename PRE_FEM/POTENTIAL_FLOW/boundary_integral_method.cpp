#include "boundary_integral_method.h"

DoubleArr calcK1K2(const Point2D &Target, const Point2D &Ref1, const Point2D &Ref2, const NormVec &ElemNormVec) {
    double FirstK1Integral = 0.0, SecondK1Integral = 0.0, FirstK2Integral = 0.0, SecondK2Integral = 0.0;
    double z1 = Ref1.z, r1 = Ref1.r;
    double z2 = Ref2.z, r2 = Ref2.r;
    double ElemLength = std::sqrt((z1 - z2) * (z1 - z2) + (r1 - r2) * (r1 - r2));

    double Scale = ElemLength / 2.0;

    for (std::size_t i = 0; i < 8; ++i) {
        double xi = Gauss8Point.Point[i].first;
        double Weight = Gauss8Point.Point[i].second;
        double z = z1 * (1 - xi) / 2.0 + z2 * (1 + xi) / 2.0;
        double r = r1 * (1 - xi) / 2.0 + r2 * (1 + xi) / 2.0;
        Point2D Ref{z, r};
        DoubleArr K1K2 = calcGreenAndDev(Target, Ref, ElemNormVec);
        FirstK1Integral += Weight * (1 - xi) / 2.0 * K1K2[1];
        SecondK1Integral += Weight * (1 + xi) / 2.0 * K1K2[1];
        FirstK2Integral += Weight * (1 - xi) / 2.0 * K1K2[0];
        SecondK2Integral += Weight * (1 + xi) / 2.0 * K1K2[0];
    }

    FirstK1Integral *= Scale;
    SecondK1Integral *= Scale;
    FirstK2Integral *= Scale;
    SecondK2Integral *= Scale;

    DoubleArr Result;
    Result.push_back(FirstK1Integral);
    Result.push_back(SecondK1Integral);
    Result.push_back(FirstK2Integral);
    Result.push_back(SecondK2Integral);
    return Result;
}

DoubleArr calcK1K2Singular(const Point2D &Target, const Point2D &Ref1, const Point2D &Ref2, const double SolidAngle) {
    double FirstK1Integral = 0.0, SecondK1Integral = 0.0, FirstK2Integral = 0.0, SecondK2Integral = 0.0;
    double z1 = Ref1.z, r1 = Ref1.r;
    double z2 = Ref2.z, r2 = Ref2.r;
    double ElemLength = std::sqrt((z1 - z2) * (z1 - z2) + (r1 - r2) * (r1 - r2));

    double Scale = ElemLength / 2.0;

    for (std::size_t i = 0; i < 8; ++i) {
        double xi = Gauss8Point.Point[i].first;
        double Weight = Gauss8Point.Point[i].second;
        double z = z1 * (1 - xi) / 2.0 + z2 * (1 + xi) / 2.0;
        double r = r1 * (1 - xi) / 2.0 + r2 * (1 + xi) / 2.0;
        Point2D RegularRef{z, r};
        DoubleArr RegularK1K2 = calcGreenAndDevRegularPart(Target, RegularRef, SolidAngle);
        FirstK1Integral += Weight * (1 - xi) / 2.0 * RegularK1K2[1];
        SecondK1Integral += Weight * (1 + xi) / 2.0 * RegularK1K2[1];
        FirstK2Integral += Weight * (1 - xi) / 2.0 * RegularK1K2[0];
        SecondK2Integral += Weight * (1 + xi) / 2.0 * RegularK1K2[0];

        xi = LogGauss8Point.Point[i].first;
        Weight = Gauss8Point.Point[i].second;
        z = z1 * (1 - xi) / 2.0 + z2 * (1 + xi) / 2.0;
        r = r1 * (1 - xi) / 2.0 + r2 * (1 + xi) / 2.0;
        Point2D LogRef{z, r};
        DoubleArr LogK1K2 = calcGreenAndDevLogPart(Target, LogRef, SolidAngle);
        FirstK1Integral += Weight * (1 - xi) / 2.0 * LogK1K2[1];
        SecondK1Integral += Weight * (1 + xi) / 2.0 * LogK1K2[1];
        FirstK2Integral += Weight * (1 - xi) / 2.0 * LogK1K2[0];
        SecondK2Integral += Weight * (1 + xi) / 2.0 * LogK1K2[0];
    }

    FirstK1Integral *= Scale;
    SecondK1Integral *= Scale;
    FirstK2Integral *= Scale;
    SecondK2Integral *= Scale;

    DoubleArr Result;
    Result.push_back(FirstK1Integral);
    Result.push_back(SecondK1Integral);
    Result.push_back(FirstK2Integral);
    Result.push_back(SecondK2Integral);
    return Result;
}

void buildSystemOfEquation(const Point2DArr &NodeCoord, const ElemArr &ElemConn, const NormVecArr &OutwardNormVec, const DoubleArr &SolidAngle, const DoubleArr &DevPhi, Mat &A, VecCol &B) {
    A.clear();
    B.clear();
    A.resize(NodeCoord.size());
    B.resize(NodeCoord.size());
    for (auto &row : A)
        row.resize(NodeCoord.size());

    for (Ind NodeInd = 0; NodeInd < NodeCoord.size(); ++NodeInd) {
        Point2D TargetPoint = NodeCoord[NodeInd];
        double SolidAngleAtNode = SolidAngle[NodeInd];

        for (Ind ElemInd = 0; ElemInd < ElemConn.size(); ++ElemInd) {
            Ind FirstPointInd = ElemConn[ElemInd].first;
            Ind SecondPointInd = ElemConn[ElemInd].second;
            Point2D FirstPointOfElem = NodeCoord[FirstPointInd];
            Point2D SecondPointOfElem = NodeCoord[SecondPointInd];
            NormVec OutwardNormVecElem = OutwardNormVec[ElemInd];

            bool isSingular = isNodeBelongingToElem(NodeInd, ElemInd, ElemConn);

            DoubleArr K1K2;
            if (isSingular) {
                K1K2 = calcK1K2Singular(TargetPoint, FirstPointOfElem, SecondPointOfElem, SolidAngleAtNode);
            } else {
                K1K2 = calcK1K2(TargetPoint, FirstPointOfElem, SecondPointOfElem, OutwardNormVecElem);
            }

            A[NodeInd][FirstPointInd] += K1K2[1];
            A[NodeInd][SecondPointInd] += K1K2[2];
            A[NodeInd][NodeInd] = SolidAngleAtNode;
            B[NodeInd] += K1K2[3] * DevPhi[FirstPointInd] + K1K2[4] * DevPhi[SecondPointInd];
        }
    }
}

inline bool isNodeBelongingToElem(Ind NodeInd, Ind ElemInd, const ElemArr &ElemConn) {
    const auto& Element = ElemConn[ElemInd];
    return (NodeInd == Element.first || NodeInd == Element.second);
}