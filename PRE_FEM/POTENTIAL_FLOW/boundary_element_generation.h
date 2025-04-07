#ifndef BOUNDARY_ELEMENT_GENERATION_H
#define BOUNDARY_ELEMENT_GENERATION_H

#include "helper.h"

class BoundaryElementGeneration {
private:
    double L1, L2;
    double H1, H2, H3;
    void initVerticalPhysicalDimension(double L1, double L2);
    void initHorizontalPhysicalDimension(double H1, double H2, double H3);

private:
    Size NumElemFirstBound;
    Size NumElemSecondBound;
    Size NumElemThirdBound;
    Size NumElemFourthBound;
    Size NumElemFifthBound;
    Size NumElemSixthBound;
    Size NumElemTot;
    Size NumNodeTot;

    Point2DArr NodeCoord;
    NormVecArr OutwardNormVec;
    DoubleArr SolidAngle;
    DoubleArr Flux;
    ElemArr ElemConn;

private:
    void initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3);
    void initBoundaryElementDimension(Size NumElemFirstBound, Size NumElemSecondBound, Size NumElemThirdBound, Size NumElemFourBound, Size NumElemFifthBound, Size NumElemSixthBound);

private:
    void buildElemConn();
    void buildNodeCoord();
    void buildOutwardNormVec();
    void buildSolidAngle();

public:
    BoundaryElementGeneration() =  default;
    ~BoundaryElementGeneration() { std::cout << "DESTRUCTION OF BOUNDARY ELEMENT GENERATION OBJECT!!!" << std::endl; }

public:
    void initBoundaryElementGeneration(double L1, double L2, double H1, double H2, double H3, Size NumElemFirstBound, Size NumElemSecondBound, Size NumElemThirdBound, Size NumElemFourBound, Size NumElemFifthBound, Size NumElemSixthBound);
    void initVelocity(double InVec);
public:
    void dispElemConn();
    void dispNodeCoord();
    void dispOutwardNormVec();
    void dispSolidAngle();
    void dispFlux();

public:
    void getNodeCoord(Point2DArr &NodeCoord);
    void getOutwardNormVec(NormVecArr &OutwardNormVec);
    void getSolidAngle(DoubleArr &SolidAngle);
    void getElemConn(ElemArr &ElemConn);
    void getFlux(DoubleArr &Flux);

};

#endif