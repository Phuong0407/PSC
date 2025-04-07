#ifndef GRID_GENERATION_H
#define GRID_GENERATION_H

#include "helper.h"
#include <iostream>

class GridGeneration {
private:
    double L1, L2;
    double H1, H2, H3;
    void initVerticalPhysicalDimension(double L1, double L2);
    void initHorizontalPhysicalDimension(double H1, double H2, double H3);

private:
    Size NumElemFirstDom, NumElemSecondDom, NumElemThirdDom;
    Size NumElemVer, NumElemHor;
    Size NumNodeVer, NumNodeHor;
    Size NumElemTol, NumNodeTol;
    void initHorizontalGridDimension(Size NumElemFirstDom, Size NumElemSecondDom, Size NumElemThirdDom);
    void initVerticalGridDimension(Size NumElemVer);
    void initTotalNumberOfNode();
    void initTotalNumberOfElem();

private:
    double InVel, OutVel;

private:
    CoordArr NodeCoord;
    void buildNodeCoord();

private:
    ElemConnArr GridConn;
    void buildGridConnection();

private:
    DirichletBoundArr DirichletBoundCondition;
    void buildAxisBoundaryCondition();
    void buildSpinnerBoundaryCondition();
    void buildExitBoundaryCondition();
    void buildWallBoundaryCondition();
    void buildInletBoundaryCondition();

private:
    IndArr AxisNode;
    void buildIndexOfDirichletBoundaryNode();

private:
    IndArr InletNode;
    IndArr ExitNode;
    IndArr WallNode;
    IndArr SpinnerNode;

    void buildIndexOfNeumanBoundaryNode();
    void buildIndexOfSpinnerNeumannBoundaryNode();
    void buildIndexOfInletAndExitNode();
    void buildIndexOfWallNeumannNode();

public:
    GridGeneration() =  default;
    ~GridGeneration() {
        std::cout << "DESTRUCTION OF GRID GENERATION OBJECT!!!" << std::endl;
    }

    void initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3);
    void initGridDomainDimension(Size NumElemVer, Size NumElemFirstDom, Size NumElemSecondDom, Size NumElemThirdDom);
    void buildGridGenerationData();
    void buildBoundaryVelocity(double InVel);
    void buildBoundaryCondition();
};



#endif