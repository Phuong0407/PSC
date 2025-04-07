#ifndef STIFFNESS_SYSTEM_H
#define STIFFNESS_SYSTEM_H

#include "../include/helper.h"
#include <iostream>
#include <cmath>

class StiffnessSystem {
private:
    void buildAxisBoundaryCondition(double InVel, const IndArr &AxisInd, const CoordArr &NodeCoord);
    void buildSpinnerBoundaryCondition(const IndArr &SpinnerIndNeumann);
    void buildWallBoundaryCondition(const IndArr &WallIndNewmann);
    void buildInletBoundaryCondition(double InVel, const IndArr &InletNode, double z);
    void buildExitBoundaryCondition(double OutVel, const IndArr &ExitNode, double z);
    void buildInletBoundaryCondition(double InVel, const IndArr &InletNode, double z);
    void buildExitBoundaryCondition(double OutVel, const IndArr &ExitNode, double z);
    void buildBoundaryCondition(double InVel, double ExitVel, double H, const CoordArr &NodeCoord, const IndArr &AxisNode, const IndArr &SpinnerNode, const IndArr &ExitNode, const IndArr &WallNode, const IndArr &InletNode);

private:
    bool isElemNeumannInletBoundary(Ind IndElem, const ElemConnArr &GridConn, const IndSet &InletElemIndNeumann);
    bool isElemNeumannExitBoundary(Ind IndElem, const ElemConnArr &GridConn, const IndSet &ExitElemIndNeumann);

private:
    Mat GloStiffnessMatrix;
    ColVec ForceMat;
    ColVec ValAxis;
    ColVec ValExit;
    ColVec ValInlet;
    ColVec ValWall;
    ColVec ValSpinner;

private:
    IndArr AxisNode;
    IndArr SpinnerNode;
    IndArr ExitNode;
    IndArr WallNode;
    IndArr InletNode;

private:
    void calcLocalAxisymmetricElem(double z1, double r1, double z2, double r2, double z3, double r3, Mat &LocalStiffnessMatrix);
    void buildGlobalStiffnessMatrix(const ElemConnArr &GridConn, const CoordArr &NodeCoord);

public:
    StiffnessSystem() = default;
    ~StiffnessSystem() {
        std::cout << "DESTRUCTION OF STIFFNESS SYSTEM OBJECT!!!" << std::endl;
    }
    void initGlobalStiffnessSystem(double InletVel, const ElemConnArr &GridConn, double H, double L1, double L2, const CoordArr &NodeCoord, const IndArr &AxisNode, const IndArr &SpinnerNode, const IndArr &ExitNode, const IndArr &WallNode, const IndArr &InletNode);
    void modifyGlobalStiffnessSystemUsingDirichletBoundary(const IndArr &IndDirichletBound, const ColVec &ValDirichletBound);
    void modifyGlobalStiffnessSystem();
};

#endif