#ifndef GRADIENT_RECOVERY_H
#define GRADIENT_RECOVERY_H

#include "recovery_patch_conn_manager.h"
#include "../include/dense_system_solver.h"

#include <cmath>
#include <vector>
#include <utility>
#include <iostream>

class GradientRecovery : protected RecoveryPatchConnManager {
private:
    inline void allocateRecoveredDevArr() { RecoveredDev.resize(NumNodeTol); }

private:
    void genGradientRecovery();
    void getCoordinateForPatch(const Ind IndNode, CoordArr &PatchCoord);
    void buildRawFittingMat(const CoordArr &LocPatch, Mat &RawFittingMat);
    void buildSquaredFittingMatAndRHS(const Mat &RawFittingMat, const ColVect &PatchNodalSolution, Mat &SquaredFittingMat, ColVect &FittingRHS);
    void calcGradientRecoveryFromPatch(const Ind IndNode, const CoordArr &LocCoord, const NodeSol &PatchNodalSolution);
    void genLocalCoordAndSolution(const Ind IndNode, CoordArr &LocCoord, NodeSol &SolutionForPatch);
    double calcLongestLengthByPatch(const Ind IndPatch);
    void calcLongestLengthByPatch();
    void solveForFittingCoefficient(Mat &SquaredFittingMat, const ColVect &FittingRHS, ColVect &FittingCoefficient);
    void buildDerivativeFromFittingCoefficient(const Ind IndNode, const ColVect &FittingCoefficient);

// private:
    // RecoveryPatchConnManager RecoveryPatchConn;

protected:
    // Size NumNodeTol;
    ColVect LongestLengthForFirstLevelPatch;
    // ElemConn ElemConnData;
    // NodeConnSetArr FirstLevelNodeRecoveryConn;
    // NodeConnSetArr FullNodeRecoveryConn;
    // NodeConnArr NodeRecoveryConnArr;

private:
    DenseSystemSolver DenseSolver;

protected:
    CoordArr Coord;
    NodeSol FieldVal;
    DevArr RecoveredDev;

public:
    GradientRecovery() = default;
    ~GradientRecovery() { std::cout << "DESTROY GRADIENT RECOVERY" << std::endl; }
    void initGradientRecovery(const ElemConn &ElemConnData, const CoordArr &Coord, const NodeSol &NodalSolution);

public:
    void dispPatchCoord(const CoordArr &PatchCoord) const;
    void getRecoveryDev(DevArr &RecoveredDev);
};


#endif