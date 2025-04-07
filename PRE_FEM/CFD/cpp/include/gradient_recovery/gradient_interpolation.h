#ifndef GRADIENT_INTERPOLATION_H
#define GRADIENT_INTERPOLATION_H

#include "../include/common_type_and_method.h"
#include "../include/gradient_recovery/recovery_patch_conn_manager.h"
#include "../include/gradient_recovery/gradient_recovery.h"

class GradientInterpolation{
private:
    GradientRecovery RecoveryHandler;

private:
    Ind getIndexOfBelongingElem(double x, double y);
    bool isPointInTriangle(double x, double y, const TriConn TriElem);
    void calcAreaElem();

private:
    ElemConn ElemConnData;
    CoordArr Coord;
    NodeSol FieldVal;
    DevArr RecoveredDev;

protected:
    ColVect AreaPerElem;

public:
    GradientInterpolation() = default;
    ~GradientInterpolation() { std::cout << "DESTROY GRADIENT INTERPOLATION" << std::endl; }

public:
    void initGradientInterpolation(const ElemConn &ElemConnData, const CoordArr &Coord, const NodeSol &NodalSolution);
    void interpolateDev(double x, double y, PointDev &InterpolatedDev);
};

#endif