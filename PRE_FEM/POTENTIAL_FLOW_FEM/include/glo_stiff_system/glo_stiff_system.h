#ifndef GLO_STIFF_MAT_H
#define GLO_STIFF_MAT_H

#include "../include/common_type_and_method.h"
#include "../include/glo_stiff_system/loc_axisym_stiff_mat.h"
#include "../include/glo_stiff_system/dirichlet_bound.h"
#include "../include/glo_stiff_system/mod_glo_stiff_system.h"

#include <array>
#include <cmath>
#include <vector>
#include <utility>
#include <iostream>

class GloStiffSystem :  protected LocAxisymStiffMat,
                        protected DirichletBound,
                        protected ModificationGloStiffSystem
{
private:
    double L1, L2;
    double InVel, OutVel;
    unsigned int NumNode, NumNodeHor, NumNodeVer;
    CoordArr Coord;
    ElemConn ElemConnData;
    IndArr IndDirichletBound;

private:
    Mat GloStiffMat;
    ColVect ValDirichletBound;
    ColVect ForceMat;

private:
    void getElemCoord(const Ind IndElem, CoordArr &ElemNodeCoord);

private:
    void assembleGloStiffMat(unsigned int IndElem, Mat &LocStiffMat);

private:
    void calcGloStiffMat();

    void storeParameter(
        double L1, double L2,
        double InVel, double OutVel,
        unsigned int NumNode, unsigned int NumNodeHor, unsigned int NumNodeVer,
        const CoordArr &Coord, const ElemConn &ElemConnData, IndArr &IndDirichletBound
    );
    void initializeMatrice();
    void setupDirichletBoundary();
    void applyDirichletBoundaryCondition();

public:
    GloStiffSystem() = default;
    void initGloStiffSystem(
        double L1, double L2,
        double InVel, double OutVel,
        unsigned int NumNode, unsigned int NumNodeHor, unsigned int NumNodeVer,
        const CoordArr &Coord, const ElemConn &ElemConnData, IndArr &IndDirichletBound
    );

public:
    void getGloStiffSystem(Mat &GloStiffMat, ColVect &ValDirichletBound, ColVect &ForceMat);
};

#endif