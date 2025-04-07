#ifndef LAPLACE_SOLVER_H
#define LAPLACE_SOLVER_H

#include "common_type_and_method.h"
#include "grid_generation/grid_generation.h"
#include "glo_stiff_system/glo_stiff_system.h"
#include "csr_sparse_solver.h"

#include <cmath>
#include <vector>
#include <iostream>

class LaplaceSolver {
private:
    double L1, L2;
    double H1, H2, H3;
    unsigned int N, M1, M2, M3;
    double InVel, OutVel;

    NodeSol Solution;

    GridGeneration GridDat;
    GloStiffSystem StiffMat;
    CSRSparseSolver Solver;

private:
    CoordArr Coord;
    ElemConn ElemConnData;
    IndArr IndDirichletBound;
    Mat GloStiffMat;
    NodeSol ValDirichletBound;
    ColVect ForceMat;

private:
    void calcOutputVel();

public:
    LaplaceSolver() = default;
    LaplaceSolver(
        double L1, double L2,
        double H1, double H2, double H3,
        unsigned int N, unsigned int M1, unsigned int M2, unsigned int M3,
        double InVel
    );
    void solveForNodalVal();
    void getNodalSolution(NodeSol &Solution);
    void getElemConn(ElemConn &ElemConnData) const;
    void getCoord(CoordArr &Coord) const;
};

#endif