#include "../include/laplace_solver.h"

LaplaceSolver::LaplaceSolver(
                                double L1, double L2,
                                double H1, double H2, double H3,
                                unsigned int N, unsigned int M1, unsigned int M2, unsigned int M3,
                                double InVel
                            )
                            : L1(L1), L2(L2),
                            H1(H1), H2(H2), H3(H3),
                            N(N), M1(M1), M2(M2), M3(M3),
                            InVel(InVel)
{
    GridDat.initGridGeneration(L1, L2, H1, H2, H3, N, M1, M2, M3);
    GridDat.getGridDat(Coord, ElemConnData, IndDirichletBound);
    calcOutputVel();

    unsigned int NumNode = GridDat.getTolNumNode();
    unsigned int NumNodeHor = GridDat.getNumNodeHor();
    unsigned int NumNodeVer = GridDat.getNumNodeVer();

    StiffMat.initGloStiffSystem(L1, L2, InVel, OutVel, NumNode, NumNodeHor, NumNodeVer, Coord, ElemConnData, IndDirichletBound);
    StiffMat.getGloStiffSystem(GloStiffMat, ValDirichletBound, ForceMat);
}

void LaplaceSolver::solveForNodalVal() {
    Solver.initSparseSolver(GloStiffMat);
    Solver.solveSparseSystem(ForceMat, Solution);
}

void LaplaceSolver::getNodalSolution(NodeSol &Solution){
    Solution = this->Solution;
}

void LaplaceSolver::calcOutputVel() {
    OutVel = InVel / (1 - L2 * L2 / (L1 * L1));
}

void LaplaceSolver::getElemConn(ElemConn &ElemConn) const {
    ElemConn = this->ElemConnData;
}

void LaplaceSolver::getCoord(CoordArr &Coord) const {
    Coord = this->Coord;
}
