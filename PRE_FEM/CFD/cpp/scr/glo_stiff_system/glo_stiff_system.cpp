#include "../include/glo_stiff_system/glo_stiff_system.h"

void GloStiffSystem::initGloStiffSystem(
    double L1, double L2,
    double InVel, double OutVel,
    unsigned int NumNode, unsigned int NumNodeHor, unsigned int NumNodeVer,
    const CoordArr &Coord, const ElemConn &ElemConnData, IndArr &IndDirichletBound
) {
    storeParameter(L1, L2, InVel, OutVel, NumNode, NumNodeHor, NumNodeVer, Coord, ElemConnData, IndDirichletBound);
    initializeMatrice();
    calcGloStiffMat();
    setupDirichletBoundary();
    applyDirichletBoundaryCondition();
}

void GloStiffSystem::storeParameter(
    double L1, double L2,
    double InVel, double OutVel,
    unsigned int NumNode, unsigned int NumNodeHor, unsigned int NumNodeVer,
    const CoordArr &Coord, const ElemConn &ElemConnData, IndArr &IndDirichletBound
) {
    this->L1 = L1;
    this->L2 = L2;
    this->InVel = InVel;
    this->OutVel = OutVel;
    this->NumNode = NumNode;
    this->NumNodeHor = NumNodeHor;
    this->NumNodeVer = NumNodeVer;
    this->Coord = Coord;
    this->ElemConnData = ElemConnData;
    this->IndDirichletBound = IndDirichletBound;
}

void GloStiffSystem::initializeMatrice() {
    allocateMat(GloStiffMat, NumNode, NumNode);
    allocateColVect(ForceMat, NumNode);
}

void GloStiffSystem::setupDirichletBoundary() {
    genDirichletBound(L1, L2, InVel, OutVel, NumNodeHor, NumNodeVer, Coord, IndDirichletBound, ValDirichletBound);
}

void GloStiffSystem::applyDirichletBoundaryCondition() {
    modifyGloStiffSystem(IndDirichletBound, ValDirichletBound, GloStiffMat, ForceMat);
}

void GloStiffSystem::calcGloStiffMat() {
    Mat LocStiffMat;
    for (std::size_t IndElem = 0; IndElem < ElemConnData.size(); ++IndElem) {
        CoordArr ElemNodeCoord(NumNodePerElem);
        getElemCoord(IndElem, ElemNodeCoord);
        
        Point2D FirstPoint = ElemNodeCoord[0];
        Point2D SecondPoint = ElemNodeCoord[1];
        Point2D ThirdPoint = ElemNodeCoord[2];

        calcLocAxisymStiffMat(FirstPoint, SecondPoint, ThirdPoint, LocStiffMat);
        assembleGloStiffMat(IndElem, LocStiffMat);
    }
}

void GloStiffSystem::assembleGloStiffMat(unsigned int IndElem, Mat &LocStiffMat) {
    for (unsigned int LocRow = 0; LocRow < 3; ++LocRow) {
        unsigned int GloRow = ElemConnData[IndElem][LocRow];
        for (unsigned int LocCol = 0; LocCol < 3; ++LocCol) {
            unsigned int GloCol = ElemConnData[IndElem][LocCol];
            GloStiffMat[GloRow][GloCol] += LocStiffMat[LocRow][LocCol];
        }
    }
}

void GloStiffSystem::getGloStiffSystem(Mat &GloStiffMat, ColVect &ValDirichletBound, ColVect &ForceMat) {
    GloStiffMat = this->GloStiffMat;
    ValDirichletBound = this->ValDirichletBound;
    ForceMat = this->ForceMat;
}

void GloStiffSystem::getElemCoord(const Ind IndElem, CoordArr &ElemNodeCoord) {
    IndArr IndNode(NumNodePerElem);
    CoordArr ElemNode(NumNodePerElem);
    for (std::size_t i = 0; i < NumNodePerElem; ++i) {
        IndNode[i] = ElemConnData[IndElem][i];
        ElemNodeCoord[i] = {Coord[IndNode[i]].first, Coord[IndNode[i]].second};
    }
}
