#include "../include/helper.h"

void PhysicalDomainDimension::initVerticalDimension(double L1, double L2) {
    this->L1 = L1;
    this->L2 = L2;
}

void PhysicalDomainDimension::initHorizontalDimension(double H1, double H2, double H3) {
    this->H1 = H1;
    this->H2 = H2;
    this->H3 = H3;
}

void PhysicalDomainDimension::initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3) {
    initVerticalDimension(L1, L2);
    initHorizontalDimension(H1, H2, H3);
}
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
void GridDimension::initNumberOfVerticalElements(unsigned int N) {
    NumElemVer = N;
}

void GridDimension::initNumberOfHorizontalElements(unsigned int M1, unsigned int M2, unsigned int M3) {
    NumElemFirstDom = M1;
    NumElemSecondDom = M2;
    NumElemThirdDom = M3;
    NumElemHor = NumElemFirstDom + NumElemSecondDom + NumElemThirdDom;
}

void GridDimension::initNumberOfElements() {
    NumElemTol = 2 * NumElemVer * NumElemHor;
}

void GridDimension::initNumberOfNodes() {
    NumNodeHor = NumElemHor + 1;
    NumNodeVer = NumElemVer + 1;
    NumNodeTol = NumNodeHor * NumNodeVer;
    NumElemTol = NumElemVer * NumElemHor;
}

void GridDimension::initGridDimension(unsigned int N, unsigned int M1, unsigned int M2, unsigned int M3) {
    initNumberOfVerticalElements(N);
    initNumberOfHorizontalElements(M1, M2, M3);
    initNumberOfElements();
    initNumberOfNodes();
}
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
void GridGeneration::initGridGeneration(double L1, double L2, double H1, double H2, double H3, unsigned int N, unsigned int M1, unsigned int M2, unsigned int M3, double InVel) {
    initPhysicalDomainDimension(L1, L2, H1, H2, H3);
    initGridDimension(N, M1, M2, M3);
    initGridStep();
    genGridCoord();
    genElemConn();
    InFlux = - InVel;
    OutFlux = 
    // genIndDirichletBound();
}

void GridGeneration::initGridStep() {
    dx1 = H1 / NumElemFirstDom;
    dx2 = H2 / NumElemSecondDom;
    dx3 = H3 / NumElemThirdDom;
    dy1 = L1 / NumElemVer;
    dy2 = L2 / NumElemVer;
}

void GridGeneration::genGridCoord(){
    for (unsigned int i = 0; i <= NumElemVer; i++) {
        for (unsigned j = 0; j <= NumElemFirstDom; j++) {
            double x = dx1 * j;
            double y = dy1 * i;
            Coord.emplace_back(x, y);
        }

        double a = ((dy2 - dy1) * i + L1 - L2) / H1;
        double b = dy1 * i - a * H1;
        for (unsigned int j = NumElemFirstDom + 1; j <= NumElemFirstDom + NumElemSecondDom; ++j) {
            double x = H1 + dx2 * (j - NumElemFirstDom);
            double y = a * x + b;
            Coord.emplace_back(x, y);
        }

        for (unsigned j = NumElemFirstDom + NumElemSecondDom + 1; j <= NumElemFirstDom + NumElemSecondDom + NumElemThirdDom; ++j) {
            double x = H1 + H2 + dx3 * (j - NumElemFirstDom - NumElemSecondDom);
            double y = L1 - L2 + dy2 * i;
            Coord.emplace_back(x, y);
        }
    }
}

void GridGeneration::genElemConn() {
    ElemConnData.resize(NumElemTol, ElemConn(3));
    for (unsigned int i = 0; i < NumElemVer; ++i) {
        for (unsigned int j = 0; j < NumElemHor; ++j) {
            unsigned int IndElem = i * NumNodeHor + j;
            unsigned int FirstNodeInd = IndElem;
            unsigned int SecondNodeInd = FirstNodeInd + 1;
            unsigned int ThirdNodeInd = SecondNodeInd + NumNodeHor;
            unsigned int FourthNodeInd = FirstNodeInd + NumNodeHor;
            ElemConn LowerTriangleElem(3);
            LowerTriangleElem[0] = FirstNodeInd;
            LowerTriangleElem[1] = SecondNodeInd;
            LowerTriangleElem[2] = ThirdNodeInd;
            ElemConnData.push_back(LowerTriangleElem);
            ElemConn UpperTriangleElem(3);
            UpperTriangleElem[0] = FirstNodeInd;
            UpperTriangleElem[1] = ThirdNodeInd;
            UpperTriangleElem[2] = FourthNodeInd;
            ElemConnData.push_back(UpperTriangleElem);
            // push_back({FirstNodeInd, SecondNodeInd, ThirdNodeInd});
            // ElemConnData.push_back({FirstNodeInd, ThirdNodeInd, FourthNodeInd});
        }
    }
}

void GridGeneration::genNeumannBoundary() {

}
// void GridGeneration::genIndDirichletBound() {
//     std::set<unsigned int> OrderedIndSet;

//     for (unsigned int i = 0; i < NumNodeHor; ++i) {
//         unsigned int LowerIndBound = i;
//         unsigned int UpperIndBound = (NumNodeVer - 1) * NumNodeHor + i;
//         OrderedIndSet.insert(LowerIndBound);
//         OrderedIndSet.insert(UpperIndBound);
//     }

//     for (unsigned int i = 1; i < NumNodeVer - 1; ++i) {
//         unsigned int LeftIndBound = i * NumNodeHor;
//         unsigned int RightIndBound = (i + 1) * NumNodeHor - 1;
//         OrderedIndSet.insert(LeftIndBound);
//         OrderedIndSet.insert(RightIndBound);
//     }

//     // IndDirichletBound.assign(OrderedIndSet.begin(), OrderedIndSet.end());
// }



void GridGeneration::dispGridCoord() const {
    for (std::size_t i = 0; i < Coord.size(); ++i) {
        std::cout << "(" << Coord[i].z << ", " << Coord[i].r << ")" << std::endl;
    }
}

void GridGeneration::dispElemConn() const {
    for (std::size_t i = 0; i < ElemConnData.size(); ++i) {
        std::cout << "[" << ElemConnData[i][0] << ", " << ElemConnData[i][1] << ", " << ElemConnData[i][2] << "]" << std::endl;
    }
}

void GridGeneration::dispIndDirichletBound() const {
    std::cout << "THE INDEX OF DIRICHLET BOUNDARY NODE" << std::endl;
    // for (std::size_t i = 0; i < IndDirichletBound.size(); ++i) {
    //     std::cout << "[" << IndDirichletBound[i] << "]" << std::endl;
    // }
}

void GridGeneration::dispNeumannBound() const {
    std::cout << "THE NEUMANN BOUNDARY ELEMENT INDEX - FIRSTNODEIND - SECONDNODEIND - FLUX1 - FLUX2" << std::endl;
}


void GridGeneration::getGridDat(Point2DArr &Coord, ElemConnArr &ElemConnData, EdgeFluxArr &NeumannBound) {
    Coord = this->Coord;
    ElemConnData = this->ElemConnData;
    NeumannBound = this->NeumannBound;
}
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
