#include "../include/grid_generation.h"

void GridGeneration::initVerticalPhysicalDimension(double L1, double L2) {
    this->L1 = L1;
    this->L2 = L2;
}
void GridGeneration::initHorizontalPhysicalDimension(double H1, double H2, double H3) {
    this->H1 = H1;
    this->H2 = H2;
    this->H3 = H3;
}

void GridGeneration::initHorizontalGridDimension(Size NumElemFirstDom, Size NumElemSecondDom, Size NumElemThirdDom) {
    this->NumElemFirstDom = NumElemFirstDom;
    this->NumElemSecondDom = NumElemSecondDom;
    this->NumElemThirdDom = NumElemThirdDom;
    this->NumElemHor = this->NumElemFirstDom
                         + this->NumElemSecondDom
                         + this->NumElemThirdDom;
    this->NumNodeHor = this->NumElemHor + 1;
}

void GridGeneration::initVerticalGridDimension(Size NumElemVer) {
    this->NumElemVer = NumElemVer;
    this->NumNodeVer = this->NumElemVer + 1;
}

void GridGeneration::initTotalNumberOfNode() {
    this->NumNodeTol = this->NumNodeVer * this->NumNodeVer;
}

void GridGeneration::initTotalNumberOfElem() {
    this->NumElemTol = this->NumElemHor * this->NumElemVer;
}

void GridGeneration::initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3) {
    initVerticalPhysicalDimension(L1, L2);
    initHorizontalPhysicalDimension(H1, H2, H3);
}

void GridGeneration::initGridDomainDimension(Size NumElemVer, Size NumElemFirstDom, Size NumElemSecondDom, Size NumElemThirdDom) {
    initVerticalGridDimension(NumElemVer);
    initHorizontalGridDimension(NumElemFirstDom, NumElemSecondDom, NumElemThirdDom);
    initTotalNumberOfElem();
    initTotalNumberOfNode();
}

void GridGeneration::buildNodeCoord(){
    const double dx1 = H1 / NumElemFirstDom;
    const double dx2 = H2 / NumElemSecondDom;
    const double dx3 = H3 / NumElemThirdDom;
    const double dy1 = L1 / NumElemVer;
    const double dy2 = L2 / NumElemVer;

    for (unsigned int i = 0; i <= NumElemVer; i++) {
        for (unsigned j = 0; j <= NumElemFirstDom; j++) {
            double x = dx1 * j;
            double y = dy1 * i;
            NodeCoord.push_back(std::make_pair(x, y));
        }

        double a = ((dy2 - dy1) * i + L1 - L2) / H1;
        double b = dy1 * i - a * H1;
        for (unsigned int j = NumElemFirstDom + 1; j <= NumElemFirstDom + NumElemSecondDom; ++j) {
            double x = H1 + dx2 * (j - NumElemFirstDom);
            double y = a * x + b;
            NodeCoord.push_back(std::make_pair(x, y));
        }

        for (unsigned j = NumElemFirstDom + NumElemSecondDom + 1; j <= NumElemFirstDom + NumElemSecondDom + NumElemThirdDom; ++j) {
            double x = H1 + H2 + dx3 * (j - NumElemFirstDom - NumElemSecondDom);
            double y = L1 - L2 + dy2 * i;
            NodeCoord.push_back(std::make_pair(x, y));
        }
    }
}

void GridGeneration::buildGridConnection() {
    for (unsigned int i = 0; i < NumElemVer; ++i) {
        for (unsigned int j = 0; j < NumElemHor; ++j) {
            unsigned int IndElem = i * NumNodeHor + j;
            unsigned int FirstNodeInd = IndElem;
            unsigned int SecondNodeInd = FirstNodeInd + 1;
            unsigned int ThirdNodeInd = SecondNodeInd + NumNodeHor;
            unsigned int FourthNodeInd = FirstNodeInd + NumNodeHor;
            GridConn.push_back({FirstNodeInd, SecondNodeInd, ThirdNodeInd});
            GridConn.push_back({FirstNodeInd, ThirdNodeInd, FourthNodeInd});
        }
    }
}

void GridGeneration::buildBoundaryVelocity(double InVel) {
    this->InVel = InVel;
    this->OutVel = InVel / (1 - std::pow(L2/L1, 2.0));
}

void GridGeneration::buildAxisBoundaryCondition() {
    for (Ind IndNode = 0; IndNode <= NumElemFirstDom; ++IndNode) {
        double Phi = 
        DirichletBoundCondition.emplace_back(IndNode, Phi);
    }
}


// void GridGeneration::buildIndexOfDirichletBoundaryNode() {
//     for (std::size_t IndNode = 0; IndNode < NumElemFirstDom; ++IndNode) {
//         AxisNode.push_back(IndNode);
//     }
// }

// void GridGeneration::buildIndexOfNeumanBoundaryNode() {
//     buildIndexOfSpinnerNeumannBoundaryNode();
//     buildIndexOfInletAndExitNode();
//     buildIndexOfWallNeumannNode();
// }

// void GridGeneration::buildAxisBoundaryCondition() { 

// }

// void GridGeneration::buildInletBoundaryCondition() {
//     for (std::size_t i = 1; i < NumElemVer - 1; ++i) {
//         Ind InletElemInd = NumNodeHor * i;
//         Ind OutletElemInd = NumNodeHor * (i + 1) - 1;
//         InletNode.push_back(InletElemInd);
//         ExitNode.push_back(OutletElemInd);
//     }
// }

// void GridGeneration::buildWallBoundaryCondition() {
//     for (std::size_t i = 0; i < NumNodeHor; ++i) {
//         Ind WallInd = NumElemVer * NumNodeHor + i;
//         WallNode.push_back(WallInd);
//     }
// }

// void GridGeneration::buildSpinnerBoundaryCondition() {
//     for (std::size_t IndNode = NumElemFirstDom + 1; IndNode < NumElemFirstDom + NumNodeVer; ++IndNode) {
//         SpinnerNode.push_back(IndNode);

//     }
// }

void GridGeneration::buildGridGenerationData() {
    buildNodeCoord();
    buildGridConnection();
    buildIndexOfDirichletBoundaryNode();
    buildIndexOfNeumanBoundaryNode();
}