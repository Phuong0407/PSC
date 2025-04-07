#include "../include/stiffness_system.h"

void StiffnessSystem::initGlobalStiffnessSystem(double InletVel, const ElemConnArr &GridConn, double H, double L1, double L2, const CoordArr &NodeCoord, const IndArr &AxisNode, const IndArr &SpinnerNode, const IndArr &ExitNode, const IndArr &WallNode, const IndArr &InletNode) {
    buildGlobalStiffnessMatrix(GridConn, NodeCoord);
    double ExitVel = InletVel / (1 - std::pow(L2 / L1, 2.0));
    buildBoundaryCondition(InletVel, ExitVel, H, NodeCoord, AxisNode, SpinnerNode, ExitNode, WallNode, InletNode);

}

bool StiffnessSystem::isElemNeumannInletBoundary(Ind IndElem, const ElemConnArr &GridConn, const IndSet &InletElemIndNeumann) {
    Ind IndFirstNode = GridConn[IndElem][0];
    Ind IndSecondNode = GridConn[IndElem][3];
    return (InletElemIndNeumann.find(IndFirstNode) != InletElemIndNeumann.end()) && (InletElemIndNeumann.find(IndSecondNode) != InletElemIndNeumann.end());
}

bool StiffnessSystem::isElemNeumannExitBoundary(Ind IndElem, const ElemConnArr &GridConn, const IndSet &ExitElemIndNeumann) {
    Ind IndFirstNode = GridConn[IndElem][0];
    Ind IndSecondNode = GridConn[IndElem][3];
    return (ExitElemIndNeumann.find(IndFirstNode) != ExitElemIndNeumann.end()) && (ExitElemIndNeumann.find(IndSecondNode) != ExitElemIndNeumann.end());
}

void StiffnessSystem::buildGlobalStiffnessMatrix(const ElemConnArr &GridConn, const CoordArr &NodeCoord) {
    Size NumberOfNode = NodeCoord.size();
    allocateMat(GloStiffnessMatrix, NumberOfNode, NumberOfNode);
    for (Ind IndElem = 0; IndElem < GridConn.size(); ++ IndElem) {
        Ind IndFirstNode = GridConn[IndElem][0];
        Ind IndSecondNode = GridConn[IndElem][1];
        Ind IndThirdNode = GridConn[IndElem][3];
        double z1 = NodeCoord[IndFirstNode].first, r1 = NodeCoord[IndFirstNode].second;
        double z2 = NodeCoord[IndSecondNode].first, r2 = NodeCoord[IndSecondNode].second;
        double z3 = NodeCoord[IndThirdNode].first, r3 = NodeCoord[IndThirdNode].second;
        Mat LocalStiffnessMatrix;
        calcLocalAxisymmetricElem(z1, r1, z2, r2, z3, r3, LocalStiffnessMatrix);
        for (std::size_t LocRow; LocRow < NumberOfNodePerElem; ++LocRow) {
            unsigned int GloRow = GridConn[IndElem][LocRow];
            for (unsigned int LocCol = 0; LocCol < 3; ++LocCol) {
                Ind GloCol = GridConn[IndElem][LocCol];
                GloStiffnessMatrix[GloRow][GloCol] += LocalStiffnessMatrix[LocRow][LocCol];
            }
        }
    }
}

void StiffnessSystem::buildAxisBoundaryCondition(double InVel, const IndArr &IndDirichlet, const CoordArr &NodeCoord) {
    for (std::size_t IndNode = 0; IndNode < IndDirichlet.size(); ++IndNode) {
        Ind IndexOfDirichletBoundary = IndDirichlet[IndNode];
        double z = NodeCoord[IndexOfDirichletBoundary].first;
        double Phi = InVel * z;
        ValAxis.push_back(Phi);
    }
}

void StiffnessSystem::buildSpinnerBoundaryCondition(const IndArr &SpinnerIndNeumann) {
    for (std::size_t IndNode = 0; IndNode < SpinnerIndNeumann.size(); ++IndNode) {
        ValSpinner.push_back(0.0);
    }
}

void StiffnessSystem::buildWallBoundaryCondition(const IndArr &WallIndNewmann) {
    for (std::size_t IndNode = 0; IndNode < WallIndNewmann.size(); ++IndNode) {
        ValWall.push_back(0.0);
    }
}

void StiffnessSystem::buildBoundaryCondition(double InVel, double ExitVel, double H, const CoordArr &NodeCoord, const IndArr &AxisNode, const IndArr &SpinnerNode, const IndArr &ExitNode, const IndArr &WallNode, const IndArr &InletNode) {
    buildAxisBoundaryCondition(InVel, AxisNode, NodeCoord);
    buildSpinnerBoundaryCondition(SpinnerNode);
    buildExitBoundaryCondition(ExitVel, ExitNode, H);
    buildWallBoundaryCondition(WallNode);
    buildInletBoundaryCondition(InVel, InletNode, 0);
}

void StiffnessSystem::calcLocalAxisymmetricElem(double z1, double r1, double z2, double r2, double z3, double r3, Mat &LocalStiffnessMatrix) {
    allocateMat(LocalStiffnessMatrix, NumberOfNodePerElem, NumberOfNodePerElem);

    ColVec c(3), d(3);
    c[0] = z2 - z3, c[1] = z3 - z1, c[2] = z1 - z2;
    d[0] = r3 - r1, d[1] = r1 - r3, d[2] = r2 - r1;
    
    double Area = calcTriangularArea(z1, r1, z2, r2, z3, r3);
    double RCentroid = (r1 + r2 + r3) / 3.0;

    for(std::size_t i = 0; i < NumberOfNodePerElem; ++i) {
        for(std::size_t j = 0; j < NumberOfNodePerElem; ++j) {
            LocalStiffnessMatrix[i][j] = M_PI * RCentroid * (c[i] * c[j] + d[i] * d[j]) / (2 * Area);
        }
    }
}

void StiffnessSystem::buildInletBoundaryCondition(double InVel, const IndArr &InletNode, double z) {
    for (std::size_t i = 0; i < InletNode.size(); ++i) {
        Ind IndNode = InletNode[i];
        double Phi = - InVel * z;
        ValInlet.push_back(Phi);
    }
}

void StiffnessSystem::buildExitBoundaryCondition(double ExitVel, const IndArr &ExitNode, double z) {
    for (std::size_t i = 0; i < ExitNode.size(); ++i) {
        Ind IndNode = ExitNode[i];
        double Phi = ExitVel * z;
        ValInlet.push_back(Phi);
    }
}

void StiffnessSystem::modifyGlobalStiffnessSystem() {
    for (std::size_t i = 0; i < AxisNode.size(); ++i) {
        Ind IndNode = AxisNode[i];

    }
}

void StiffnessSystem::modifyGlobalStiffnessSystemUsingDirichletBoundary(const IndArr &IndDirichletBound, const ColVec &ValDirichletBoundArr) {
    for (std::size_t i = 0; i < IndDirichletBound.size(); ++i) {
        Ind ModifyingCol = IndDirichletBound[i];
        double ValDirichletBoundNode = ValDirichletBoundArr[i];
        for (std::size_t row = 0; row < ForceMat.size(); ++row) {
            ForceMat[row] -= GloStiffnessMatrix[row][ModifyingCol] * ValDirichletBoundNode;
        }
        for (std::size_t row = 0; row < GloStiffnessMatrix.size(); ++row) {
            GloStiffnessMatrix[row][ModifyingCol] = 0.0;
            GloStiffnessMatrix[ModifyingCol][row] = 0.0;
        }
        GloStiffnessMatrix[ModifyingCol][ModifyingCol] = 1.0;
    }
}

