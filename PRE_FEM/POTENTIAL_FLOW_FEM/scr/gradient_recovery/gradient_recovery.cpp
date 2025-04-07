#include "../include/gradient_recovery/gradient_recovery.h"

void GradientRecovery::initGradientRecovery(const ElemConn &ElemConnData, const CoordArr &Coord, const NodeSol &FieldVal) {
    NumNodeTol = Coord.size();
    this->Coord.assign(Coord.begin(), Coord.end());
    this->FieldVal.assign(FieldVal.begin(), FieldVal.end());
    this->initRecoveryPatchConnManager(NumNodeTol, ElemConnData);
    this->genGradientRecovery();
}

void GradientRecovery::genGradientRecovery() {
    for (Ind IndNode = 0; IndNode < NumNodeTol; ++IndNode) {
        CoordArr LocCoord;
        NodeSol SolutionForPatch;
        calcLongestLengthByPatch();
        genLocalCoordAndSolution(IndNode, LocCoord, SolutionForPatch);
        calcGradientRecoveryFromPatch(IndNode, LocCoord, SolutionForPatch);
    }
}

void GradientRecovery::getCoordinateForPatch(const Ind IndNode, CoordArr &PatchCoord) {
    PatchCoord.push_back(Coord[IndNode]);
    for (const auto &IndConnectedNode : FullNodeRecoveryConn[IndNode]) {
        if (IndConnectedNode != IndNode) {
            PatchCoord.push_back(Coord[IndConnectedNode]);
        }
    }
}

void GradientRecovery::genLocalCoordAndSolution(const Ind IndNode, CoordArr &LocCoord, NodeSol &SolutionForPatch) {
    LocCoord.clear();
    SolutionForPatch.clear();

    double x = Coord[IndNode].first;
    double y = Coord[IndNode].second;

    double MaxLengthOfPatch = LongestLengthForFirstLevelPatch[IndNode];

    for (const auto &IndConnectedNode : FullNodeRecoveryConn[IndNode]) {
        double xi = (Coord[IndConnectedNode].first - x) / MaxLengthOfPatch;
        double eta = (Coord[IndConnectedNode].second - y) / MaxLengthOfPatch;
        LocCoord.emplace_back(xi, eta);
        SolutionForPatch.push_back(FieldVal[IndConnectedNode]);
    }
}

void GradientRecovery::calcGradientRecoveryFromPatch(const Ind IndNode, const CoordArr &LocCoord, const NodeSol &PatchFieldVal) {
    ColVect FittingRHS;
    ColVect FittingCoefficient;
    Mat RawFittingMat, SquaredFittingMat;
    buildRawFittingMat(LocCoord, RawFittingMat);
    buildSquaredFittingMatAndRHS(RawFittingMat, PatchFieldVal, SquaredFittingMat, FittingRHS);
    solveForFittingCoefficient(SquaredFittingMat, FittingRHS, FittingCoefficient);
    buildDerivativeFromFittingCoefficient(IndNode, FittingCoefficient);
}

void GradientRecovery::buildRawFittingMat(const CoordArr &LocCoord, Mat &RawFittingMat) {
    Size OrderOfFittingPolynomial = 6;
    Size NumOfPatchPoint = LocCoord.size();
    allocateMat(RawFittingMat, NumOfPatchPoint, OrderOfFittingPolynomial);
    for (Ind IndNode = 0; IndNode < NumOfPatchPoint; ++IndNode) {
        double xi = LocCoord[IndNode].first;
        double eta = LocCoord[IndNode].second;
        RawFittingMat[IndNode][0] = 1.0;
        RawFittingMat[IndNode][1] = xi;
        RawFittingMat[IndNode][2] = eta;
        RawFittingMat[IndNode][3] = xi * xi;
        RawFittingMat[IndNode][4] = xi * eta;
        RawFittingMat[IndNode][5] = eta * eta;
    }
}

void GradientRecovery::buildSquaredFittingMatAndRHS(const Mat &RawFittingMat, const ColVect &PatchFieldVal, Mat &SquaredFittingMat, ColVect &FittingRHS) {
    Size Row = RawFittingMat.size();
    Size Col = RawFittingMat[0].size();
    allocateMat(SquaredFittingMat, Col, Col);
    allocateColVect(FittingRHS, Col);

    for (std::size_t k = 0; k < Row; ++k) {
        for (std::size_t i = 0; i < Col; ++i) {
            FittingRHS[i] += RawFittingMat[k][i] * PatchFieldVal[k];
            for (std::size_t j = 0; j < Col; ++j) {
                SquaredFittingMat[i][j] += RawFittingMat[k][i] * RawFittingMat[k][j];
            }
        }
    }
}

void GradientRecovery::solveForFittingCoefficient(Mat &SquaredFittingMat, const ColVect &FittingRHS, ColVect &FittingCoefficient) {
    DenseSolver.solveDenseMatrixSystem(SquaredFittingMat, FittingRHS, FittingCoefficient);
}

void GradientRecovery::buildDerivativeFromFittingCoefficient(const Ind IndNode, const ColVect &FittingCoefficient) {
    allocateRecoveredDevArr();
    double LongestLengthForPatch = LongestLengthForFirstLevelPatch[IndNode];
    RecoveredDev[IndNode][0] = FittingCoefficient[1] / LongestLengthForPatch;
    RecoveredDev[IndNode][1] = FittingCoefficient[2] / LongestLengthForPatch;
    RecoveredDev[IndNode][2] = 2.0 * FittingCoefficient[3] / std::pow(LongestLengthForPatch, 2.0);
    RecoveredDev[IndNode][3] = FittingCoefficient[4] / std::pow(LongestLengthForPatch, 2.0);
    RecoveredDev[IndNode][4] = 2.0 * FittingCoefficient[5] / std::pow(LongestLengthForPatch, 2.0);
}

double GradientRecovery::calcLongestLengthByPatch(const Ind IndPatch) {
    double MaxLength = 0.0;
    for (const auto &IndFirstLevelNode : FirstLevelNodeRecoveryConn[IndPatch]) {
        double dx = Coord[IndFirstLevelNode].first - Coord[IndPatch].first;
        double dy = Coord[IndFirstLevelNode].second - Coord[IndPatch].second;
        double TempMaxLength = std::sqrt(dx * dx + dy * dy);
        if (TempMaxLength > MaxLength)
            MaxLength = TempMaxLength;
    }
    return MaxLength;
}

void GradientRecovery::calcLongestLengthByPatch() {
    for (Ind IndPatch = 0; IndPatch < FirstLevelNodeRecoveryConn.size(); ++IndPatch) {
        double LongestLengthByPatch = calcLongestLengthByPatch(IndPatch);
        LongestLengthForFirstLevelPatch.push_back(LongestLengthByPatch);
    }
}

void GradientRecovery::dispPatchCoord(const CoordArr &PatchCoord) const {
    std::cout << "PatchCoord.size()" << PatchCoord.size() << std::endl;
    for (Ind IndNode = 0; IndNode < PatchCoord.size(); ++IndNode) {
        std::cout << "IndNode = " << IndNode << " (" << PatchCoord[IndNode].first << ", " << PatchCoord[IndNode].second << ")" << std::endl;
    }
}

void GradientRecovery::getRecoveryDev(DevArr &RecoveredDev) {
    RecoveredDev.assign(this->RecoveredDev.begin(), this->RecoveredDev.end());
}