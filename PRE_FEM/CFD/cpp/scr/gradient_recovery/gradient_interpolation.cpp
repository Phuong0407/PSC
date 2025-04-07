#include "../include/gradient_recovery/gradient_interpolation.h"

void GradientInterpolation::initGradientInterpolation(const ElemConn &ElemConnData, const CoordArr &Coord, const NodeSol &FieldVal) {
    RecoveryHandler.initGradientRecovery(ElemConnData, Coord, FieldVal);
    RecoveryHandler.getRecoveryDev(this->RecoveredDev);
    this->Coord.assign(Coord.begin(), Coord.end());
    this->FieldVal.assign(FieldVal.begin(), FieldVal.end());
    this->ElemConnData.assign(ElemConnData.begin(), ElemConnData.end());
    calcAreaElem();
}

Ind GradientInterpolation::getIndexOfBelongingElem(double x, double y) {
    for (Ind IndElem = 0; IndElem < ElemConnData.size(); ++IndElem) {
        TriConn TriElem = ElemConnData[IndElem];
        if (isPointInTriangle(x, y, TriElem)) {
            return IndElem;
        } else {
            continue;
        }
    }
    return INVALID_ELEMENT_INDEX;
}

bool GradientInterpolation::isPointInTriangle(double z, double r, const TriConn TriElem) {
    double z1 = Coord[TriElem[0]].first, r1 = Coord[TriElem[0]].second;
    double z2 = Coord[TriElem[1]].first, r2 = Coord[TriElem[1]].second;
    double z3 = Coord[TriElem[2]].first, r3 = Coord[TriElem[2]].second;

    double Denominator = (r2 - r3) * (z1 - z3) + (z3 - z2) * (r1 - r3);

    double a = ((r2 - r3) * (z - z3) + (z3 - z2) * (r - r3)) / Denominator;
    double b = ((r3 - r1) * (z - z3) + (z1 - z3) * (r - r3)) / Denominator;
    double c = 1.0 - a - b;

    if (a >= 0 && a <= 1 && b >= 0 && b <= 1 && c >= 0 && c <= 1) {
        return true;
    } else {
        return false;
    }
}

void GradientInterpolation::interpolateDev(double z, double r, PointDev &InterpolatedDev) {
    Ind IndElem = getIndexOfBelongingElem(z, r);
    Ind IndFirstNode = ElemConnData[IndElem][0]; 
    Ind IndSecondNode = ElemConnData[IndElem][1]; 
    Ind IndThirdNode = ElemConnData[IndElem][2]; 
    
    double z1 = Coord[IndFirstNode].first;
    double r1 = Coord[IndFirstNode].second;
    double z2 = Coord[IndSecondNode].first;
    double r2 = Coord[IndSecondNode].second;
    double z3 = Coord[IndThirdNode].first;
    double r3 = Coord[IndThirdNode].second;

    double Area = AreaPerElem[IndElem];

    double FirstCoefficient = ((z2 * r3 - z3 * r2) + (r2 - r3) * z + (z3 - z2) * r) / (2 * Area);
    double SecondCoefficient = ((z3 * r1 - z1 * r3) + (r3 - r1) * z + (z1 - z3) * r) / (2 * Area);
    double ThirdCoefficient = 1.0 - FirstCoefficient - SecondCoefficient;

    double zDev = FirstCoefficient * RecoveredDev[IndFirstNode][0]
                    + SecondCoefficient * RecoveredDev[IndSecondNode][0]
                    + ThirdCoefficient * RecoveredDev[IndThirdNode][0];
    double rDev = FirstCoefficient * RecoveredDev[IndFirstNode][1]
                    + SecondCoefficient * RecoveredDev[IndSecondNode][1]
                    + ThirdCoefficient * RecoveredDev[IndThirdNode][1];
    double zzDev = FirstCoefficient * RecoveredDev[IndFirstNode][2]
                    + SecondCoefficient * RecoveredDev[IndSecondNode][2]
                    + ThirdCoefficient * RecoveredDev[IndThirdNode][2];
    double zrDev = FirstCoefficient * RecoveredDev[IndFirstNode][3]
                    + SecondCoefficient * RecoveredDev[IndSecondNode][3]
                    + ThirdCoefficient * RecoveredDev[IndThirdNode][3];
    double rrDev = FirstCoefficient * RecoveredDev[IndFirstNode][4]
                    + SecondCoefficient * RecoveredDev[IndSecondNode][4]
                    + ThirdCoefficient * RecoveredDev[IndThirdNode][4];
    InterpolatedDev = {zDev, rDev, zzDev, zrDev, rrDev};
}

void GradientInterpolation::calcAreaElem() {
    for (const auto & Elem : ElemConnData) {
        Ind IndFirstNode = Elem[0]; 
        Ind IndSecondNode = Elem[1];
        Ind IndThirdNode = Elem[2];
        double z1 = Coord[IndFirstNode].first;
        double r1 = Coord[IndFirstNode].second;
        double z2 = Coord[IndSecondNode].first;
        double r2 = Coord[IndSecondNode].second;
        double z3 = Coord[IndThirdNode].first;
        double r3 = Coord[IndThirdNode].second;
        double Area = calcTriangularArea(z1, r1, z2, r2, z3, r3);
        AreaPerElem.push_back(Area);
    }
}