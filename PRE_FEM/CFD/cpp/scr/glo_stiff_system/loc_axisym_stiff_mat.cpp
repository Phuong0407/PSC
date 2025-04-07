#include "../include/glo_stiff_system/loc_axisym_stiff_mat.h"

void LocAxisymStiffMat::calcLocAxisymStiffMat(Point2D FirstPoint, Point2D SecondPoint, Point2D ThirdPoint, Mat &AxisymStiffMat) {
    std::vector<double> c, d;
    allocateColVect(c, NumNodePerElem);
    allocateColVect(d, NumNodePerElem);
    allocateMat(AxisymStiffMat, NumNodePerElem, NumNodePerElem);

    double z1 = FirstPoint.first;
    double r1 = FirstPoint.second;
    double z2 = SecondPoint.first;
    double r2 = SecondPoint.second;
    double z3 = ThirdPoint.first;
    double r3 = ThirdPoint.second;

    c[0] = z2 - z3;
    c[1] = z3 - z1;
    c[2] = z1 - z2;
    
    d[0] = r3 - r1;
    d[1] = r1 - r3;
    d[2] = r2 - r1;

    double Area = calcTriangularArea(z1, r1, z2, r2, z3, r3);
    double RCentroid = (r1 + r2 + r3) / 3.0;

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            AxisymStiffMat[i][j] = M_PI * RCentroid * (c[i] * c[j] + d[i] * d[j]) / (2 * Area);
        }
    }
}