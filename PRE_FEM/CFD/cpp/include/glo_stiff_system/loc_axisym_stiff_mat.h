#ifndef LOC_AXISYM_STIFFNESS_MAT_H
#define LOC_AXISYM_STIFFNESS_MAT_H

#include "../include/common_type_and_method.h"

#include <cmath>
#include <vector>

class LocAxisymStiffMat {
protected:
    LocAxisymStiffMat() = default;
    void calcLocAxisymStiffMat(Point2D FirstPoint, Point2D SecondPoint, Point2D ThirdPoint, Mat &AxisymStiffMat);
};



#endif