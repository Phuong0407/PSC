#include "../include/helper.h"

void allocateVecCol(ColVec& RawVect, Size Col) {
    if (RawVect.size() == Col) {
        return;
    }
    RawVect.assign(Col, 0.0);
    return;
}

void allocateMat(Mat& RawMat, Size Row, Size Col) {
    if (RawMat.size() == Row) {
        bool isAllRowCorrectSize = true;
        for (const auto& r : RawMat) {
            if (r.size() != Col) {
                isAllRowCorrectSize = false;
                break;
            }
        }
        if (isAllRowCorrectSize) {
            return;
        }
    }
    RawMat.assign(Row, std::vector<double>(Col, 0.0));
    return;
}