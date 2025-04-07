#ifndef DENSE_SYSTEM_SOLVER_H
#define DENSE_SYSTEM_SOLVER_H

#include "helper.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

class DenseSystemSolver {
private:
    static constexpr double TINY = 1.5e-16;
    static void decomposeLU(Mat& A, IndArr& ind, int& parity, bool& isSingular);
    static void substituteBackwardLU(const Mat& A, const IndArr& ind, std::vector<double>& B);

public:
    // static std::vector<double> solveDenseMatrixSystem(Mat& A, const ColVect& B, ColVect &Solution);
    void solveDenseMatrixSystem(Mat& A, const VecCol& B, VecCol &Solution);
};

#endif
