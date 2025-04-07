#ifndef CSR_SPARSE_SOLVER_H
#define CSR_SPARSE_SOLVER_H

#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <limits>

class CSRSparseSolver {
private:
    double ZeroComparedTol;
    double MatNorm;
    long unsigned int SizeMat;
    long unsigned int NonZeroElem;

    std::vector<double> Val;
    std::vector<long unsigned int> ColInd;
    std::vector<long unsigned int> RowPtr;

    const double CONV_TOL = 1.0E-6;
    const int MAX_ITER = 100000;

public:
    CSRSparseSolver() = default;
    void initSparseSolver(const std::vector<std::vector<double>>& SPARSE_MAT);
    void multSparseMat(const std::vector<double>& X, std::vector<double>& Y) const;
    void solveSparseSystem(const std::vector<double>& RHS_VECT, std::vector<double>& X) const;
};

#endif