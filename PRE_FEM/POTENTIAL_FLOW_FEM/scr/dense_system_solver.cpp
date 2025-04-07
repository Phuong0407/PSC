#include "../include/dense_system_solver.h"

void DenseSystemSolver::decomposeLU(Mat& A, IndArr& ind, int& parity, bool& isSingular) {
    std::size_t matSize = A.size();
    std::vector<double> implicitScalingFactor(matSize);
    parity = 1;
    isSingular = false;

    for (std::size_t i = 0; i < matSize; ++i) {
        double rowMax = 0.0;
        for (std::size_t j = 0; j < matSize; ++j) {
            rowMax = std::max(rowMax, std::abs(A[i][j]));
        }
        if (rowMax < TINY) {
            isSingular = true;
            return;
        }
        implicitScalingFactor[i] = 1.0 / rowMax;
    }

    for (std::size_t j = 0; j < matSize; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            double sum = A[i][j];
            for (std::size_t k = 0; k < i; ++k) {
                sum -= A[i][k] * A[k][j];
            }
            A[i][j] = sum;
        }

        double pivot = 0.0;
        std::size_t iMax = j;

        for (std::size_t i = j; i < matSize; ++i) {
            double sum = A[i][j];
            for (std::size_t k = 0; k < j; ++k) {
                sum -= A[i][k] * A[k][j];
            }
            A[i][j] = sum;

            double scaledSum = implicitScalingFactor[i] * std::abs(sum);
            if (scaledSum > pivot) {
                pivot = scaledSum;
                iMax = i;
            }
        }

        if (j != iMax) {
            std::swap(A[j], A[iMax]);
            parity = -parity;
            std::swap(implicitScalingFactor[j], implicitScalingFactor[iMax]);
        }

        ind[j] = iMax;

        if (std::abs(A[j][j]) < TINY) {
            A[j][j] = TINY;
        }

        if (j != matSize - 1) {
            double divisor = 1.0 / A[j][j];
            for (std::size_t i = j + 1; i < matSize; ++i) {
                A[i][j] *= divisor;
            }
        }
    }
}

void DenseSystemSolver::substituteBackwardLU(const Mat& A, const IndArr& ind, ColVect& B) {
    std::size_t matSize = A.size();
    int firstNonNullInd = -1;

    for (std::size_t i = 0; i < matSize; ++i) {
        std::size_t rowPermInd = ind[i];
        double sum = B[rowPermInd];
        B[rowPermInd] = B[i];

        if (firstNonNullInd != -1) {
            for (std::size_t j = firstNonNullInd; j < i; ++j) {
                sum -= A[i][j] * B[j];
            }
        } else if (sum != 0.0) {
            firstNonNullInd = static_cast<int>(i);
        }
        B[i] = sum;
    }

    for (int i = matSize - 1; i >= 0; --i) {
        double sum = B[i];
        for (std::size_t j = i + 1; j < matSize; ++j) {
            sum -= A[i][j] * B[j];
        }
        B[i] = sum / A[i][i];
    }
}

void DenseSystemSolver::solveDenseMatrixSystem(Mat& A, const ColVect& B, ColVect& Solution) {
    std::size_t matSize = A.size();
    if (matSize == 0 || A[0].size() != matSize || B.size() != matSize) {
        throw std::invalid_argument("Matrix dimensions must be consistent.");
    }

    IndArr ind(matSize);
    int parity = 1;
    bool isSingular = false;

    decomposeLU(A, ind, parity, isSingular);
    if (isSingular) {
        throw std::runtime_error("Matrix is singular.");
    }

    Solution = B;
    substituteBackwardLU(A, ind, Solution);
}
