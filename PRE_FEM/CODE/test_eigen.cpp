#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

int main() {
    // Dense matrix example
    Eigen::MatrixXd mat(2, 2);
    mat(0, 0) = 3;
    mat(1, 0) = 2.5;
    mat(0, 1) = -1;
    mat(1, 1) = mat(1, 0) + mat(0, 1);
    std::cout << "Dense matrix:\n" << mat << "\n\n";

    // Sparse matrix example
    Eigen::SparseMatrix<double> spMat(3, 3);
    spMat.insert(0, 0) = 1;
    spMat.insert(1, 1) = 2;
    spMat.insert(2, 2) = 3;
    std::cout << "Sparse matrix:\n" << Eigen::MatrixXd(spMat) << "\n\n";

    // Iterative solver example
    Eigen::SparseMatrix<double> A(3, 3);
    A.insert(0, 0) = 4; A.insert(0, 1) = -1; A.insert(1, 1) = 4; A.insert(1, 2) = -1; A.insert(2, 2) = 4;
    Eigen::VectorXd b(3);
    b << 1, 2, 3;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(b);
    std::cout << "Solution x:\n" << x << std::endl;

    return 0;
}