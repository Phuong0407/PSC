// #include "grid.cpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <iomanip>
#include <fstream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#define index(i, j, col) ((i) * (col) + (j))
#define transpose_index(i, j, row) ((j) * (row) + (i))

using namespace std;
struct matrix{
    unsigned int row = 0, col = 0;
    double* elements;
};

matrix initializeMatrix(unsigned int row, int col) {
    matrix _matrix;
    _matrix.row = row;
    _matrix.col = col;
    _matrix.elements = new double[row * col];
    for (unsigned int i = 0; i < row * col; ++i)
        _matrix.elements[i] = 0.0;
    return _matrix;
}


void printMatrix(const matrix mat) {
    for (unsigned int i = 0; i < mat.row; i++) {
        for (unsigned int j = 0; j < mat.col; j++) {
            std::cout << mat.elements[index(i, j, mat.col)] << " ";
        }
        std::cout << std::endl;
    }
}


matrix addMatrix(const matrix A, const matrix B) {
    if (A.row != B.row || A.col != B.col)
        throw std::invalid_argument("Matrix addition not possible: Dimensions of A and B must match");

    matrix _matrix = initializeMatrix(A.row, A.col);

    for (unsigned int i = 0; i < A.row; i++) {
        for (unsigned int j = 0; j < A.col; j++) {
            _matrix.elements[index(i, j, A.col)] =
                A.elements[index(i, j, A.col)] + B.elements[index(i, j, A.col)];
        }
    }
    return _matrix;
}


matrix multiplyScalarMatrix(const double Scalar, const matrix B) {
    matrix _matrix = initializeMatrix(B.col, B.col);
    for (unsigned int i = 0; i < B.col; i++) {
        for (unsigned int j = 0; j < B.col; j++) {
            _matrix.elements[index(i, j, _matrix.col)] = Scalar * B.elements[index(i, j, B.col)];
        }
    }
    return _matrix;
}


matrix multiplyMatrix(const matrix A, const matrix B) {
    if (A.col != B.row)
        throw std::invalid_argument("Matrix multiplication not possible: A.col must equal B.row");

    matrix _matrix = initializeMatrix(A.row, B.col);
    for (unsigned int i = 0; i < A.row; i++) {
        for (unsigned int j = 0; j < B.col; j++) {
            for (unsigned int k = 0; k < A.col; k++) {
                _matrix.elements[index(i, j, B.col)] += A.elements[index(i, k, A.col)] * B.elements[index(k, j, B.col)];
            }
        }
    }
    return _matrix;
}

matrix multiplyTransposeMatrixMatrix(const matrix B) {
    matrix _matrix = initializeMatrix(B.col, B.col);
    for (unsigned int i = 0; i < B.col; i++) {
        for (unsigned int j = 0; j < B.col; j++) {
            for (unsigned int k = 0; k < B.row; k++) {
                _matrix.elements[index(i, j, B.col)] += B.elements[index(k, i, B.col)] * B.elements[index(k, j, B.col)];
            }
        }
    }
    return _matrix;
}

matrix initializeGradientMatrix(double xi, double eta){
    matrix result = initializeMatrix(2, 4);
    
    result.elements[index(0, 0, result.col)] = (eta - 1)/ 4.0;
    result.elements[index(0, 1, result.col)] = (1 - eta)/ 4.0;
    result.elements[index(0, 2, result.col)] = (1 + eta)/ 4.0;
    result.elements[index(0, 3, result.col)] = (-eta - 1)/ 4.0;

    result.elements[index(1, 0, result.col)] = (xi - 1)/ 4.0;
    result.elements[index(1, 1, result.col)] = (-xi - 1)/ 4.0;
    result.elements[index(1, 2, result.col)] = (1 + xi)/ 4.0;
    result.elements[index(1, 3, result.col)] = (1 - xi)/ 4.0;

    return result;
}

matrix initializeCoordinatesMatrix(double x_1, double y_1, double x_2, double y_2, double x_3, double y_3, double x_4, double y_4){
    matrix result = initializeMatrix(4, 2);

    result.elements[index(0, 0, result.col)] = x_1;
    result.elements[index(0, 1, result.col)] = y_1;

    result.elements[index(1, 0, result.col)] = x_2;
    result.elements[index(1, 1, result.col)] = y_2;

    result.elements[index(2, 0, result.col)] = x_3;
    result.elements[index(2, 1, result.col)] = y_3;

    result.elements[index(3, 0, result.col)] = x_4;
    result.elements[index(3, 1, result.col)] = y_4;

    return result;
}

matrix calculateLocalShapeFunction(double x_1, double y_1, double x_2, double y_2, double x_3, double y_3, double x_4, double y_4, double xi, double eta) {

    matrix GradientMatrix = initializeGradientMatrix(xi, eta);
    matrix CoordinatesMatrix = initializeCoordinatesMatrix(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4);
    matrix Jacobian = multiplyMatrix(GradientMatrix, CoordinatesMatrix);

    double DetJacobian = Jacobian.elements[index(0, 0, Jacobian.col)] * Jacobian.elements[index(1, 1, Jacobian.col)]
                        - Jacobian.elements[index(0, 1, Jacobian.col)] * Jacobian.elements[index(1, 0, Jacobian.col)];

    matrix InverseJacobian = initializeMatrix(2, 2);

    double J_11 = Jacobian.elements[index(0, 0, InverseJacobian.col)];
    InverseJacobian.elements[index(0, 0, InverseJacobian.col)] = Jacobian.elements[index(1, 1, Jacobian.col)] / DetJacobian;
    InverseJacobian.elements[index(0, 1, InverseJacobian.col)] = - Jacobian.elements[index(0, 1, Jacobian.col)] / DetJacobian;
    InverseJacobian.elements[index(1, 0, InverseJacobian.col)] = - Jacobian.elements[index(1, 0, Jacobian.col)] / DetJacobian;
    InverseJacobian.elements[index(1, 1, InverseJacobian.col)] = J_11 / DetJacobian;

    matrix ShapeFunctionMatrix = multiplyMatrix(InverseJacobian, GradientMatrix);
    
    matrix LocalShapeFunctionMatrix = multiplyTransposeMatrixMatrix(ShapeFunctionMatrix);
    LocalShapeFunctionMatrix = multiplyScalarMatrix(DetJacobian, LocalShapeFunctionMatrix);
    
    return LocalShapeFunctionMatrix;
}


matrix calculateLocalShapeFunctionRectangular(double x, double y){
    matrix LocalStiffnessMatrix = initializeMatrix(4, 4);
    
    LocalStiffnessMatrix.elements[index(0, 0, LocalStiffnessMatrix.col)] = (x/y + y/x)/3.0;
    LocalStiffnessMatrix.elements[index(0, 1, LocalStiffnessMatrix.col)] = (x/y - 2.0*y/x)/6.0;
    LocalStiffnessMatrix.elements[index(0, 2, LocalStiffnessMatrix.col)] = -(x/y + y/x)/6.0;
    LocalStiffnessMatrix.elements[index(0, 3, LocalStiffnessMatrix.col)] = (y/x - 2.0*x/y)/6.0;

    LocalStiffnessMatrix.elements[index(1, 0, LocalStiffnessMatrix.col)] = (x/y - 2.0*y/x)/6.0;
    LocalStiffnessMatrix.elements[index(1, 1, LocalStiffnessMatrix.col)] = (x/y + y/x)/3.0;
    LocalStiffnessMatrix.elements[index(1, 2, LocalStiffnessMatrix.col)] = (y/x - 2.0*x/y)/6.0;
    LocalStiffnessMatrix.elements[index(1, 3, LocalStiffnessMatrix.col)] = -(x/y + y/x)/6.0;

    LocalStiffnessMatrix.elements[index(2, 0, LocalStiffnessMatrix.col)] = -(x/y + y/x)/6.0;
    LocalStiffnessMatrix.elements[index(2, 1, LocalStiffnessMatrix.col)] = (y/x - 2.0*x/y)/6.0;
    LocalStiffnessMatrix.elements[index(2, 2, LocalStiffnessMatrix.col)] = (x/y + y/x)/3.0;
    LocalStiffnessMatrix.elements[index(2, 3, LocalStiffnessMatrix.col)] = (x/y - 2.0*y/x)/6.0;

    LocalStiffnessMatrix.elements[index(3, 0, LocalStiffnessMatrix.col)] = (y/x - 2.0*x/y)/6.0;
    LocalStiffnessMatrix.elements[index(3, 1, LocalStiffnessMatrix.col)] = -(x/y + y/x)/6.0;
    LocalStiffnessMatrix.elements[index(3, 2, LocalStiffnessMatrix.col)] = (x/y - 2.0*y/x)/6.0;
    LocalStiffnessMatrix.elements[index(3, 3, LocalStiffnessMatrix.col)] = (x/y + y/x)/3.0;

    return LocalStiffnessMatrix;
}


matrix calculateLocalStiffnessMatrix(double x_1, double y_1, double x_2, double y_2, double x_3, double y_3, double x_4, double y_4) {
    if(x_1 == x_4 && x_2 == x_3 && y_1 == y_2 && y_3 == y_4)
        return calculateLocalShapeFunctionRectangular(x_2 - x_1, y_3 - y_2);
    
    matrix ShapeMatrix_11 = calculateLocalShapeFunction(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, -1.0/sqrt(3), -1.0/sqrt(3));
    matrix ShapeMatrix_12 = calculateLocalShapeFunction(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, -1.0/sqrt(3), 1.0/sqrt(3));
    matrix ShapeMatrix_21 = calculateLocalShapeFunction(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, 1.0/sqrt(3), -1.0/sqrt(3));
    matrix ShapeMatrix_22 = calculateLocalShapeFunction(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4, 1.0/sqrt(3), 1.0/sqrt(3));

    matrix LocalStiffnessMatrix = addMatrix(ShapeMatrix_11, ShapeMatrix_12);
    LocalStiffnessMatrix = addMatrix(LocalStiffnessMatrix, ShapeMatrix_21);
    LocalStiffnessMatrix = addMatrix(LocalStiffnessMatrix, ShapeMatrix_22);

    return LocalStiffnessMatrix;
}

class StructuredGrid {
public:
    double H_1, H_2, H_3;           // Widths of the regions
    double L_1, L_2;                // Heights of the regions
    unsigned int N;                 // Vertical Number of Elements divisions
    // Horizontal Number of Elements for each region
    unsigned int M_1;               // M_1 = First Rectangular Regions, two verticals boundary included
    unsigned int M_2;               // M_2 = Oblique Regions, two verticals boundary excluded
    unsigned int M_3;               // M_3 = Second Rectangular Regions, two verticals boundary included

    matrix x, y;

    StructuredGrid(double H_1, double H_2, double H_3, double L_1, double L_2, int N, int M_1, int M_2, int M_3) : H_1(H_1), H_2(H_2), H_3(H_3), L_1(L_1), L_2(L_2), N(N), M_1(M_1), M_2(M_2), M_3(M_3) {
        x = initializeMatrix(N + 1, M_1 + M_2 + M_3 + 1);
        y = initializeMatrix(N + 1, M_1 + M_2 + M_3 + 1);
    }

    void generateGrid() {
        for (unsigned int i = 0; i <= N; ++i) {
            for (unsigned int j = 0; j <= M_1; ++j) {
                x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = H_1/M_1 * j;
                y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = L_1/N * i;
            }
            for (unsigned int j = M_1 + 1; j < M_1 + M_2; ++j) {
                double x_2 = H_1 + (H_2/M_2) * (j - M_1);
                x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = x_2;
                y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = (L_2 - L_1)/(H_2*N) * i * x_2 + (L_1/N) * i - (L_2-L_1)/(H_2*N) * i * H_1;
            }
            for (unsigned int j = M_1 + M_2; j <= M_1 + M_2 + M_3; ++j) {
                x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = H_1 + H_2 + (H_3/M_3) * (j - M_1 - M_2);
                y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = (L_2 / N) * i;
            }
        }
    }

    void displayGrid() const {
        std::cout << "x Grid Coordinates:" << std::endl;
        std::cout << N << endl;
        for (unsigned int i = 0; i <= M_1 + M_2 + M_3; ++i) {
            for (unsigned int j = 0; j <= N; ++j) {
                std::cout << std::setw(2) << std::fixed << std::setprecision(2) << "(" << x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] << ", " << y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] << ") ";
            }
            std::cout << std::endl;
        }
        // Save the grid coordinates to a text file
        std::ofstream file("structured_grid_coordinates_cpp.txt");
        if (file.is_open()) {
            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j <= M_1 + M_2 + M_3; ++j) {
                    file << "(" << std::fixed << std::setprecision(2) << x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] << ", " << y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] << ")-";
                }
                file << "\n";
            }
            file.close();
            std::cout << "Coordinates have been saved to 'structured_grid_coordinates.txt'." << std::endl;
        } else {
            std::cerr << "Unable to open file!" << std::endl;
        }

    }
};


inline unsigned int ZeroNodeIndex(unsigned int i, unsigned int j, unsigned int hNode){
    return i * hNode + j;
}
inline unsigned int FirstNodeIndex(unsigned int i, unsigned int j, unsigned int hNode){
    return i * hNode + j + 1;
}
inline unsigned int SecondNodeIndex(unsigned int i, unsigned int j, unsigned int hNode){
    return (i + 1) * hNode + j + 1;
}
inline unsigned int ThirdNodeIndex(unsigned int i, unsigned int j, unsigned int hNode){
    return (i + 1) * hNode + j;
}



class GlobalStiffnessMatrix{
    public:
        unsigned int nElem; // number of elements
        unsigned int hElem; // number of horizontal elements
        unsigned int vElem; // number of vertical elements
        unsigned int nNode; // number of nodes
        unsigned int hNode; // number of horizontal nodes
        unsigned int vNode; // number of vertical nodes
        
        unsigned  N;                    // Horizontal divisions
        unsigned int M_1, M_2, M_3;     // Vertical divisions for each region
        double H_1, H_2, H_3;           // Widths of the regions
        double L_1, L_2;                // Heights of the regions
        double U_1, U_2;                // Velocity
        
        StructuredGrid GridCoordinates;
        matrix K_g;
        matrix RHS;
        matrix Node;
        matrix Vx;
        matrix Vy;

        GlobalStiffnessMatrix(double H_1, double H_2, double H_3, double L_1, double L_2, double U_1, double U_2, unsigned int N, unsigned int M_1, unsigned int M_2, unsigned int M_3) : hNode(hNode), vNode(vNode), H_1(H_1), H_2(H_2), H_3(H_3), L_1(L_1), L_2(L_2), U_1(U_1), U_2(U_2), N(N), M_1(M_1), M_2(M_2), M_3(M_3), GridCoordinates(H_1, H_2, H_3, L_1, L_2, N, M_1, M_2, M_3) {
            this->hNode = M_1 + M_2 + M_3 + 1;
            this->vNode = N + 1;
            this->nNode = this->hNode * this->vNode;
            this->hElem = this->hNode - 1;
            this->vElem = this->vNode - 1;
            this->nElem = this->hElem * this->vElem;
            GridCoordinates.generateGrid();
            K_g = assembleGlobalStiffnessMatrix();
            RHS = calculateRightHandSideMatrix(U_1, U_2, L_1 / N, L_2 / N);
            Node = solveNodalValue();
            calculateVelocity();
        }

        matrix assembleGlobalStiffnessMatrix(){
            matrix result;

            // Initialize K_g as a zero matrix of size (N x N)
            result = initializeMatrix(nNode, nNode);

            // FOR each element e = 1 to n:
            //     Retrieve element stiffness matrix K_e (m x m)
            //     Retrieve connectivity for element e, C_e (list of m global node indices)
                
            //     FOR i = 1 to m:
            //         FOR j = 1 to m:
            //             global_i = C_e[i]  # Map local node i to global node
            //             global_j = C_e[j]  # Map local node j to global node
                        
            //             K_g[global_i, global_j] += K_e[i, j]  # Add local stiffness contribution
            //         END FOR
            //     END FOR
            // END FOR

            // Return K_g

            for (unsigned int i = 0; i < vElem; i++) {
                for (unsigned int j = 0; j < hElem; j++) {
                    // retrieve the zero node coordinates
                    unsigned int ZeroIndex = ZeroNodeIndex(i, j, GridCoordinates.M_1 + GridCoordinates.M_2 + GridCoordinates.M_3 + 1);
                    double x_0 = GridCoordinates.x.elements[ZeroIndex];
                    double y_0 = GridCoordinates.y.elements[ZeroIndex];
                    
                    unsigned int FirstIndex = FirstNodeIndex(i, j, GridCoordinates.M_1 + GridCoordinates.M_2 + GridCoordinates.M_3 + 1);
                    double x_1 = GridCoordinates.x.elements[FirstIndex];
                    double y_1 = GridCoordinates.y.elements[FirstIndex];
                    
                    unsigned int SecondIndex = SecondNodeIndex(i, j, GridCoordinates.M_1 + GridCoordinates.M_2 + GridCoordinates.M_3 + 1);
                    double x_2 = GridCoordinates.x.elements[SecondIndex];
                    double y_2 = GridCoordinates.y.elements[SecondIndex];
                    
                    unsigned int ThirdIndex = ThirdNodeIndex(i, j, GridCoordinates.M_1 + GridCoordinates.M_2 + GridCoordinates.M_3 + 1);
                    double x_3 = GridCoordinates.x.elements[ThirdIndex];
                    double y_3 = GridCoordinates.y.elements[ThirdIndex];

                    // calculate the local stiffness matrix for the current element

                    matrix A = calculateLocalStiffnessMatrix(x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3);
                    // unsigned int a_10 = index(FirstIndex, SecondIndex, result.col);
                    result.elements[index(ZeroIndex, ZeroIndex, result.col)] += A.elements[index(0, 0, A.col)];
                    result.elements[index(ZeroIndex, FirstIndex, result.col)] += A.elements[index(0, 1, A.col)];
                    result.elements[index(ZeroIndex, SecondIndex, result.col)] += A.elements[index(0, 2, A.col)];
                    result.elements[index(ZeroIndex, ThirdIndex, result.col)] += A.elements[index(0, 3, A.col)];

                    result.elements[index(FirstIndex, ZeroIndex, result.col)] += A.elements[index(1, 0, A.col)];
                    result.elements[index(FirstIndex, FirstIndex, result.col)] += A.elements[index(1, 1, A.col)];
                    result.elements[index(FirstIndex, SecondIndex, result.col)] += A.elements[index(1, 2, A.col)];
                    result.elements[index(FirstIndex, ThirdIndex, result.col)] += A.elements[index(1, 3, A.col)];
                    
                    result.elements[index(SecondIndex, ZeroIndex, result.col)] += A.elements[index(2, 0, A.col)];
                    result.elements[index(SecondIndex, FirstIndex, result.col)] += A.elements[index(2, 1, A.col)];
                    result.elements[index(SecondIndex, SecondIndex, result.col)] += A.elements[index(2, 2, A.col)];
                    result.elements[index(SecondIndex, ThirdIndex, result.col)] += A.elements[index(2, 3, A.col)];
                    
                    result.elements[index(ThirdIndex, ZeroIndex, result.col)] += A.elements[index(3, 0, A.col)];
                    result.elements[index(ThirdIndex, FirstIndex, result.col)] += A.elements[index(3, 1, A.col)];
                    result.elements[index(ThirdIndex, SecondIndex, result.col)] += A.elements[index(3, 2, A.col)];
                    result.elements[index(ThirdIndex, ThirdIndex, result.col)] += A.elements[index(3, 3, A.col)];

                    // for (unsigned int k = 0; k < A.row; k++){
                    //     for (unsigned int l = 0; l < A.col; l++){
                    //         result.elements[index(k, l, A.col)] += A.elements[index(k, l, A.col)];
                    //         result.elements[transpose_index(k, l, A.col)] += A.elements[index(k, l, A.col)];
                    //     }
                    // }
                    // printMatrix(A);
                    // cout << "Local Stiffness Matrix of the " << i * hElem + j << " elements:" << endl;
                    // printMatrix(result);
                    // // printMatrix(A);
                    // cout << "End of Local Stiffness Matrix of the " << i * hElem + j << " elements:" << endl << endl;

                    //         std::ofstream file("global_stiffness_matrix_output.txt", std::ios::app);
                    //         std::ofstream file_2("local_stiffness_matrix_output.txt", std::ios::app);
                    //         if (file.is_open()) {
                    //             // Write a message about the matrix
                    //             file << "Local Stiffness Matrix of the " << i * hElem + j << " elements:" << std::endl;
                    //             file_2 << "Local Stiffness Matrix of the " << i * hElem + j << " elements:" << std::endl;

                    //             // Write the matrix to the file using nested loops
                    //             for (unsigned int row = 0; row < result.row; row++) {
                    //                 for (unsigned int col = 0; col < result.col; col++) {
                    //                     file << result.elements[index(row, col, result.col)] << " ";  // Write each element
                    //                 }
                    //                 file << std::endl;  // New line after each row
                    //             }

                    //             for (unsigned int row = 0; row < A.row; row++) {
                    //                 for (unsigned int col = 0; col < A.col; col++) {
                    //                     file_2 << A.elements[index(row, col, A.col)] << " ";  // Write each element
                    //                 }
                    //                 file_2 << std::endl;  // New line after each row
                    //             }


                    //             // End message
                    //             file << "End of Local Stiffness Matrix of the " << i * hElem + j << " elements" << std::endl << std::endl;
                    //             file_2 << "End of Local Stiffness Matrix of the " << i * hElem + j << " elements" << std::endl << std::endl;


                    //             file.close();
                    //             file_2.close();
                    //             std::cout << "Matrix data has been saved to 'stiffness_matrix_output.txt'." << std::endl;
                    //         }
                    //         else {
                    //                 std::cerr << "Unable to open file for writing!" << std::endl;
                    //         }

                }
            }
            return result;
        }

        matrix calculateRightHandSideMatrix(double U_1, double U_2, double dy_1, double dy_2) {
            matrix result = initializeMatrix(nNode, 1);

            for(unsigned int i = 0; i < vElem; ++i) {
                if(i != 0 && i != vElem - 1){
                    result.elements[ZeroNodeIndex(i, 0, hNode)] -= U_1 * dy_1 / 2.0;
                    result.elements[ThirdNodeIndex(i, 0, hNode)] -= U_1 * dy_1 / 2.0;

                    result.elements[FirstNodeIndex(i, hElem, hNode)] += U_2 * dy_2 / 2.0;
                    result.elements[SecondNodeIndex(i, hElem, hNode)] += U_2 * dy_2 / 2.0;
                }
                else if(i == 0) {
                    result.elements[ThirdNodeIndex(i, 0, result.col)] -= U_1 * dy_1 / 2.0;
                    result.elements[SecondNodeIndex(i, hElem, hNode)] += U_2 * dy_2 / 2.0;
                }
                else {
                    result.elements[ZeroNodeIndex(i, 0, result.col)] -= U_1 * dy_1 / 2.0;
                    result.elements[ZeroNodeIndex(i, hElem, result.col)] += U_2 * dy_2 / 2.0;
                }
            }
            return result;
        }

    matrix solveNodalValue(){
        matrix result = initializeMatrix(nNode, 1);

        Eigen::MatrixXd eK_g(K_g.row, K_g.col);
        Eigen::VectorXd eRHS(RHS.row);
        for (unsigned int i = 0; i < K_g.row; ++i) {
            eRHS(i) = RHS.elements[index(i, 0, RHS.col)];
            for (unsigned int j = 0; j < K_g.col; ++j) {
                eK_g(i, j) = K_g.elements[index(i, j, K_g.col)];
            }
        }

        Eigen::SparseMatrix<double> sparseK_g = eK_g.sparseView();

        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
        solver.compute(sparseK_g); // Use the sparse matrix
        if (sparseK_g.rows() != sparseK_g.cols()) {
            std::cerr << "Matrix is not square!" << std::endl;
            return result;
        }
        if (!sparseK_g.isApprox(sparseK_g.transpose())) {
            std::cerr << "Matrix is not symmetric!" << std::endl;
            return result;
        }
        
        if (solver.info() != Eigen::Success) {
            std::cerr << "Decomposition failed!" << std::endl;
            return result;
        }


        for (int i = 0; i < sparseK_g.rows(); ++i) {
            if (sparseK_g.row(i).norm() == 0) {
                std::cerr << "Row " << i << " is all zeros!" << std::endl;
            }
        }


        Eigen::VectorXd u = solver.solve(eRHS);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Solve failed!" << std::endl;
            return result;
        }

        for (unsigned int i = 0; i < K_g.row; ++i) {
            result.elements[index(i, 0, result.col)] = u(i);
        }
        return result;
    }

    void calculateVelocity() {
        matrix Vx_result = initializeMatrix(nNode, 1);
        matrix Vy_result = initializeMatrix(nNode, 1);

        for (unsigned int i = 0; i < vNode; ++i) {
            for (unsigned int j = 0; j < hNode; ++j) {
                if (i != 0 && i != vNode - 1 && j != 0 && j != hNode - 1) {
                    Vx_result.elements[index(i, j, hNode)] = (Node.elements[index(i, j + 1, hNode)] - Node.elements[index(i, j, hNode)])
                                                                / (GridCoordinates.x.elements[index(i, j + 1, hNode)] - GridCoordinates.x.elements[index(i, j, hNode)]);
                    Vy_result.elements[index(i, j, hNode)] = (Node.elements[index(i + 1, j, hNode)] - Node.elements[index(i, j, hNode)])
                                                                / (GridCoordinates.y.elements[index(i + 1, j, hNode)] - GridCoordinates.y.elements[index(i, j, hNode)]);
                }
                else if (j == 0) {
                    Vx_result.elements[index(i, j, hNode)] = U_1;
                }
                else if (j == hNode - 1) {
                    Vx_result.elements[index(i, j, hNode)] = U_2;
                }
                else if (i == 0) {
                    Vy_result.elements[index(i, j, hNode)] = 0.0;
                }
                else if (i == vNode - 1) {
                    Vy_result.elements[index(i, j, hNode)] = U_2;
                }
            }
        }
    }
};







int main(){
    // matrix A = calculateLocalStiffnessMatrix(0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5);
    // printMatrix(A);
    // matrix B = calculateLocalStiffnessMatrix(0.0, 0.0, 5.0, 0.0, 5.0, 1.0, 0.0, 1.0);
    // printMatrix(B);
    // cout << "(" << ZeroNodeIndex(3, 3, 4) << ", "
    //             << FirstNodeIndex(3, 3, 4) << ", "
    //             << SecondNodeIndex(3, 3, 4) << ", "
    //             << ThirdNodeIndex(3, 3, 4) << ")" << endl;


    double H_1 = 1.0, H_2 = 1.0, H_3 = 1.0, L_1 = 1.0, L_2 = 2.0, U_1 = 10.0, U_2 = 5.0; 
    unsigned int N = 5, M_1 = 5, M_2 = 5, M_3 = 5;

    // StructuredGrid A(H_1, H_2, H_3, L_1, L_2, N, M_1, M_2, M_3);
    // A.generateGrid();
    // A.displayGrid();


    GlobalStiffnessMatrix globalMatrix(H_1, H_2, H_3, L_1, L_2, U_1, U_2, N, M_1, M_2, M_3);

    // matrix A = globalMatrix.assembleGlobalStiffnessMatrix();

    // printMatrix(globalMatrix.Node);
    printMatrix(globalMatrix.Vx);
    // printMatrix(globalMatrix.Vy);
    
    
    // std::ofstream file("stiffness_matrix_output.txt", std::ios::app);
    // if (file.is_open()) {
    //     // Write a message about the matrix
    //     file << "Global Stiffness Matrix:" << std::endl;

    //     // Write the matrix to the file using nested loops
    //     for (unsigned int row = 0; row < A.row; row++) {
    //         for (unsigned int col = 0; col < A.col; col++) {
    //             file << A.elements[index(row, col, A.col)] << " ";  // Write each element
    //         }
    //         file << std::endl;  // New line after each row
    //     }

    //     // End message
    //     file << "End of Global Stiffness Matrix." << std::endl << std::endl;

    //     file.close();
    //     std::cout << "Matrix data has been saved to 'stiffness_matrix_output.txt'." << std::endl;
    // }
    // else {
    //         std::cerr << "Unable to open file for writing!" << std::endl;
    // }


    // printMatrix(calculateLocalStiffnessMatrix(3, 0, 3.6, 0, 3.6, 0.72, 3, 0.8));
    // printMatrix(calculateLocalStiffnessMatrix(0, 1, 0, 0, 2, 0.5, 2, 1));
    // printMatrix(calculateLocalStiffnessMatrix(0, 0, 2, 0, 1, 1.5, 0, 1.1));


}
