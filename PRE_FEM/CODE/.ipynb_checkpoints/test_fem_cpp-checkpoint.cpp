// #include "grid.cpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <iomanip>

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

    InverseJacobian.elements[index(0, 0, InverseJacobian.col)] = Jacobian.elements[index(1, 1, Jacobian.col)] / DetJacobian;
    InverseJacobian.elements[index(0, 1, InverseJacobian.col)] = - Jacobian.elements[index(0, 1, Jacobian.col)] / DetJacobian;
    InverseJacobian.elements[index(1, 0, InverseJacobian.col)] = - Jacobian.elements[index(1, 0, Jacobian.col)] / DetJacobian;
    InverseJacobian.elements[index(1, 1, InverseJacobian.col)] = Jacobian.elements[index(0, 0, Jacobian.col)] / DetJacobian;

    matrix ShapeFunctionMatrix = multiplyMatrix(InverseJacobian, GradientMatrix);
    
    matrix LocalShapeFunctionMatrix = multiplyTransposeMatrixMatrix(ShapeFunctionMatrix);
    LocalShapeFunctionMatrix = multiplyScalarMatrix(DetJacobian, LocalShapeFunctionMatrix);
    
    return LocalShapeFunctionMatrix;
}


matrix calculateLocalStiffnessMatrixRectangular(double x, double y){
    matrix LocalStiffnessMatrix = initializeMatrix(4, 4);
    
    LocalStiffnessMatrix.elements[index(0, 0, LocalStiffnessMatrix.col)] = (x/y + y/x)/3.0;
    LocalStiffnessMatrix.elements[index(0, 1, LocalStiffnessMatrix.col)] = (x/y - 2.0*y/x)/6.0;
    LocalStiffnessMatrix.elements[index(0, 2, LocalStiffnessMatrix.col)] = -(x/y + y/x)/6.0;
    LocalStiffnessMatrix.elements[index(0, 3, LocalStiffnessMatrix.col)] = (y/x - 2.0*y/x)/6.0;

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
        return calculateLocalStiffnessMatrixRectangular(x_2 - x_1, y_3 - y_2);
    
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
        // if (M_2 == 0) {
        //     L_2 = L_1;
        //     H_2 = 0;
        // }
        // if (M_1 == 0) {
        //     H_1 = 0;
        // }
        // if (M_3 == 0) {
        //     H_3 = 0;
        // }
    }

    void generateGrid() {
        // if (!M_1 == 0) {
            for (unsigned int i = 0; i <= M_1; ++i) {
                for (unsigned int j = 0; j <= N; ++j) {
                    x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = H_1/M_1 * i;
                    y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = L_1/N * j;
                }
            }
        // }

        // if (!M_2 == 0) {
            for (unsigned int i = M_1 + 1; i < M_1 + M_2; ++i) {
                for (unsigned int j = 0; j <= N; ++j) {
                    double x_2 = H_1 + (H_2/M_2) * (i - M_1);
                    x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = x_2;
                    y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = (L_2 - L_1)/(H_2*N) * j * x_2 + (L_1/N) * j - (L_2-L_1)/(H_2*N) * j * H_1;
                }
            }
        // }

        for (unsigned int i = M_1 + M_2; i <= M_1 + M_2 + M_3; ++i) {
            for (unsigned int j = 0; j <= N; ++j) {
                x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = H_1 + H_2 + (H_3/M_3) * (i - M_1 - M_2);
                y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] = (L_2 / N) * j;
            }
        }
    }

    void displayGrid() const {
        std::cout << "x Grid Coordinates:" << std::endl;
        std::cout << N << endl;
        for (unsigned int i = 0; i <= M_1 + M_2 + M_3; ++i) {
            for (unsigned int j = 0; j <= N; ++j) {
                std::cout << std::setw(5) << std::fixed << std::setprecision(2) << x.elements[index(i, j, M_1 + M_2 + M_3 + 1)] << "  ";
            }
            std::cout << std::endl;
        }

        std::cout << "y Grid Coordinates:" << std::endl;
        for (unsigned int i = 0; i <= M_1 + M_2 + M_3; ++i) {
            for (unsigned int j = 0; j <= N; ++j) {
                std::cout << std::setw(5) << std::fixed << std::setprecision(2) << y.elements[index(i, j, M_1 + M_2 + M_3 + 1)] << "  ";
            }
            std::cout << std::endl;
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
        
        unsigned  N;              // Horizontal divisions
        unsigned int M_1, M_2, M_3;     // Vertical divisions for each region
        double H_1, H_2, H_3;  // Widths of the regions
        double L_1, L_2;      // Heights of the regions
        
        StructuredGrid GridCoordinates;
        
        GlobalStiffnessMatrix(double H_1, double H_2, double H_3, double L_1, double L_2, unsigned int N, unsigned int M_1, unsigned int M_2, unsigned int M_3) : hNode(hNode), vNode(vNode), H_1(H_1), H_2(H_2), H_3(H_3), L_1(L_1), L_2(L_2), N(N), M_1(M_1), M_2(M_2), M_3(M_3), GridCoordinates(H_1, H_2, H_3, L_1, L_2, N, M_1, M_2, M_3) {
            this->hNode = M_3 + 1;
            this->vNode = N + 1;
            this->nNode = this->hNode * this->vNode;
            this->hElem = this->hNode - 1;
            this->vElem = this->vNode - 1;
            this->nElem = this->hElem * this->vElem;
            GridCoordinates.generateGrid();

            GridCoordinates.displayGrid();

        }

        // hNode = number of columns
        // vNode = number of rows
        matrix assembleGlobalStiffnessMatrix(){
            matrix result;

            // Initialize K_g as a zero matrix of size (N x N)
            result = initializeMatrix(vNode, hNode);

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

            for (unsigned int i = 0; i < hElem; i++) {
                for (unsigned int j = 0; j < vElem; j++) {
                    // retrieve the zero node coordinates
                    unsigned int ZeroIndex = ZeroNodeIndex(i, j, GridCoordinates.M_3);
                    double x_0 = GridCoordinates.x.elements[ZeroIndex];
                    double y_0 = GridCoordinates.x.elements[ZeroIndex];
                    
                    unsigned int FirstIndex = FirstNodeIndex(i, j, GridCoordinates.M_3);
                    double x_1 = GridCoordinates.x.elements[FirstIndex];
                    double y_1 = GridCoordinates.x.elements[FirstIndex];
                    
                    unsigned int SecondIndex = SecondNodeIndex(i, j, GridCoordinates.M_3);
                    double x_2 = GridCoordinates.x.elements[SecondIndex];
                    double y_2 = GridCoordinates.x.elements[SecondIndex];
                    
                    unsigned int ThirdIndex = ThirdNodeIndex(i, j, GridCoordinates.M_3);
                    double x_3 = GridCoordinates.x.elements[ThirdIndex];
                    double y_3 = GridCoordinates.x.elements[ThirdIndex];

                    // calculate the local stiffness matrix for the current element

                    matrix A = calculateLocalStiffnessMatrix(x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3);
                    for (unsigned int k = 0; k < A.row; k++){
                        for (unsigned int l = 0; l < A.col; l++){
                            result.elements[index(k, l, A.col)] += A.elements[index(k, l, A.col)];
                            result.elements[transpose_index(k, l, A.col)] += A.elements[index(k, l, A.col)];
                        }
                    }
                    
                }
            }
            return result;
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
    double H_1 = 3.0, H_2 = 3.0, H_3 = 3.0, L_1 = 4.0, L_2 = 2.0;
    unsigned int N = 5, M_1 = 5, M_2 = 5, M_3 = 5;

    GlobalStiffnessMatrix globalMatrix(H_1, H_2, H_3, L_1, L_2, N, M_1, M_2, M_3);

    matrix A = globalMatrix.assembleGlobalStiffnessMatrix();

    printMatrix(A);

}
