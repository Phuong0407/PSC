#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>
#include <iomanip>

class StructuredGrid {
public:
    // Parameters for the grid
    double H1, H2, H3;  // Widths of the regions
    double L1, L2;      // Heights of the regions
    unsigned  N;              // Horizontal divisions
    unsigned int M1, M2, M3;     // Vertical divisions for each region

    // 2D vectors to store the grid coordinates
    std::vector<std::vector<double>> X;
    std::vector<std::vector<double>> Y;

    // Constructor
    StructuredGrid(double h1, double h2, double h3, double l1, double l2, int n, int m1, int m2, int m3)
        : H1(h1), H2(h2), H3(h3), L1(l1), L2(l2), N(n), M1(m1), M2(m2), M3(m3),
          X(m3 + 1, std::vector<double>(n + 1, 0.0)),
          Y(m3 + 1, std::vector<double>(n + 1, 0.0)) {}

    // Generate the grid
    void generateGrid() {
        // Region 1: Left rectangular domain
        for (int i = 0; i <= M1; ++i) {
            for (int j = 0; j <= N; ++j) {
                X[i][j] = (H1 / M1) * i;
                Y[i][j] = (L1 / N) * j;
            }
        }

        // Region 2: Oblique domain
        for (int i = M1 + 1; i <= M2; ++i) {
            for (int j = 0; j <= N; ++j) {
                double h_middle = H1 + (H2 / (M2 - M1)) * (i - M1);
                X[i][j] = h_middle;
                Y[i][j] = ((L2 - L1) / (N * H2)) * j * h_middle +
                          (L1 / N) * j -
                          ((L2 - L1) / N) * (H1 / H2) * j;
            }
        }

        // Region 3: Right rectangular domain
        for (int i = M2 + 1; i <= M3; ++i) {
            for (int j = 0; j <= N; ++j) {
                X[i][j] = H1 + H2 + (H3 / (M3 - M2)) * (i - M2);
                Y[i][j] = (L2 / N) * j;
            }
        }
    }

    // Display the grid points
    void displayGrid() const {
        std::cout << "X Grid Coordinates:" << std::endl;
        for (const auto& row : X) {
            for (const auto& val : row) {
                std::cout << std::setw(8) << std::fixed << std::setprecision(2) << val << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "\nY Grid Coordinates:" << std::endl;
        for (const auto& row : Y) {
            for (const auto& val : row) {
                std::cout << std::setw(8) << std::fixed << std::setprecision(2) << val << " ";
            }
            std::cout << std::endl;
        }
    }
};

#endif

// int main() {
//     // Parameters
//     double H1 = 3.0, H2 = 4.0, H3 = 8.0;
//     double L1 = 8.0, L2 = 4.0;
//     int N = 10, M1 = 10, M2 = 20, M3 = 30;

//     // Create the grid object
//     StructuredGrid grid(H1, H2, H3, L1, L2, N, M1, M2, M3);

//     // Generate the grid
//     grid.generateGrid();

//     // Display the grid coordinates
//     grid.displayGrid();

//     return 0;
// }                            result.elements[index(k, l, A.col)] = A.elements[index(k, l, A.col)];
