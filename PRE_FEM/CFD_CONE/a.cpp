#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

// Define a 2D vector for grid representation
using Grid = std::vector<std::vector<double>>;


// Perform elliptic mesh generation
void generateEllipticMesh(int nx, int ny, double Ermax, bool show) {
    Grid x(nx, std::vector<double>(ny, 0.0));
    Grid y(nx, std::vector<double>(ny, 0.0));
    Grid new_x = x, new_y = y;
    Grid alpha(nx, std::vector<double>(ny, 0.0));
    Grid beta(nx, std::vector<double>(ny, 0.0));
    Grid gamma(nx, std::vector<double>(ny, 0.0));

    std::size_t max_iter = 1000;

    std::vector<double> Er1(max_iter, 0.0);
    std::vector<double> Er2(max_iter, 0.0);

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        for (std::size_t i = 1; i < N - 1; ++i) {
            for (std::size_t j = 1; j < M - 1; ++j) {
                alpha[i][j] = 0.25 * (pow(x[i][j + 1] - x[i][j - 1], 2.0) + pow(y[i][j + 1] - y[i][j - 1], 2.0));
                beta[i][j] = 0.0625 * ((x[i + 1][j] - x[i - 1][j]) * (x[i][j + 1] - x[i][j - 1]) + (y[i + 1][j] - y[i - 1][j]) * (y[i][j + 1] - y[i][j - 1]));
                gamma[i][j] = 0.25 * (pow(x[i + 1][j] - x[i - 1][j], 2.0) + pow(y[i + 1][j] - y[i - 1][j], 2.0));

                double denominator = alpha[i][j] + gamma[i][j] + 1e-9;
                new_x[i][j] = (-0.5 / denominator) * 
                             (2.0 * beta[i][j] * (x[i + 1][j + 1] - x[i - 1][j + 1] - x[i + 1][j - 1] + x[i - 1][j - 1]) -
                              alpha[i][j] * (x[i + 1][j] + x[i - 1][j]) -
                              gamma[i][j] * (x[i][j + 1] + x[i][j - 1]));

                new_y[i][j] = (-0.5 / denominator) *
                             (2.0 * beta[i][j] * (y[i + 1][j + 1] - y[i - 1][j + 1] - y[i + 1][j - 1] + y[i - 1][j - 1]) -
                              alpha[i][j] * (y[i + 1][j] + y[i - 1][j]) -
                              gamma[i][j] * (y[i][j + 1] + y[i][j - 1]));
            }
        }

        // Update errors
        double maxErrorX = 0.0, maxErrorY = 0.0;
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                maxErrorX = std::max(maxErrorX, std::abs(newx[i][j] - x[i][j]));
                maxErrorY = std::max(maxErrorY, std::abs(newy[i][j] - y[i][j]));
            }
        }
        Er1[t] = maxErrorX;
        Er2[t] = maxErrorY;

        // Neumann BC (left boundary)
        for (int j = 0; j < ny; ++j) newy[nx - 1][j] = newy[nx - 2][j];

        // Update grids
        X = newX;
        Y = newY;

        // Check convergence
        if (maxErrorX < Ermax && maxErrorY < Ermax) {
            std::cout << "Converged at iteration: " << t << std::endl;
            break;
        }

        // Display intermediate results (optional)
        if (show && t % 10 == 0) {
            std::cout << "Iteration " << t << " - Max Error X: " << maxErrorX << ", Max Error Y: " << maxErrorY << std::endl;
        }
    }

    // Final output
    std::cout << "Final Mesh Grid (X and Y):" << std::endl;
    for (int i = 0; i < nx; i += 10) {  // Print every 10th row
        for (int j = 0; j < ny; j += 10) {
            std::cout << "(" << x[i][j] << ", " << y[i][j] << ") ";
        }
        std::cout << std::endl;
    }
}
