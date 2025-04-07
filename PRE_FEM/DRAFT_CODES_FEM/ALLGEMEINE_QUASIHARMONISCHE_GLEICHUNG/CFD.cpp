#include <iostream>
#include <vector>

class CFD {
public:

std::size_t Nx, Ny;
double L, H, L_slot, L_coflow;
double U_slot, U_coflow;
double dx, dy, dt;
double Re;

std::vector<std::vector<double>> u;
std::vector<std::vector<double>> v;
std::vector<double> p;

double rho;
double Re;


public:
CFD() = default;
void construct_first_u();
void init_u_bound();
};

void CFD::init_bound() {
    for(std::size_t i = 0; i < Nx + 1; ++i) {
        if (dx * i <= L_slot) {
            u[0][i] = U_slot;
            u[Ny + 1][i] = - U_slot;
        }
        if (L_slot < dx * i <= L_coflow + L_slot) {
            u[0][i] = U_coflow;
            u[Ny + 1][i] = - U_coflow;
        }
        else {
            u[0][i] = 0.0;
            u[Ny + 1][i] = 0.0;
        }
    }
    for (std::size_t i = 0; i < Ny + 1; ++i) {
        u[i][0] = 0.0;
    }
}



void CFD::construct_first_u() {
    std::vector<std::vector<double>> A;
    std::vector<double> B;
    A.resize(Nx * Ny);
    B.resize(Nx * Ny);
    for(std::size_t i = 0; i < Nx * Ny; ++i)
        A[i].resize(Nx * Ny);

    for(std::size_t i = 1; i < Nx; ++i) {
        for(std::size_t j= 1; j < Ny; ++j) {
            std::size_t ind_ij = i * Nx * Ny + j;
            A[ind_ij][ind_ij] = (1 / dt + 2/(Re * dx * dx) + 2/(Re * dy * dy));
            A[ind_ij][ind_ij + Nx + 1] = - 1/(Re * dx * dx);
            A[ind_ij][ind_ij - Nx - 1] = - 1/(Re * dx * dx);
            A[ind_ij][ind_ij + 1] = - 1/(Re * dy * dy);
            A[ind_ij][ind_ij - 1] = - 1/(Re * dy * dy);
            B[ind_ij] = - u[i][j] * (u[i+1][j] - u[i][j]) / dx - v[i][j] * (u[i][j+1] - u[i][j]) / dy + 1 / dt * u[i][j];
        }
    }
}