/**
 * The flow domain is define as follow
 * 6------------------------------------5
 * |                                    |
 * |                                    |
 * |                                    |
 * |                                    |
 * |                  3-----------------4
 * |                 /
 * |                /
 * 1---------------2
 * The domain is symmetric about the axis 1-2 which has y = 0
 * The edge 1 - 2 - 3 - 4 has zero flux
 * The edge 5 - 6 has wall condition
 * The edge 1 - 6 has inlet conditon, velocity U_in
 * The edge 4 - 5 has outlet conditon, velocity U_out * 
 */


#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <vector>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <utility>

constexpr unsigned int OUTSIDE = 0;
constexpr unsigned int INSIDE = 1;
constexpr unsigned int INTERSECTING = 2;
constexpr unsigned int ON_BOUNDARY = 3;


class BoundaryEdge {
private:
    std::size_t first_node_ind, second_node_ind;
    double x1, x2, y1, y2;
    double nx, ny;
    double flux;

public:
    BoundaryEdge() = default;
    void initNodeIndices(std::size_t first_node_ind, std::size_t second_node_ind) {
        this->first_node_ind = first_node_ind;
        this->second_node_ind = second_node_ind;
    }
    void initEdgeCoordinates(double x1, double y1, double x2, double y2) {
        this->x1 = x1;
        this->y1 = y1;
        this->x2 = x2;
        this->y2 = y2;
    }
    void initEdgeOutwardNormalVector(double nx, double ny) {
        this->nx = nx;
        this->ny = ny;
    }
    void initFlux(double flux) {this->flux = flux;}
    double get_x1() const {return x1;}
    double get_y1() const {return y1;}
    double get_x2() const {return x2;}
    double get_y2() const {return y2;}
    double get_nx() const {return nx;}
    double get_ny() const {return ny;}
    double get_flux() const {return flux;}
    std::size_t get_first_node_ind() const {return first_node_ind;}
    std::size_t get_second_node_ind() const {return second_node_ind;}
};

class PolygonHelper {
public:
    static void compute_outward_norm_vector(BoundaryEdge& edge);
};

void PolygonHelper::compute_outward_norm_vector(BoundaryEdge& edge) {
    const double x1 = edge.get_x1(), y1 = edge.get_y1();
    const double x2 = edge.get_x2(), y2 = edge.get_y2();
    double dir_x = x2 - x1, dir_y = y2 - y1;
    double length = std::sqrt(dir_x * dir_x + dir_y * dir_y);
    double norm_x = -dir_y/length, norm_y = dir_x/length;
    double cross_product = dir_x * norm_y - dir_y * norm_x;
    if (cross_product > 0) {
        norm_x = - norm_x;
        norm_y = - norm_y;
    }
    edge.initEdgeOutwardNormalVector(norm_x, norm_y);
}


class Geometry {
private:
    double L1, L2, H1, H2, H3;
    double dx1, dx2, dx3, dy1, dy2;
    std::size_t N, M, M1, M2, M3;
    std::vector<std::vector<double>> x, y;
    std::vector<double> x_leftbound, y_leftbound;
    std::vector<double> x_rightbound, y_rightbound;
    std::vector<double> x_bottombound, y_bottombound;
    std::vector<double> x_topbound, y_topbound;
    std::vector<std::vector<std::size_t>> ElemConnData;
    std::vector<BoundaryEdge> NeumannBoundaryConditions;

public:
    Geometry() = default;
    void initPhysicalDimensions(double H1, double H2, double H3, double L1, double L2);
    void initComputationalDimensions(std::size_t M1, std::size_t M2, std::size_t M3, std::size_t N);
    void initGridStep();
    void initCoordsArrSize();
    void generateAlgebraicGrid();
    void generateEllipicGrid();
    void generateBoundaryNode();
    void transfiniteInterpolation();
    void generateStretchingGrid(double alpha, double eta1);
    void generateGridConnection();
    void generateNeumannBoundaryConditions(double Flux1, double Flux2, double Flux3, double Flux4);
    std::vector<BoundaryEdge> getNeumannBoundaryConditions() const;

public:
    std::vector<double> get_x() const;
    std::vector<double> get_y() const;
    const std::vector<std::vector<std::size_t>>& getElementConnectionData() const {return ElemConnData;}
    std::size_t getTotalHorizontalElements() const {return M;}
    std::size_t getTotalVerticalElements() const {return N;}
};

void Geometry::initPhysicalDimensions(double H1, double H2, double H3, double L1, double L2) {
    this->H1 = H1;
    this->H2 = H2;
    this->H3 = H3;
    this->L1 = L1;
    this->L2 = L2;
}

void Geometry::initComputationalDimensions(std::size_t M1, std::size_t M2, std::size_t M3, std::size_t N) {
    this->M1 = M1;
    this->M2 = M2;
    this->M3 = M3;
    this->N = N;
    this->M = M1 + M2 + M3;
}

void Geometry::initGridStep() {
    dx1 = H1 / M1;
    dx2 = H2 / M2;
    dx3 = H3 / M3;
    dy1 = L1 / N;
    dy2 = L2 / N;
}

void Geometry::initCoordsArrSize() {
    x.resize(N + 1, std::vector<double>(M + 1));
    y.resize(N + 1, std::vector<double>(M + 1));
    x_leftbound.resize(N + 1);
    y_leftbound.resize(N + 1);
    x_rightbound.resize(N + 1);
    y_rightbound.resize(N + 1);
    x_bottombound.resize(M + 1);
    y_bottombound.resize(M + 1);
    x_topbound.resize(M + 1);
    y_topbound.resize(M + 1);
}

void Geometry::generateAlgebraicGrid() {
    for (std::size_t i = 0; i <= N; ++i) {
        for (std::size_t j = 0; j <= M1; ++j)
            x[i][j] = dx1 * j, y[i][j] = dy1 * i;

        double a = ((dy2 - dy1) * i + L1 - L2) / H2;
        double b = dy1 * i - a * H1;
        for (std::size_t j = M1 + 1; j <= M1 + M2; ++j) {
            double x_1 = H1 + dx2 * (j - M1);
            double y_1 = a * x_1 + b;
            x[i][j] = x_1;
            y[i][j] = y_1;
        }

        for (std::size_t j = M1 + M2 + 1; j <= M; ++j) {
            x[i][j] = H1 + H2 + dx3 * (j - M1 - M2);
            y[i][j] = L1 - L2 + dy2 * i;
        }
    }
}

void Geometry::generateEllipicGrid() {
    if (x.empty() || y.empty())
        throw std::runtime_error("ERROR: ALGEBRAIC-GENERATING COORDINATES ARRAY IS EMPTY. GENERATE IT FIRST BY USING generateAlgebraicGrid()!");

    std::vector<std::vector<double>> new_x(N + 1, std::vector<double>(M + 1, 0.0));
    std::vector<std::vector<double>> new_y(N + 1, std::vector<double>(M + 1, 0.0));
    std::vector<std::vector<double>> alpha(N + 1, std::vector<double>(M + 1, 0.0));
    std::vector<std::vector<double>> beta(N + 1, std::vector<double>(M + 1, 0.0));
    std::vector<std::vector<double>> gamma(N + 1, std::vector<double>(M + 1, 0.0));
    new_x = x, new_y = y;

    double min_err = 1e-6;
    std::size_t max_iter = 10000;
    std::vector<double> error_1(max_iter, 0.0);
    std::vector<double> error_2(max_iter, 0.0);

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < M; ++j) {
                double Dxi_x = (x[i][j + 1] - x[i][j - 1]);
                double Dxi_y = (y[i][j + 1] - y[i][j - 1]);
                double Deta_x = (x[i + 1][j] - x[i - 1][j]);
                double Deta_y = (y[i + 1][j] - y[i - 1][j]);
                alpha[i][j] = 0.25 * (std::pow(Deta_x, 2.0) + std::pow(Deta_y, 2.0));
                gamma[i][j] = 0.25 * (std::pow(Dxi_x, 2.0) + std::pow(Dxi_y, 2.0));
                beta[i][j] = 0.0625 * (Deta_x * Dxi_x + Deta_y * Dxi_y);
                double denominator = alpha[i][j] + gamma[i][j] + 1e-9;
                new_x[i][j] = 0.5 / denominator * (-2.0 * beta[i][j] * (x[i + 1][j + 1] - x[i + 1][j - 1] - x[i - 1][j + 1] + x[i - 1][j - 1]) +
                              alpha[i][j] * (x[i][j + 1] + x[i][j - 1]) + gamma[i][j] * (x[i + 1][j] + x[i - 1][j]));
                new_y[i][j] = 0.5 / denominator * (-2.0 * beta[i][j] * (y[i + 1][j + 1] - y[i + 1][j - 1] - y[i - 1][j + 1] + y[i - 1][j - 1]) +
                              alpha[i][j] * (y[i][j + 1] + y[i][j - 1]) + gamma[i][j] * (y[i + 1][j] + y[i - 1][j]));
            }
        }
        double max_error_x = 0.0, max_error_y = 0.0;
        for (std::size_t i = 1; i < N; ++i) {
            for (std::size_t j = 1; j < M; ++j) {
                max_error_x = std::max(max_error_x, std::abs(new_x[i][j] - x[i][j]));
                max_error_y = std::max(max_error_y, std::abs(new_y[i][j] - y[i][j]));
            }
        }
        error_1[iter] = max_error_x;
        error_2[iter] = max_error_y;
        x = new_x;
        y = new_y;
        if (error_1[iter] < min_err && error_2[iter] < min_err)
            break;
    }
}

void Geometry::transfiniteInterpolation() {
    double deta = 1.0 / N;
    double dxi = 1.0 / M;
    for (std::size_t i = 0; i < N + 1; ++i) {
        double eta = deta * i;
        for (std::size_t j = 0; j < M + 1; ++j) {
            double xi = dxi * j;
            x[i][j] = (1 - xi) * x_leftbound[i] + xi * x_rightbound[i]
                     + (1 - eta) * x_bottombound[j] + eta * x_topbound[j]
                     - (1- xi) * (1 - eta) * x_bottombound[0] - (1 - xi) * eta * x_topbound[0]
                     - (1 - eta) * xi * x_bottombound[M] - xi * eta * x_topbound[M];
            y[i][j] = (1 - xi) * y_leftbound[i] + xi * y_rightbound[i] +
                      (1 - eta) * y_bottombound[j] + eta * y_topbound[j] -
                      (1- xi) * (1 - eta) * y_bottombound[0] - (1 - xi) * eta * y_topbound[0] -
                      (1 - eta) * xi * y_bottombound[M] - xi * eta * y_topbound[M];
        }        
    }
}

void Geometry::generateBoundaryNode() {
    for (std::size_t i = 0; i <= M1; ++i) {
        double x = H1/M1 * i;
        x_bottombound[i] = x;
        y_bottombound[i] = 0.0;
        x_topbound[i] = x;
        y_topbound[i] = L1;
    }
    double a = (L1 - L2) / H2;
    for (std::size_t i = M1 + 1; i <= M1 + M2; ++i) {
        double x = H1 + dx2 * (i - M1);
        double y = a * dx2 * (i - M1);
        x_bottombound[i] = x;
        y_bottombound[i] = y;
        x_topbound[i] = x;
        y_topbound[i] = L1;
    }
    for (std::size_t i = M1 + M2 + 1; i <= M; ++i) {
        double x = H1 + H2 + dx3 * (i - M1 - M2);
        x_bottombound[i] = x;
        y_bottombound[i] = L1 - L2;
        x_topbound[i] = x;
        y_topbound[i] = L1;
    }
    for (std::size_t i = 0; i <= N; ++i) {
        x_leftbound[i] = 0.0;
        x_rightbound[i] = H1 + H2 + H3;
        y_leftbound[i] = L1/N * i;
        y_rightbound[i] = L1 - L2 + L2/N * i;
    }
}

void Geometry::generateStretchingGrid(double alpha, double eta1) {
    double deta = 1.0 / N;
    double dxi = 1.0 / M;
    for (std::size_t i = 0; i < N + 1; ++i) {
        double eta = i * deta;
        for (std::size_t j = 0; j <= M; ++j) {
            double xi = j * dxi;
            x[i][j] = x_bottombound[j];
            if (eta <= eta1) {
                y[i][j] = (y_topbound[j] - y_bottombound[j]) * eta1 * (std::exp(alpha * eta / eta1) - 1.0) / (std::exp(alpha) - 1) + y_bottombound[j];
            }
            else {
                y[i][j] = (y_topbound[j] - y_bottombound[j]) * (1.0 - (1.0 - eta1) * (std::exp(alpha * (1 - eta) / (1 - eta1)) - 1.0) / (std::exp(alpha) - 1)) + y_bottombound[j];
            }
        }
    }
}


void Geometry::generateGridConnection() {
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < M; ++j) {
            unsigned int ind_elem = i * (M + 1) + j;
            unsigned int first_node_ind = ind_elem;
            unsigned int second_node_ind = first_node_ind + 1;
            unsigned int third_node_ind = second_node_ind + M + 1;
            unsigned int fourth_node_ind = first_node_ind + M + 1;
            ElemConnData.push_back({first_node_ind, second_node_ind, third_node_ind});
            ElemConnData.push_back({first_node_ind, third_node_ind, fourth_node_ind});
        }
    }
}

std::vector<double> Geometry::get_x() const {
    std::vector<double> x_coords;
    for (std::size_t i = 0; i < N + 1; ++i) {
        for (std::size_t j = 0; j < M + 1; ++j) {
            x_coords.push_back(x[i][j]);
        }
    }
    return x_coords;
}
std::vector<double> Geometry::get_y() const {
    std::vector<double> y_coords;
    for (std::size_t i = 0; i < N + 1; ++i) {
        for (std::size_t j = 0; j < M + 1; ++j) {
            y_coords.push_back(y[i][j]);
        }
    }
    return y_coords;
}

void Geometry::generateNeumannBoundaryConditions(double Flux1, double Flux2, double Flux3, double Flux4) {
    NeumannBoundaryConditions.clear();
    for (std::size_t i = 0; i < M; ++i) {
        std::size_t first_node_ind = i, second_node_ind = first_node_ind + 1;
        double x1 = x[0][i], y1 = y[0][i];
        double x2 = x[0][i + 1], y2 = y[0][i + 1];
        BoundaryEdge edge;
        edge.initNodeIndices(first_node_ind, second_node_ind);
        edge.initEdgeCoordinates(x1, y1, x2, y2);
        edge.initFlux(Flux1);
        PolygonHelper::compute_outward_norm_vector(edge);
        NeumannBoundaryConditions.push_back(edge);
    }
    for (std::size_t i = 0; i < N; ++i) {
        std::size_t first_node_ind = i * (M + 1) + M, second_node_ind = first_node_ind + M + 1;
        double x1 = x[i][M], y1 = y[i][M];
        double x2 = x[i + 1][M], y2 = y[i + 1][M];
        BoundaryEdge edge;
        edge.initNodeIndices(first_node_ind, second_node_ind);
        edge.initEdgeCoordinates(x1, y1, x2, y2);
        edge.initFlux(Flux2);
        PolygonHelper::compute_outward_norm_vector(edge);
        NeumannBoundaryConditions.push_back(edge);
    }
    for (int i = M; i > 0; --i) {
        std::size_t first_node_ind = N * (M + 1) + i, second_node_ind = first_node_ind - 1;
        double x1 = x[N][i], y1 = y[N][i];
        double x2 = x[N][i - 1], y2 = y[N][i - 1];
        BoundaryEdge edge;
        edge.initNodeIndices(first_node_ind, second_node_ind);
        edge.initEdgeCoordinates(x1, y1, x2, y2);
        edge.initFlux(Flux3);
        PolygonHelper::compute_outward_norm_vector(edge);
        NeumannBoundaryConditions.push_back(edge);
    }
    for (int i = N; i > 0; --i) {
        std::size_t first_node_ind = i * (M + 1), second_node_ind = first_node_ind - M - 1;
        double x1 = x[i][0], y1 = y[i][0];
        double x2 = x[i - 1][0], y2 = y[i - 1][0];
        BoundaryEdge edge;
        edge.initNodeIndices(first_node_ind, second_node_ind);
        edge.initEdgeCoordinates(x1, y1, x2, y2);
        edge.initFlux(Flux4);
        PolygonHelper::compute_outward_norm_vector(edge);
        NeumannBoundaryConditions.push_back(edge);
    }
}

std::vector<BoundaryEdge> Geometry::getNeumannBoundaryConditions() const {
    return NeumannBoundaryConditions;
}


#endif