#include <mesh_generation/delaunay_triangulation/delaunator.hpp>

#include <geometry/geometry.hpp>

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <random>

namespace geometry {
    using Point2D = Point<2>;
}

std::vector<geometry::Point2D> generatePoints() {
    std::vector<geometry::Point2D> points;
    points.reserve(20);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0); 

    for (int i = 0; i < 20; ++i) {
        points.emplace_back(dist(gen), dist(gen));
    }

    return points;
}
void writeGnuplotScript(const std::vector<geometry::Point2D>& coords, const std::vector<std::size_t>& triangles) {
    std::ofstream gnuplot_script("triangulation.gnu");

    if (!gnuplot_script.is_open()) {
        std::cerr << "Error: Could not create GNUplot script file.\n";
        return;
    }

    gnuplot_script << "set title 'Delaunay Triangulation'\n";
    gnuplot_script << "set xrange [-10:10]\n";
    gnuplot_script << "set yrange [-10:10]\n";
    gnuplot_script << "set size ratio -1\n";
    gnuplot_script << "plot '-' with lines lc rgb 'blue' title 'Triangulation'\n";

    for (std::size_t i = 0; i + 2 < triangles.size(); i += 3) {
        if (triangles[i] >= coords.size() || triangles[i + 1] >= coords.size() || triangles[i + 2] >= coords.size()) {
            std::cerr << "Error: Triangle index out of bounds.\n";
            continue;
        }

        const auto& p0 = coords[triangles[i]];
        const auto& p1 = coords[triangles[i + 1]];
        const auto& p2 = coords[triangles[i + 2]];

        gnuplot_script << p0.coordinate(0) << " " << p0.coordinate(1) << "\n"
                       << p1.coordinate(0) << " " << p1.coordinate(1) << "\n\n";
        gnuplot_script << p1.coordinate(0) << " " << p1.coordinate(1) << "\n"
                       << p2.coordinate(0) << " " << p2.coordinate(1) << "\n\n";
        gnuplot_script << p2.coordinate(0) << " " << p2.coordinate(1) << "\n"
                       << p0.coordinate(0) << " " << p0.coordinate(1) << "\n\n";
    }

    gnuplot_script << "e\n";
    gnuplot_script.close();
}

void runGnuplot() {
    int status = system("gnuplot -p triangulation.gnu");
    if (status == 0) {
        std::cout << "GNUplot executed successfully.\n";
    } else {
        std::cerr << "Error: Failed to execute GNUplot.\n";
    }
}

int main() {
    std::vector<geometry::Point2D> coords = generatePoints();

    mesh_generation::Delaunator d(coords);
    std::vector<std::size_t> triangles = d.getTriangle();

    writeGnuplotScript(coords, triangles);
    runGnuplot();

    return 0;
}