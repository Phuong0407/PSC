#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "../include/laplace_solver.h"
#include "../include/gradient_recovery/gradient_recovery.h"
#include "../include/gradient_recovery/gradient_interpolation.h"

// #include "../scr/laplace_solver.cpp"
// #include "../scr/gradient_recovery/gradient_interpolation.cpp"

using namespace std;

int main() {
    CoordArr Coord;
    NodeSol Solution;
    ElemConn ElemConnData;
    LaplaceSolver A(4.0, 2.0, 3.0, 3.0, 3.0, 10, 5, 10, 5, 1.0);
    A.getCoord(Coord);
    A.solveForNodalVal();
    A.getNodalSolution(Solution);
    A.getElemConn(ElemConnData);

    // std::cout << "MAIN FUNCTION" << std::endl;

    PointDev InterpolatedDev;
    GradientInterpolation B;
    B.initGradientInterpolation(ElemConnData, Coord, Solution);

    std::string filename = "GradDevRec.txt";

    for (std::size_t IndNode = 0; IndNode < Solution.size(); IndNode++) {
        B.interpolateDev(Coord[IndNode].first, Coord[IndNode].second, InterpolatedDev);
        std::cout << "(" << Coord[IndNode].first << ", " << Coord[IndNode].second << ") | zDev = " << InterpolatedDev[0] << "   rDev = " << InterpolatedDev[1] << "   zzDev = " << InterpolatedDev[2] << "   zrDev = " << InterpolatedDev[3] << "   rrDev = " << InterpolatedDev[4] << std::endl;

        std::ofstream outfile(filename , std::ios::app); // Open file in append mode

        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return -1;
        }

        // Write the output to the file
        outfile << "(" << Coord[IndNode].first << ", " << Coord[IndNode].second << ") | zDev = " 
                << InterpolatedDev[0] << "   rDev = " << InterpolatedDev[1]
                << "   zzDev = " << InterpolatedDev[2] << "   zrDev = " << InterpolatedDev[3]
                << "   rrDev = " << InterpolatedDev[4] << std::endl;

        outfile.close(); // Close the file
    }

    return 0;
}