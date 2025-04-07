#include "boundary_element_generation.cpp"
// #include "boundary_integral_method.cpp"
#include "boundary_integral_method.h"
#include "dense_system_solver.cpp"

int main() {
    ElemArr ElemConn;
    Point2DArr NodeCoord;
    DoubleArr SolidAngle;
    NormVecArr OutwardNormVec;
    // DoubleArr DevPhi;
    FluxArr DevPhi2;
    
    // BoundaryElementGeneration BoundaryElementGenerationHandler;
    // BoundaryElementGenerationHandler.initBoundaryElementGeneration(4.0, 2.0, 3.0, 3.0, 3.0, 200, 200, 200, 600, 200, 200);
    // BoundaryElementGenerationHandler.initVelocity(1.0);

    // BoundaryElementGenerationHandler.dispNodeCoord();
    // BoundaryElementGenerationHandler.dispOutwardNormVec();
    // BoundaryElementGenerationHandler.dispSolidAngle();
    // BoundaryElementGenerationHandler.dispElemConn();
    // BoundaryElementGenerationHandler.dispFlux();

    // BoundaryElementGenerationHandler.getElemConn(ElemConn);
    // BoundaryElementGenerationHandler.getNodeCoord(NodeCoord);
    // BoundaryElementGenerationHandler.getOutwardNormVec(OutwardNormVec);
    // BoundaryElementGenerationHandler.getSolidAngle(SolidAngle);
    // BoundaryElementGenerationHandler.getFlux(DevPhi);

    NodeCoord.emplace_back(10.0, 0.0);
    NodeCoord.emplace_back(10.0, 10.0);
    NodeCoord.emplace_back(5.0, 10.0);
    NodeCoord.emplace_back(0.0, 10.0);
    NodeCoord.emplace_back(0.0, 0.0);

    ElemConn.emplace_back(0,1);
    ElemConn.emplace_back(1,2);
    ElemConn.emplace_back(2,3);
    ElemConn.emplace_back(3,4);
    // ElemConn.emplace_back(4,0);

    SolidAngle.emplace_back(M_PI);
    SolidAngle.emplace_back(M_PI);
    SolidAngle.emplace_back(2 * M_PI);
    SolidAngle.emplace_back(M_PI);
    SolidAngle.emplace_back(M_PI);

    OutwardNormVec.emplace_back(1.0,0.0);
    OutwardNormVec.emplace_back(0.0,1.0);
    OutwardNormVec.emplace_back(0.0,1.0);
    OutwardNormVec.emplace_back(-1.0,0.0);
    // OutwardNormVec.emplace_back(0.0,-1.0);

    // DevPhi.push_back(1.0);
    // DevPhi.push_back(1.0);
    // DevPhi.push_back(0.0);
    // DevPhi.push_back(1.0);
    // DevPhi.push_back(1.0);

    DevPhi2.resize(ElemConn.size());
    DevPhi2[0][0] = 1.0;
    DevPhi2[0][1] = 1.0;

    DevPhi2[1][0] = 0.0;
    DevPhi2[1][1] = 0.0;

    DevPhi2[2][0] = 0.0;
    DevPhi2[2][1] = 0.0;

    DevPhi2[3][1] = -1.0;
    DevPhi2[3][1] = -1.0;

    // DevPhi2[4][0] = 1.0;
    // DevPhi2[4][1] = 1.0;
    
    // for(const auto &node : NodeCoord) {
    //     std::cout << "(" << node.z << ", " << node.r << ")" << std::endl;
    // }

    Mat A;
    VecCol B;
    buildSystemOfEquation(NodeCoord, ElemConn, OutwardNormVec, SolidAngle, DevPhi2, A, B);

    std::cout << "Matrix A:" << std::endl;
    for (const auto &row : A) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Vector B:" << std::endl;
    for (const auto &elem : B) {
        std::cout << elem << std::endl;
    }

    VecCol Solution;
    DenseSystemSolver DenseSystemSolverHandler;
    DenseSystemSolverHandler.solveDenseMatrixSystem(A, B, Solution);
    std::cout << "Solution:" << std::endl;
    for (std::size_t i = 0; i < Solution.size(); ++i) {
        std::cout << "(" << NodeCoord[i].z << ", " << NodeCoord[i].r << ") = " << Solution[i] << std::endl;
    }

    // // const Point2D &TargetPoint = std::make_pair(1, 0);
    // // const Point2D &FirstRefPoint = std::make_pair(0, 0);
    // // const Point2D &SecondRefPoint = std::make_pair(1, 0);
    // // const Vec2D &NormVec = std::make_pair(1, 0);

    // // ElemGreenFunc Result = calcWeightOfPotentialAndDevAtPointWithRespectedToAReferencePoint({0, 0}, {0, 0}, {1, 0}, {1,0});
    // // std::cout << Result[0] << " " << Result[1] << " " << Result[2] << " " << Result[3] << std::endl;

    // std::cout << std::comp_ellint_2(1.0) << std::endl;
    // double K11, K12, K21, K22;
    // calcSingularK1K2Integral(10.0, 0.0, 0.0, 0.0, 10.0, 0.0, M_PI, K11, K12, K21, K22);
    // std::cout << K11 << " " << K12 << " " << K21 << " " << K22 << std::endl;FirstPointInd

    // std::cout << std::comp_ellint_1(std::sqrt(0.9758065490037983));

    return 0;
}