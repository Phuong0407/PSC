#include "../include/grid_generation.h"

int main() {
    GridGeneration GridGenerationHandler;
    GridGenerationHandler.initPhysicalDomainDimension(4.0, 2.0, 3.0, 3.0, 3.0);
    GridGenerationHandler.initGridDomainDimension(2, 2, 2, 2);
    GridGenerationHandler.buildGridGenerationData();
    GridGenerationHandler.buildBoundaryVelocity(1.0);
    GridGenerationHandler.buildBoundaryCondition();
    return 0;
}