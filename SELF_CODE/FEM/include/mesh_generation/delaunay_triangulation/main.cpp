#include "geometry_helper.hpp"

#include <iostream>

using namespace mesh_generation;

int main() {
    std::cout << "Adding simplices...\n";

    VertexCoordinate<2> a;
    
    // Add first triangle: (1, 2, 3) with coordinates
    a.addVertex(10.0, 20.0);

    // Add second triangle: (2, 3, 4) with coordinates
    a.addVertex(10.02, 20.0);
    a.addVertex(20.3, 20.0);
    a.addVertex(10.4, 20.0);
    a.addVertex(2.2, 20.0);

    std::cout << "Stored simplices and coordinates:\n";
    a.printCoordinates();

    return 0;
}