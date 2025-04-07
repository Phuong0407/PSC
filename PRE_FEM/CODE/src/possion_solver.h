#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "matrix.h"

class NodeCoordinates{
    public:
    double x = 0.0, y = 0.0;
    NodeCoordinates(double x, double y) : x(x), y(y) {}
    ~NodeCoordinates() {}
};

class GlobalStiffnessMatrix{
    public:
    unsigned int nElem; // number of elements
    unsigned int nNOde; // number of nodes
    unsigned int hNode; // number of horizontal nodes
    unsigned int vNode; // number of vertical nodes

    NodeCoordinates* nCoordinatesMatrix; // Node Coordinates Matrix
    

};