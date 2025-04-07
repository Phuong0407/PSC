#include "../include/grid_generation/grid_dimension.h"

void GridDimension::initNumberOfVerticalElements(unsigned int N) {
    NumElemVer = N;
}

void GridDimension::initNumberOfHorizontalElements(unsigned int M1, unsigned int M2, unsigned int M3) {
    NumElemFirstDom = M1;
    NumElemSecondDom = M2;
    NumElemThirdDom = M3;
    NumElemHor = NumElemFirstDom + NumElemSecondDom + NumElemThirdDom;
}

void GridDimension::initNumberOfElements() {
    NumElemTol = NumElemVer * NumElemHor;
    NumNodeTol = 2 * NumElemTol;
}

void GridDimension::initNumberOfNodes() {
    NumNodeHor = NumElemHor + 1;
    NumNodeVer = NumElemVer + 1;
    NumNodeTol = NumNodeHor * NumNodeVer;
    NumElemTol = NumElemVer * NumElemHor;
}

void GridDimension::initGridDimension(unsigned int N, unsigned int M1, unsigned int M2, unsigned int M3) {
    initNumberOfVerticalElements(N);
    initNumberOfHorizontalElements(M1, M2, M3);
    initNumberOfElements();
    initNumberOfNodes();
}