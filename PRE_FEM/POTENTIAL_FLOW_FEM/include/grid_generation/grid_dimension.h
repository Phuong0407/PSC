#ifndef GRID_DIMENSION_H
#define GRID_DIMENSION_H

#include <vector>

class GridDimension{
private:
    void initNumberOfVerticalElements(unsigned int N);
    void initNumberOfHorizontalElements(unsigned int M1, unsigned int M2, unsigned int M3);
    void initNumberOfElements();
    void initNumberOfNodes();

protected:
    unsigned int NumElemFirstDom;
    unsigned int NumElemSecondDom;
    unsigned int NumElemThirdDom;

    unsigned int NumElemHor;
    unsigned int NumElemVer;
    unsigned int NumElemTol;
    unsigned int NumElemTri;
    
    unsigned int NumNodeTol;
    unsigned int NumNodeHor;
    unsigned int NumNodeVer;

public:
    GridDimension() = default;
    void initGridDimension(unsigned int N, unsigned int M1, unsigned int M2, unsigned int M3);
};

#endif