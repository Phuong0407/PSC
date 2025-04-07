#ifndef DIRICHLET_BOUND_H
#define DIRICHLET_BOUND_H

#include "../include/common_type_and_method.h"

#include <cmath>
#include <vector>
#include <iostream>

class DirichletBound {
protected:
    DirichletBound() = default;
    void genDirichletBound(
        const double L1, const double L2,
        const double InVel, const double OutVel,
        const unsigned int NumNodeHor, const unsigned int NumNodeVer, 
        const CoordArr &Coord, const IndArr &IndDirichletBound, ColVect &ValDirichletBound
    );

public:
    void dispDirichletBound(
                            const std::vector<unsigned int> &IndDirichletBound,
                            const std::vector<double> &ValDirichletBound
                        ) const;
};


#endif