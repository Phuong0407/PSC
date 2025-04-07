#include "../include/glo_stiff_system/dirichlet_bound.h"


void DirichletBound::genDirichletBound(
    const double L1, const double L2,
    const double InVel, const double OutVel,
    const unsigned int NumNodeHor, const unsigned int NumNodeVer, 
    const CoordArr &Coord, const IndArr &IndDirichletBound, ColVect &ValDirichletBound
)
{
    allocateColVect(ValDirichletBound, IndDirichletBound.size());

    for (std::size_t i = 0; i < IndDirichletBound.size(); ++i) {
        unsigned int IndNode = IndDirichletBound[i];
        unsigned int IndRow = IndNode / NumNodeHor;
        unsigned int IndCol = IndNode % NumNodeHor;
        
        if (IndRow == 0) {
            ValDirichletBound[i] = 0.0;
            continue;
        }
        if (IndCol == 0) {
            double r = Coord[IndNode].second;
            ValDirichletBound[i] = 0.5 * InVel * r * r;
            continue;
        }
        if (IndCol == NumNodeHor - 1) {
            double r = Coord[IndNode].second;
            
            ValDirichletBound[i] = 0.5 * OutVel * (std::pow(r, 2.0) - std::pow(L2, 2.0));
            continue;
        }
        if (IndRow == NumNodeVer - 1) {
            ValDirichletBound[i] = 0.5 * InVel * std::pow(L1, 2.0);
            continue;
        }
    }
}

void DirichletBound::dispDirichletBound(const std::vector<unsigned int> &IndDirichletBound, const std::vector<double> &ValDirichletBound) const {
    for (std::size_t i = 0; i < ValDirichletBound.size(); ++i) {
        std::cout << IndDirichletBound[i] << " Dirichlet Bound psi = " << ValDirichletBound[i] << std::endl;
    }
}
