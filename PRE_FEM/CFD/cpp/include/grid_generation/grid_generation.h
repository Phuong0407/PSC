#ifndef GRID_GENERATION_H
#define GRID_GENERATION_H

#include "../include/common_type_and_method.h"
#include "../include/grid_generation/physical_domain_dimension.h"
#include "../include/grid_generation/grid_dimension.h"

#include <set>
#include <array>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>

class GridGeneration :      private PhysicalDomainDimension,
                            private GridDimension
{
private:
    double dy1, dy2;
    double dx1, dx2, dx3;

private:
    void initGridStep();

private:
    void genGridCoord();
    void genElemConn();
    void genIndDirichletBound();

protected:
    CoordArr Coord;
    ElemConn ElemConnData;
    IndArr IndDirichletBound;

public:
    GridGeneration() = default;
    void initGridGeneration(
        double L1, double L2,
        double H1, double H2, double H3,
        unsigned int N, 
        unsigned int M1, unsigned int M2, unsigned int M3
    );

public:
    void dispGridCoord() const;
    void dispElemConn() const;
    void dispIndDirichletBound() const;

public:
    void getGridDat(CoordArr &Coord, ElemConn &ElemConnData, IndArr &IndDirichletBound);

public:
    inline unsigned int getTolNumNode() const {
        return NumNodeTol;
    }

    inline unsigned int getNumNodeHor() const {
        return NumNodeHor;
    }

    inline unsigned int getNumNodeVer() const {
        return NumNodeVer;
    }
};

#endif