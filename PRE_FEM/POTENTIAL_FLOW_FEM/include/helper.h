#ifndef LAPLACE_SOLVER_H
#define LAPLACE_SOLVER_H

#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <set>
#include <array>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>

struct Point2DStruct {
    double z, r;
    Point2DStruct(double z_, double r_) : z(z_), r(r_) {}
};
using Point2D = struct Point2DStruct;
using Point2DArr = std::vector<Point2D>;

struct Vector2DStruct {
    double nz, nr;
    Vector2DStruct(double nz_, double nr_) : nz(nz_), nr(nr_) {}
};
using Vec2D = struct Vector2DStruct;

using IndOfNode = unsigned int;
using IndOfElem = unsigned int;
using ElemConn = std::vector<IndOfNode>;
using ElemConnArr = std::vector<ElemConn>;

struct EdgeFluxStruct {
    IndOfElem ElemInd;
    IndOfNode FirstNodeInd;
    IndOfNode SecondNodeInd;
    double FirstNodeFlux;
    double SecondNodeFlux;
};
using EdgeFlux = struct EdgeFluxStruct;
using EdgeFluxArr = std::vector<EdgeFlux>;
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
class PhysicalDomainDimension{
private:
    void initVerticalDimension(double L1, double L2);
    void initHorizontalDimension(double H1, double H2, double H3);

public:
    double L1, L2;
    double H1, H2, H3;
    PhysicalDomainDimension() = default;
    void initPhysicalDomainDimension(double L1, double L2, double H1, double H2, double H3);
};
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
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
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
class GridGeneration :
                        private PhysicalDomainDimension,
                        private GridDimension {
private:
    double dy1, dy2;
    double dx1, dx2, dx3;
    double InFlux, OutFlux;
private:
    void initGridStep();

private:
    void genGridCoord();
    void genElemConn();
    // void genIndDirichletBound();
    void genNeumannBoundary();


protected:
    Point2DArr Coord;
    ElemConnArr ElemConnData;
    EdgeFluxArr NeumannBound;
    // IndArr IndDirichletBound;

public:
    GridGeneration() = default;
    void initGridGeneration(
        double L1, double L2,
        double H1, double H2, double H3,
        unsigned int N, unsigned int M1, unsigned int M2, unsigned int M3,
        double InVel
    );

public:
    void dispGridCoord() const;
    void dispElemConn() const;
    void dispIndDirichletBound() const;
    void dispNeumannBound() const;

public:
    void getGridDat(Point2DArr &Coord, ElemConnArr &ElemConnData, EdgeFluxArr &NeumannBound);

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
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */
/* ******************************************************************* */















#endif