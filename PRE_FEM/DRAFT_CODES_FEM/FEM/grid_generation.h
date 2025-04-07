#ifndef GRID_GENERATION_H
#define GRID_GENERATION_H

#include "geometry.cc"
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

class BowyerWatson {
private:
    std::vector<Triangle> triangles;

    Triangle getSuperTriangle(const std::vector<Point>& points) const;
    bool isPointInCircumcircle(const Point& p, const Triangle& t) const;
    std::vector<std::pair<Point, Point>> getPolygonalHoleEdges(const std::vector<Triangle>& badTriangles) const;
    void removeBadTriangles(const std::vector<Triangle>& badTriangles);
    void removeDuplicateEdges(std::vector<std::pair<Point, Point>>& edges) const;

public:
    void triangulate(const std::vector<Point>& points, const Polygon& polygon);
    const std::vector<Triangle>& getTriangles() const { return triangles; }
    void printTriangles() const;
};

#endif