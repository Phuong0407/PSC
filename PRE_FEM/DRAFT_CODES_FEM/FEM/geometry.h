/*
 * This is the geometric configuration for flow field domain
 * An axisymmetric configuration around z-axis: edge 1-2
 * The original project use 6 points:
 * 6---------------------------------------5
 * |                                       |
 * |                                       |
 * |                                       |
 * |                                       |
 * |                    3------------------4
 * |                   /
 * |                  /
 * 1_________________2
 * The whole problem using neumann boundary condition:
 * 1-2-3-4 has no flux
 * 4-5 has outflux
 * 5-6 has wall condition
 * 6-1 has influx
 * COUNTER-CLOCKWISE orientation for all polygon
 */
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <stdexcept>
#include <string>
#include <cmath>

using Size = unsigned int;

// constexpr int LOCATION_INSIDE = 1;
// constexpr int LOCATION_OUTSIDE = 2;
// constexpr int LOCATION_BOUNDARY = 3;

// using REL_POS_TO_POLYGON = int;

class Point2D {
private:
    double x, y;

public:
    Point2D() = default;
    Point2D(double x, double y) : x(x), y(y) {}
};

struct Point {
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}

    bool operator==(const Point& other) const {
        const double EPSILON = 1e-9;
        return (std::fabs(x - other.x) < EPSILON && std::fabs(y - other.y) < EPSILON);
    }
};

struct Vect {
    double x, y;
    Vect() : x(0), y(0) {}
    Vect(double x, double y) : x(x), y(y) {}
};

inline double vectorLength(double x, double y) {
    return std::sqrt(x * x + y * y);
}

struct Edge {
    Point start;
    Point end;
    Vect normal;
    Edge() : start(), end(), normal() {}
    Edge(const Point& start_, const Point& end_, const Vect& normal_)
        : start(start_), end(end_), normal(normal_) {}
};

struct Triangle {
    Point p1, p2, p3;

    Triangle(const Point& p1, const Point& p2, const Point& p3) : p1(p1), p2(p2), p3(p3) {}
    bool operator==(const Triangle& other) const {
        return (p1 == other.p1 && p2 == other.p2 && p3 == other.p3) ||
               (p1 == other.p2 && p2 == other.p3 && p3 == other.p1) ||
               (p1 == other.p3 && p2 == other.p1 && p3 == other.p2);
    }
    bool containsVertex(const Point& vertex) const {
        return p1 == vertex || p2 == vertex || p3 == vertex;
    }
};


class InvalidPolygonException : public std::runtime_error {
public:
    explicit InvalidPolygonException(const std::string& message)
        : std::runtime_error(message) {}
};


class Polygon {
private:
    std::vector<Point> vertices;
    std::vector<Edge> edges;
private:
    bool isPointOnSegment(const Point& p, const Point& a, const Point& b) const;
    inline bool isValidPolygon() const { return vertices.size() >= 3; }
    bool isPointInPolygon(const Point& p) const;
    bool isTriangleInPolygon(const Triangle& triangle) const {
        return isPointInPolygon(triangle.p1) &&
               isPointInPolygon(triangle.p2) &&
               isPointInPolygon(triangle.p3);
    }
public:
    REL_POS_TO_POLYGON getPointRelLocToPolygon(const Point& p, Point* edgeStart = nullptr, Point* edgeEnd = nullptr) const;
    void generateEdgesWithNormals();

public:
    inline void add_vertex(double x, double y) { vertices.emplace_back(x, y); }

public:
    void filterTriangles(std::vector<Triangle>& triangles) const;
    void draw() const;
    void dispNormalVector() const;
};

#endif