#include "geometry.h"

bool Polygon::isPointInPolygon(const Point& p) const {
    int n = vertices.size();
    if (n < 3) return false; // Not a valid polygon

    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Point& pi = vertices[i];
        const Point& pj = vertices[j];

        if (((pi.y > p.y) != (pj.y > p.y)) &&
            (p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x)) {
            inside = !inside;
        }
    }
    return inside;
}

bool Polygon::isPointOnSegment(const Point& p, const Point& a, const Point& b) const
{
    double cross = (p.y - a.y) * (b.x - a.x) - (p.x - a.x) * (b.y - a.y);
    if (std::abs(cross) > 1e-9) return false;

    double dot = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);
    if (dot < 0) return false;

    double lenSq = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
    return dot <= lenSq;
}

void Polygon::filterTriangles(std::vector<Triangle>& triangles) const {
    triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
                                    [&](const Triangle& t) {
                                        return !isTriangleInPolygon(t);
                                    }),
                    triangles.end());
}

void Polygon::draw() const {
    std::cout << "DRAWING POLYGON WITH VERTICES:\n";
    for (const auto& vertex : vertices) {
        std::cout << "(" << vertex.x << ", " << vertex.y << ")\n";
    }
}

REL_POS_TO_POLYGON Polygon::getPointRelLocToPolygon(const Point& p, Point* edgeStart, Point* edgeEnd) const {
    if (!isValidPolygon()) {
        throw InvalidPolygonException("INVALID POLYGON: A POLYGON MUST HAVE AT LEAST 3 VERTICES.");
    }

    int n = vertices.size();
    int intersections = 0;

    for (int i = 0; i < n; ++i) {
        const Point& a = vertices[i];
        const Point& b = vertices[(i + 1) % n];

        if (isPointOnSegment(p, a, b)) {
            if (edgeStart) *edgeStart = a;
            if (edgeEnd) *edgeEnd = b;
            return LOCATION_BOUNDARY;
        }

        if ((a.y <= p.y && b.y > p.y) || (b.y <= p.y && a.y > p.y)) {
            double xIntersect = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if (xIntersect > p.x) {
                intersections++;
            }
        }
    }

    if (intersections % 2 == 1) {
        return LOCATION_INSIDE;
    } else {
        return LOCATION_OUTSIDE;
    }
}

void Polygon::generateEdgesWithNormals() {
    edges.clear();
    int n = vertices.size();

    for (int i = 0; i < n; ++i) {
        const Point& a = vertices[i];
        const Point& b = vertices[(i + 1) % n];

        double dx = b.x - a.x;
        double dy = b.y - a.y;

        double nx = -dy;
        double ny = dx;

        double length = vectorLength(nx, ny);
        if (length != 0) {
            nx /= length;
            ny /= length;
        }

        Point midpoint((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
        Point testPoint(midpoint.x + nx * 0.1, midpoint.y + ny * 0.1);

        if (getPointRelLocToPolygon(testPoint) == LOCATION_INSIDE) {
            nx = -nx;
            ny = -ny;
        }
        edges.emplace_back(a, b, Vect(nx, ny));
    }
}

void Polygon::dispNormalVector() const {
    if (edges.empty()) {
        throw std::runtime_error("NO EDGES AVAILABLE. PLEASE BUILD EDGES FIRST.");
    }

    std::cout << "EDGE NORMALS:\n";
    for (const auto& edge : edges) {
        std::cout << "EDGE: (" << edge.start.x << ", " << edge.start.y << ") -> ("
                  << edge.end.x << ", " << edge.end.y << "), "
                  << "NORMAL: (" << edge.normal.x << ", " << edge.normal.y << ")\n";
    }
}
