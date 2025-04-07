#include "grid_generation.h"

bool BowyerWatson::isPointInCircumcircle(const Point& p, const Triangle& t) const {
    double ax = t.p1.x - p.x;
    double ay = t.p1.y - p.y;
    double bx = t.p2.x - p.x;
    double by = t.p2.y - p.y;
    double cx = t.p3.x - p.x;
    double cy = t.p3.y - p.y;

    double det = (ax * ax + ay * ay) * (bx * cy - cx * by) -
                    (bx * bx + by * by) * (ax * cy - cx * ay) +
                    (cx * cx + cy * cy) * (ax * by - bx * ay);

    return det > 0;
}

Triangle BowyerWatson::getSuperTriangle(const std::vector<Point>& points) const {
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();

    for (const auto& p : points) {
        minX = std::min(minX, p.x);
        minY = std::min(minY, p.y);
        maxX = std::max(maxX, p.x);
        maxY = std::max(maxY, p.y);
    }

    double dx = maxX - minX;
    double dy = maxY - minY;
    double deltaMax = std::max(dx, dy);
    double midX = (minX + maxX) / 2;
    double midY = (minY + maxY) / 2;

    Point p1(midX - 2 * deltaMax, midY - deltaMax);
    Point p2(midX + 2 * deltaMax, midY - deltaMax);
    Point p3(midX, midY + 2 * deltaMax);

    return Triangle(p1, p2, p3);
}

std::vector<std::pair<Point, Point>> BowyerWatson::getPolygonalHoleEdges(const std::vector<Triangle>& badTriangles) const {
    std::vector<std::pair<Point, Point>> edges;

    for (const auto& triangle : badTriangles) {
        edges.emplace_back(triangle.p1, triangle.p2);
        edges.emplace_back(triangle.p2, triangle.p3);
        edges.emplace_back(triangle.p3, triangle.p1);
    }

    return edges;
}

void BowyerWatson::removeDuplicateEdges(std::vector<std::pair<Point, Point>>& edges) const {
    auto edgeComparator = [](const std::pair<Point, Point>& e1, const std::pair<Point, Point>& e2) {
        return (e1.first == e2.first && e1.second == e2.second) ||
               (e1.first == e2.second && e1.second == e2.first);
    };

    edges.erase(std::unique(edges.begin(), edges.end(), edgeComparator), edges.end());
}

void BowyerWatson::removeBadTriangles(const std::vector<Triangle>& badTriangles) {
    for (const auto& bad : badTriangles) {
        triangles.erase(std::remove(triangles.begin(), triangles.end(), bad), triangles.end());
    }
}

void BowyerWatson::triangulate(const std::vector<Point>& points, const Polygon& polygon) {
    Triangle superTriangle = getSuperTriangle(points);
    triangles.push_back(superTriangle);

    for (const auto& point : points) {
        std::vector<Triangle> badTriangles;

        for (const auto& triangle : triangles) {
            if (isPointInCircumcircle(point, triangle)) {
                badTriangles.push_back(triangle);
            }
        }

        removeBadTriangles(badTriangles);

        std::vector<std::pair<Point, Point>> edges = getPolygonalHoleEdges(badTriangles);
        removeDuplicateEdges(edges);

        for (const auto& edge : edges) {
            triangles.emplace_back(edge.first, edge.second, point);
        }
    }

    triangles.erase(
        std::remove_if(triangles.begin(), triangles.end(),
                       [&](const Triangle& t) {
                           return superTriangle.containsVertex(t.p1) ||
                                  superTriangle.containsVertex(t.p2) ||
                                  superTriangle.containsVertex(t.p3);
                       }),
        triangles.end());

    polygon.filterTriangles(triangles);
}


void BowyerWatson::printTriangles() const {
    std::cout << triangles.size() << std::endl;
    for (const auto& triangle : triangles) {
        std::cout << "Triangle: (" << triangle.p1.x << ", " << triangle.p1.y << ") -> ("
                    << triangle.p2.x << ", " << triangle.p2.y << ") -> ("
                    << triangle.p3.x << ", " << triangle.p3.y << ")\n";
    }
}

int main() {
    // Define points
    std::vector<Point> points = {
        Point(0, 0), Point(3, 0), Point(6, 2), Point(9, 2),
        Point(9, 4), Point(0, 4)
    };

    Polygon polygon;
    polygon.add_vertex(0,0);
    polygon.add_vertex(3,0);
    polygon.add_vertex(6,2);
    polygon.add_vertex(9,2);
    polygon.add_vertex(9,4);
    polygon.add_vertex(0,4);
    // Perform Delaunay triangulation
    BowyerWatson bowyerWatson;
    bowyerWatson.triangulate(points, polygon);

    // Print the resulting triangles
    bowyerWatson.printTriangles();
    // DelaunayVisualizer visualizer(points, bowyerWatson.triangles);
    // visualizer.run();
    return 0;
}
