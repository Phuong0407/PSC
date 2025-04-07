#ifndef GEOMETRY_CIRCLE_2D_HPP
#define GEOMETRY_CIRCLE_2D_HPP

#include <geometry/geometry.hpp>
#include <geometry/point.hpp>

#include <array>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace geometry {

    class Circle2D {
    private:
        Point2D center{};
        double radius{};

    public:
        constexpr Circle2D() noexcept = default;

        constexpr Circle2D(double cx, double cy, double r) noexcept
            : center{cx, cy}, radius(std::max(0.0, r)) {}

        explicit constexpr Circle2D(const std::array<double, 2>& c, double r) noexcept
            : center(c[0], c[1]), radius(std::max(0.0, r)) {}

        [[nodiscard]] constexpr const Point2D& getCenter() const noexcept {
            return center;
        }

        [[nodiscard]] constexpr double getRadius() const noexcept {
            return radius;
        }

        void set_center(double cx, double cy) noexcept {
            center = {cx, cy};
        }

        void set_radius(double r) {
            if (r < 0) throw std::invalid_argument("Radius cannot be negative.");
            radius = r;
        }

        [[nodiscard]] constexpr double distance_to_point(double x, double y) const noexcept {
            return std::hypot(center.coordinate(0) - x, center.coordinate(1) - y);
        }

        [[nodiscard]] constexpr bool contains(double x, double y) const noexcept {
            return distance_to_point(x, y) <= radius;
        }

        [[nodiscard]] constexpr bool intersects(const Circle2D& other) const noexcept {
            double d = std::hypot(center.coordinate(0) - other.center.coordinate(0), center.coordinate(1) - other.center.coordinate(1));
            return d <= (radius + other.radius);
        }

        void print(std::ostream& os = std::cout) const {
            os << "Circle2D[Center: (" << std::fixed << std::setprecision(6)
               << center.coordinate(0) << ", " << center.coordinate(1) << "), Radius: " << radius << "]";
        }

        static Circle2D squardRadiusCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3) {
            double dx = x2 - x1;
            double dy = y2 - y1;
            double ex = x3 - x1;
            double ey = y3 - y1;

            double bl = dx * dx + dy * dy;
            double cl = ex * ex + ey * ey;
            double d = dx * ey - dy * ex;

            double Cx = x1 + (ey * bl - dy * cl) * 0.5 / d;
            double Cy = y1 + (dx * cl - ex * bl) * 0.5 / d;

            double x = (ey * bl - dy * cl) *0.5 / d;
            double y = (dx * cl - ex * bl) *0.5 / d;

            double R = 0.0;
            if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) && (d > 0.0 || d < 0.0)) {
                R = x * x + y * y;
            } else {
                R = std::numeric_limits<double>::max();
            }

            return Circle2D(Cx, Cy, R);
        }

        static Circle2D squardRadiusCircumcircle(const Point2D& P1, const Point2D& P2, const Point2D& P3) {
            double x1 = P1.coordinate(0), y1 = P1.coordinate(1);
            double x2 = P2.coordinate(0), y2 = P2.coordinate(1);
            double x3 = P3.coordinate(0), y3 = P3.coordinate(1);
            return squardRadiusCircumcircle(x1, y1, x2, y2, x3, y3);
        }

        bool isInSquaredRadiusCircle(const Point2D& P) const {
            return center.squaredDistTo(P) <= radius + std::numeric_limits<double>::epsilon();
        }

        static bool isInSquaredRadiusCircle(const Point2D& P1, const Point2D& P2, const Point2D& P3, const Point2D& P) {
            const Circle2D& circle = squardRadiusCircumcircle(P1, P2, P3);
            return circle.isInSquaredRadiusCircle(P);
        }
    };

} // namespace geometry

#endif // GEOMETRY_CIRCLE_2D_HPP