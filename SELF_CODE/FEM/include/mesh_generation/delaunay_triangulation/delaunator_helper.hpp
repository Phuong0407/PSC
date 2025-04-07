#ifndef MESH_GENERATION_DELAUNATOR_HELPER_HPP
#define MESH_GENERATION_DELAUNATOR_HELPER_HPP

#include <math_utils/integer.hpp>
#include <math_utils/floating_point.hpp>

#include <geometry/geometry.hpp>
#include <geometry/point.hpp>

#include <cmath>

namespace mesh_generation
{

class Delaunator_helper {
public:
    static inline bool orient(const geometry::Point2D& P1, const geometry::Point2D& P2, const geometry::Point2D& P3) {
        const double px = P1.coordinate(0), py = P1.coordinate(1);
        const double qx = P2.coordinate(0), qy = P2.coordinate(1);
        const double rx = P3.coordinate(0), ry = P3.coordinate(1);
        return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) > 0.0;
    }

    static inline double pseudoAngle(const double dx, const double dy) {
        const double p = dx / (std::abs(dx) + std::abs(dy));
        return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0;
    }

    struct compare {
        const std::vector<geometry::Point2D>& coords;
        const geometry::Point2D& center;
    
        bool operator()(std::size_t i, std::size_t j) {
            const geometry::Point2D& P1 = coords[i];
            const geometry::Point2D& P2 = coords[j];
            const double d1 = center.squaredDistTo(P1);
            const double d2 = center.squaredDistTo(P2);
            const double diff1 = d1 - d2;
            const double diff2 = P1.coordinate(0) - P2.coordinate(0);
            const double diff3 = P1.coordinate(1) - P2.coordinate(1);
    
            if (diff1 > 0.0 || diff1 < 0.0) {
                return diff1 < 0;
            } else if (diff2 > 0.0 || diff2 < 0.0) {
                return diff2 < 0;
            } else {
                return diff3 < 0;
            }
        }
    };
};

} // namespace mesh_generation

#endif // MESH_GENERATION_DELAUNATOR_HELPER_HPP