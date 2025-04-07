#ifndef GEOMETRY_DELAUNAY_TRIANGULATION_CONVEX_HULL_HPP
#define GEOMETRY_DELAUNAY_TRIANGULATION_CONVEX_HULL_HPP

#include <geometry/geometry.hpp>

#include <vector>
#include <algorithm>

namespace geometry
{

    template<unsigned int spdim = 2, typename real_t = double, typename index_t = unsigned int>
    class convex_hull {
        private:
            static real_t cross_product(
                const point<spdim, real_t>& O,
                const point<spdim, real_t>& A,
                const point<spdim, real_t>& B)
            {
                return (A(0) - O(0)) * (B(1) - O(1)) - (A(1) - O(1)) * (B(0) - O(0));
            }

        public:
            static std::vector<point<spdim, real_t>> compute_convex_hull(std::vector<point<spdim, real_t>> points) {
                if (points.size() <= 1)
                    throw std::invalid_argument("Number of points is smaller or equal to 1.");
    
                std::sort(points.begin(), points.end(), [](const point<spdim, real_t>& a, const point<spdim, real_t>& b) {
                    return (a(0) < b(0)) || (a(0) == b(0) && a(1) < b(1));
                });
    
                std::vector<point<spdim, real_t>> lower;
                for (const auto& p : points) {
                    while (lower.size() >= 2 && cross_product(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
                        lower.pop_back();
                    }
                    lower.push_back(p);
                }
    
                std::vector<point<spdim, real_t>> upper;
                for (auto it = points.rbegin(); it != points.rend(); ++it) {
                    while (upper.size() >= 2 && cross_product(upper[upper.size() - 2], upper[upper.size() - 1], *it) <= 0) {
                        upper.pop_back();
                    }
                    upper.push_back(*it);
                }
    
                lower.pop_back();
                upper.pop_back();
    
                lower.insert(lower.end(), upper.begin(), upper.end());
    
                return lower;
            }

    };

} // namespace geometry


#endif