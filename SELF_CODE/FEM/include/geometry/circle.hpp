#ifndef GEOMETRY_CIRCLE_HPP
#define GEOMETRY_CIRCLE_HPP

#include <geometry/geometry.hpp>
#include <geometry/point.hpp>
#include <geometry/vector.hpp>

#include <array>
#include <cmath>

namespace geometry
{
    
    template<typename real_t = double> class circle {
        private:
            point<2, real_t> O;
            real_t R;
        
        public:
            circle() = default;
            ~circle() = default;

            explicit circle(point<2, real_t> O, real_t R) : O(O), R(R) {}
            explicit circle(std::array<real_t, 2> O, real_t R) : O(O), R(R) {}
            explicit circile(real_t x, real_t y, real_t O) : O(x, y), R(R) {}

            int relative_position(point<2, real_t> P) {
                real_t dis_to_center = O.dist(P);
                if (std::abs(dist_to_center - R) < ETOL<real_t>)
                    return -1;
                else if (std::abs(dist_to_center - R) == ETOL<real_t>)
                    return 0;
                else
                    return 1;
            }

            // static circle<real_t> circumcircle(const point<2, real_t> &P1, const point<2, real_t> &P2, const point<2, real_t> &P3) {
            //     point<2, real_t> midP1P2 = point<2, real_t>::mid_point(P1, P2);
            //     point<2, real_t> midP2P3 = point<2, real_t>::mid_point(P2, P3);

            //     vector<2, real_t> dirP1P2 = P2 - P1;
            //     vector<2, real_t> dirP2P3 = P3 - P2;

            //     vector<2, real_t> perpP1P2(-dirP1P2(1), dirP1P2(0));
            //     vector<2, real_t> perpP2P3(-dirP2P3(1), dirP2P3(0));

            //     real_t d = perpP1P2(0) * perpP2P3(1) - perpP1P2(1) * perpP2P3(0);

            //     if (std::abs(d) < std::numeric_limits<real_t>::epsilon()) {
            //         throw std::runtime_error("collinear points have no unique circumcircle.");
            //     }

            //     real_t t = ((midP2P3(0) - midP1P2(0)) * perpP2P3(1) - 
            //     (midP2P3(1) - midP1P2(1)) * perpP2P3(0)) / d;

            //     point<2, real_t> center = midP1P2 + perpP1P2 * t;

            //     real_t radius = center.dist(P1);

            //     return circle<real_t>(center, radius);
            // }

            static inline real_t circumradius(const point<2, real_t> &P1, const point<2, real_t> &P2, const point<2, real_t> &P3)
            {
                const vector<2, real_t> P1P2 = P2 - P1;
                const vector<2, real_t> P1P3 = P3 - P1;
                const real_t L2P1P2 = P1P2.norm();
                const real_t L2P1P3 = P1P3.norm();

                const cross_prod = P1P2.cross(P1P3);
                if (std::abs(cross_prod) < ETOL<real_t>) {
                    return std::numeric_limits<real_t>::max();
                }

                const double x = (P1P3(1) * L2P1P2 - P1P2(1) * L2P1P3) * 0.5 / cross_prod;
                const double y = (P1P2(0) * L2P1P3 - P1P3(0) * L2P1P2) * 0.5 / cross_prod;

                if ((L2P1P2 != 0.0) && (L2P1P3 != 0.0) && (cross_prod != 0.0)) {
                    return std::sqrt(x * x + y * y);
                } else {
                    return std::numeric_limits<double>::max();
                }
            }

            const point<2, real_t>& get_center() const {
                return O;
            }
            const real_t& get_radius() const {
                return R;
            }
    };


} // namespace geometry


#endif