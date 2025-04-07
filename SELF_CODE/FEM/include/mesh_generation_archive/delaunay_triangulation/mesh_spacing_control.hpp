#ifndef MESH_GENERATION_DELAUNAY_TRIANGULATION_MESH_SPACING_CONTROL_HPP
#define MESH_GENERATION_DELAUNAY_TRIANGULATION_MESH_SPACING_CONTROL_HPP

#include <geometry/point.hpp>

#include <cmath>
#include <limits>
#include <algorithm>

namespace mesh_generation
{

    template<unsigned int spdim = 2, typename real_t = double, typename index_t = unsigned int>
    class mesh_control_point_t {
        private:
            std::vector<geometry::point<spdim, real_t>> positions;
            std::vector<real_t> amplifications;
            std::vector<real_t> decays;

        public:
            mesh_control_point_t() = default;
            ~mesh_control_point_t() = default;

            mesh_control_point_t(
                std::vector<geometry::point<spdim, real_t>> positions,
                std::vector<real_t> amplifications,
                std::vector<real_t> decays
            ) : positions(positions),
                amplifications(amplifications),
                decays(decays)
            {
                if (positions.size() != amplifications.size())
                    throw std::invalid_argument("Number of control point does not match the number of amplifications data");
                if (positions.size() != decays.size())
                    throw std::invalid_argument("Number of control point does not match the number of decays data");
            }

            real_t compute_control_point_spacing(const geometry::point<spdim, real_t>& P) const {
                real_t dpi = std::numeric_limits<real_t>::infinity();
                index_t num_control_points = positions.size();
                for (index_t i = 0; i < num_control_points; ++i) {
                    real_t distance = positions[i].dist(P);
                    real_t ndpi = amplifications[i] * std::exp(decays[i] * distance);
                    dpi = std::min(dpi, ndpi);
                }
                return dpi;
            }
    };

} // namespace mesh_generation


#endif