#ifndef MESH_GENERATION_DELAUNAY_DELAUNAY_TRIANGULATION_HPP
#define MESH_GENERATION_DELAUNAY_DELAUNAY_TRIANGULATION_HPP

#include <geometry/polygon.hpp>
#include <mesh_generation/delaunay_triangulation/mesh_spacing_control.hpp>

namespace mesh_generation
{
    template<unsigned int spdim = 2, typename real_t = double, typename index_t = unsigned int>
    class delaunay_triangulation{
        private:
            std::vector<geometry::point<spdim, real_t>> boundary;
            std::vector<real_t> dpi;
            mesh_control_point_t<spdim, real_t, index_t> mesh_control_point;

            void compute_point_distribution_on_boundary() {
                unsigned int N = boundary.size();
                for (unsigned int i = 0; i < N; ++i) {
                    const geometry::point<spdim, real_t>& Pi = boundary[i];
                    real_t distancePP1 = Pi.dist(boundary[(i+1) % N]);
                    real_t distancePM1 = Pi.dist(boundary[(i + N - 1) % N]);
                    dpi[i] = 0.5 * (distancePP1 + distancePM1);
                }
            }

            std::vector<Triangle> delaunay_triangulation_boundary() {
                return triangulate(boundary);
            }
            

        public:
            delaunay_triangulation() = default;
            ~delaunay_triangulation() = default;

            void init_mesh_control_point(
                std::vector<geometry::point<spdim, real_t>> positions,
                std::vector<real_t> amplifications,
                std::vector<real_t> decays)
            {
                mesh_control_point(positions, amplifications, decays);
            }

            void compute_delaunay_triangulation(const std::vector<geometry::point<spdim, real_t>> &P) {
                // step 1: compute point distribution function dpi for each boundary
                compute_point_distribution_on_boundary();
                // step 2: perform triangulation for the boundary

                // step 3: initiate j = 0

                // step 4: for each triangle Tm within the domain, perform the following:

                // 4(a) Define a point at centroid Tm with node n1, n2, n3
                // with Pc = 1/3 (rn1 + rn2 + rn3)

                // 4(b) Derive the point distribution dpc by interpolation the point distribution function from the node n1, n2, n3
                // dpc = 1/3 (dpn1 + dpn2 + dpn3)

                // 4(c) 


                // step 5: if j = 0, go to 7

                // step 6: perform delaunay triangulation of the derived points Pj, j = 1 --> N. Go to 3

                // step 7: smooth the grid, using Laplace filter


            }

    };

} // namespace mesh_generation


#endif