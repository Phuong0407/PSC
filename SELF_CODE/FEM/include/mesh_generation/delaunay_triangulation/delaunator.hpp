#ifndef MESH_GENERATION_DELAUNAY_TRIANGULATION_DELAUNATOR_HPP
#define MESH_GENERATION_DELAUNAY_TRIANGULATION_DELAUNATOR_HPP

#include <mesh_generation/delaunay_triangulation/delaunator_helper.hpp>

#include <geometry/geometry.hpp>
#include <geometry/point.hpp>
#include <geometry/circle2D.hpp>

#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>

namespace mesh_generation
{

template<unsigned int spdim = 2, typename index_t = std::size_t>
class Delaunator {
private:
    std::vector<geometry::Point2D> coords;
    std::vector<index_t> triangles;
    std::vector<index_t> half_edges;
    std::vector<index_t> hull_prev;
    std::vector<index_t> hull_next;
    std::vector<index_t> hull_tri;
    index_t hull_start;

    geometry::Point2D center;
    index_t hash_size;
    std::vector<index_t> hash;
    std::vector<index_t> edge_stack;

    index_t legalizeEdge(index_t a) {
        index_t i = 0;
        index_t ar = 0;
        edge_stack.clear();

        while (true) {
            const index_t b = half_edges[a];
            const index_t a0 = 3 * (a / 3);
            ar = a0 + (a + 2) % 3;

            if (b == math_utils::integer<std::size_t>::INVALID_INDEX) {
                if (i > 0) {
                    i--;
                    a = edge_stack[i];
                    continue;
                } else {
                    break;
                }
            }

            const index_t b0 = 3 * (b / 3);
            const index_t al = a0 + (a + 1) % 3;
            const index_t bl = b0 + (b + 2) % 3;
    
            const index_t p0 = triangles[ar];
            const index_t pr = triangles[a];
            const index_t pl = triangles[al];
            const index_t p1 = triangles[bl];

            const bool illigal = geometry::Circle2D::isInSquaredRadiusCircle(coords[p0], coords[pr], coords[pl], coords[p1]);

            if (illigal) {
                triangles[a] = p1;
                triangles[b] = p0;
    
                auto hbl = half_edges[bl];
                
                if (hbl == math_utils::integer<std::size_t>::INVALID_INDEX) {
                    std::size_t e = hull_start;
                    do {
                        if (hull_tri[e] == bl) {
                            hull_tri[e] = a;
                            break;
                        }
                        e = hull_next[e];
                    } while (e != hull_start);
                }
                linkEdges(a, hbl);
                linkEdges(b, half_edges[ar]);
                linkEdges(ar, bl);
                std::size_t br = b0 + (b + 1) % 3;
    
                if (i < edge_stack.size()) {
                    edge_stack[i] = br;
                } else {
                    edge_stack.push_back(br);
                }
                i++;
            } else {
                if (i > 0) {
                    i--;
                    a = edge_stack[i];
                    continue;
                } else {
                    break;
                }
            }
        }
        return ar;
    }

    index_t getHashKey(const geometry::Point2D& P) const {
        const double dx = P.coordinate(0) - center.coordinate(0);
        const double dy = P.coordinate(1) - center.coordinate(1);
        return math_utils::integer<std::size_t>::mod(
            static_cast<index_t>(
                std::llround(std::floor(Delaunator_helper::pseudoAngle(dx, dy) * static_cast<double>(hash_size)))),
                hash_size
            );
    }

    index_t addTriangle(index_t i0, index_t i1, index_t i2, index_t a, index_t b, index_t c) {
        index_t t = triangles.size();
        triangles.push_back(i0);
        triangles.push_back(i1);
        triangles.push_back(i2);
        linkEdges(t, a);
        linkEdges(t + 1, b);
        linkEdges(t + 2, c);
        return t;
    }

    void linkEdges(index_t a, index_t b) {
        index_t size_half_edges = half_edges.size();
        if (a == size_half_edges) {
            half_edges.push_back(b);
        } else if (a < size_half_edges) {
            half_edges[a] = b;
        } else {
            throw std::runtime_error("Cannot link edge");
        }
    }

    const geometry::Point2D computeBoundingBoxCentroid() const {
        index_t N = coords.size();
        double max_x = std::numeric_limits<double>::min();
        double max_y = std::numeric_limits<double>::min();
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        for (index_t i = 0; i < N; i++) {
            const double x = coords[i].coordinate(0);
            const double y = coords[i].coordinate(1);
            
            if (x < min_x) min_x = x;
            if (y < min_y) min_y = y;
            if (x > max_x) max_x = x;
            if (y > max_y) max_y = y;
        }
        return geometry::Point2D((min_x + max_x)/2.0, (min_y + max_y)/2.0);
    }

    index_t computeSeedPointID(const geometry::Point2D& bounding_centroid) const {
        index_t N = coords.size();
        index_t seed_point_ID = 0;
        double min_dist = std::numeric_limits<double>::max();
        for (index_t i = 0; i < N; i++) {
            const geometry::Point2D point_iter = coords[i];
            const double distance = bounding_centroid.squaredDistTo(point_iter);
            if (distance < min_dist) {
                seed_point_ID = i;
                min_dist = distance;
            }
        }
        return seed_point_ID;
    }

    index_t computeFirstPointToSeedPointID(const index_t seed_point_ID) const {
        index_t N = coords.size();
        index_t first_point_to_seed_point_ID = 0;
        double min_dist = std::numeric_limits<double>::max();
        const geometry::Point2D seed_point = coords[seed_point_ID];
        for (index_t i = 0; i < N; ++i) {
            if (i == seed_point_ID)
                continue;
            const geometry::Point2D point_iter = coords[i];
            const double distance = seed_point.squaredDistTo(point_iter);
            if (distance < min_dist) {
                min_dist = distance;
                first_point_to_seed_point_ID = i;
            }
        }
        return first_point_to_seed_point_ID;
    }

    index_t computeSecondPointToSeedPointID(const index_t seed_point_ID, const index_t first_point_ID) const {
        index_t N = coords.size();
        index_t second_point_to_seed_point_ID = 0;
        double min_radius = std::numeric_limits<double>::max();
        const geometry::Point2D seed_point = coords[seed_point_ID];
        const geometry::Point2D first_point = coords[first_point_ID];
        for (index_t i = 0; i < N; ++i) {
            if (i == seed_point_ID || i == first_point_ID)
                continue;
            const geometry::Point2D point_iter = coords[i];
            const double radius = geometry::Circle2D::squardRadiusCircumcircle(seed_point, first_point, point_iter).getRadius();
            if (radius < min_radius) {
                min_radius = radius;
                second_point_to_seed_point_ID = i;
            }
        }
        if (!(min_radius < std::numeric_limits<double>::max()))
            throw std::runtime_error("Can not form a triangulation from seed point.");
        return second_point_to_seed_point_ID;
    }

    void initHashTable() {
        index_t N = coords.size();
        hash_size = static_cast<index_t>(std::llround(std::ceil(std::sqrt(N))));
        hash.resize(hash_size);
        std::fill(hash.begin(), hash.end(), math_utils::integer<std::size_t>::INVALID_INDEX);
    }

    void initHull(index_t i0, index_t i1, index_t i2) {
        index_t N = coords.size();
        hull_prev.resize(N);
        hull_next.resize(N);
        hull_tri.resize(N);

        hull_start = i0;

        const geometry::Point2D seed_point = coords[i0];
        const geometry::Point2D first_point = coords[i1];
        const geometry::Point2D second_point = coords[i2];


        hull_next[i0] = hull_prev[i2] = i1;
        hull_next[i1] = hull_prev[i0] = i2;
        hull_next[i2] = hull_prev[i1] = i0;
    
        hull_tri[i0] = 0;
        hull_tri[i1] = 1;
        hull_tri[i2] = 2;

        hash[getHashKey(seed_point)] = i0;
        hash[getHashKey(first_point)] = i1;
        hash[getHashKey(second_point)] = i2;
    }

public:

    Delaunator(const std::vector<geometry::Point2D>& coords) :
    coords(coords),
    triangles(), half_edges(),
    hull_prev(), hull_next(), hull_tri(), hull_start(),
    center(),
    hash_size(), hash(), edge_stack()
    {
        index_t N = coords.size();
        std::vector<index_t> ids(N);
        std::iota(ids.begin(), ids.end(), 0);

        const geometry::Point2D bounding_centroid = computeBoundingBoxCentroid();
        index_t seed_point_ID = computeSeedPointID(bounding_centroid);
        index_t first_point_to_seed_point_ID = computeFirstPointToSeedPointID(seed_point_ID);
        index_t second_point_to_seed_point_ID = computeSecondPointToSeedPointID(seed_point_ID, first_point_to_seed_point_ID);
        geometry::Point2D seed_point = coords[seed_point_ID];
        geometry::Point2D first_point = coords[first_point_to_seed_point_ID];
        geometry::Point2D second_point = coords[second_point_to_seed_point_ID];


        if(Delaunator_helper::orient(seed_point, first_point, second_point)) {
            std::swap(first_point_to_seed_point_ID, second_point_to_seed_point_ID);
            double x1 = first_point.coordinate(0);
            double y1 = first_point.coordinate(1);
            double x2 = second_point.coordinate(0);
            double y2 = second_point.coordinate(1);
            first_point.set({x2, y2});
            second_point.set({x1, y1});
        }
        center = geometry::Circle2D::squardRadiusCircumcircle(seed_point, first_point, second_point).getCenter();

        std::sort(ids.begin(), ids.end(), Delaunator_helper::compare{ coords, center });

        initHashTable();
        initHull(seed_point_ID, first_point_to_seed_point_ID, second_point_to_seed_point_ID);
        
        index_t hull_size = 3;

        index_t max_nutriangles = N < 3 ? 1 : 2 * N - 5;
        triangles.reserve(max_nutriangles * 3);
        half_edges.reserve(max_nutriangles * 3);

        addTriangle(
            seed_point_ID,
            first_point_to_seed_point_ID,
            second_point_to_seed_point_ID, math_utils::integer<std::size_t>::INVALID_INDEX, math_utils::integer<std::size_t>::INVALID_INDEX, math_utils::integer<std::size_t>::INVALID_INDEX
        );

        geometry::Point2D near_duplicate_point(
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
        );

        for (index_t k = 0; k < N; ++k) {
            const index_t i = ids[k];
            const geometry::Point2D point_iter = coords[i];
            if (k > 0 && near_duplicate_point == point_iter) {
                std::cout << "Fuck you "<< k << std::endl;
                continue;
            }
            
            near_duplicate_point.set(point_iter.coordinate());

            if (point_iter == seed_point ||
                point_iter == first_point ||
                point_iter == second_point)
                continue;
            
            index_t start = 0;
            index_t key = getHashKey(point_iter);
            for (index_t j = 0; j < hash_size; j++) {
                start = hash[math_utils::integer<std::size_t>::mod(key + j, hash_size)];
                if (start != math_utils::integer<std::size_t>::INVALID_INDEX && start != hull_next[start])
                    break;
            }

            start = hull_prev[start];
            index_t e = start;
            index_t q;

            while (q = hull_next[e], !Delaunator_helper::orient(point_iter, coords[e], coords[q])) {
                e = q;
                if (e == start) {
                    e = math_utils::integer<std::size_t>::INVALID_INDEX;
                    break;
                }
            }

            if (e == math_utils::integer<std::size_t>::INVALID_INDEX)
                continue;
            
            index_t t = addTriangle(e, i, hull_next[e], 
                math_utils::integer<std::size_t>::INVALID_INDEX, math_utils::integer<std::size_t>::INVALID_INDEX,
                hull_tri[e]);
            
            hull_tri[i] = legalizeEdge(t + 2);
            hull_tri[e] = t;
            hull_size++;

            index_t next = hull_next[e];
            while (
                q = hull_next[next],
                Delaunator_helper::orient(point_iter, coords[next], coords[q])
            ) {
                t = addTriangle(next, i, q, hull_tri[i], math_utils::integer<std::size_t>::INVALID_INDEX, hull_tri[next]);
                hull_tri[i] = legalizeEdge(t + 2);
                hull_next[next] = next;
                hull_size--;
                next = q;
            }

            if (e == start) {
                while (
                    q = hull_prev[e],
                    Delaunator_helper::orient(point_iter, coords[q], coords[e])
                ) {
                    t = addTriangle(q, i, e, math_utils::integer<std::size_t>::INVALID_INDEX, hull_tri[e], hull_tri[q]);
                    legalizeEdge(t + 2);
                    hull_tri[q] = t;
                    hull_next[e] = e;
                    hull_size--;
                    e = q;
                }
            }

            hull_prev[i] = e;
            hull_start = e;
            hull_prev[next] = i;
            hull_next[e] = i;
            hull_next[i] = next;
    
            hash[getHashKey(point_iter)] = i;
            hash[getHashKey(coords[e])] = e;
        }
    }

    double getHullArea() {
        std::vector<double> hull_area;
        size_t e = hull_start;
        do {
            hull_area.push_back((coords[e].coordinate(0) - coords[hull_prev[e]].coordinate(0)) * (coords[e].coordinate(1) + coords[hull_prev[e]].coordinate(1)));
            e = hull_next[e];
        } while (e != hull_start);
        return math_utils::floating_point<std::size_t>::kahan_sum(hull_area);
    }

    const std::vector<index_t>& getTriangle() const {
        return triangles;
    }

};

} // namespace mesh_generation

#endif // MESH_GENERATION_DELAUNAY_TRIANGULATION_DELAUNATOR_HPP