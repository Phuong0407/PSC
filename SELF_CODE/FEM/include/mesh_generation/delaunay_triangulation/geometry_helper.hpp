#ifndef MESH_GENERATION_DELAUNAY_TRIANGULATION_GEOMETRY_HELPER_HPP
#define MESH_GENERATION_DELAUNAY_TRIANGULATION_GEOMETRY_HELPER_HPP

#include <vector>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <type_traits>
#include <iostream>

namespace mesh_generation {

constexpr unsigned int INVALID_COORDINATE_INDEX = std::numeric_limits<unsigned int>::max();
constexpr double COORDINATE_TOLERANCE = 1e-9;

template<unsigned int dim>
class VertexCoordinate {
    static_assert(dim == 2 || dim == 3, "The discretization supports only 2D or 3D domains.");

private:
    std::vector<double> coordinates;

    inline unsigned int getVertexCoordinatePosition(unsigned int vertexIndex, unsigned int CoordinatOffset) {
        if (coordinates.empty())
            return INVALID_COORDINATE_INDEX;

        unsigned int coordinate_index = dim * vertexIndex + CoordinatOffset;
        if (coordinate_index >= coordinates.size())
            return INVALID_COORDINATE_INDEX;
        return coordinate_index;
    }

public:
    VertexCoordinate() = default;
    ~VertexCoordinate() = default;

    explicit VertexCoordinate(std::vector<double> coordinates)
    : coordinates(std::move(coordinates)) {}
    
    explicit VertexCoordinate(unsigned int estimated_vertices)
    : coordinates() {
        coordinates.reserve(estimated_vertices);
    }

    VertexCoordinate(const VertexCoordinate&) = delete;
    VertexCoordinate(VertexCoordinate&&) noexcept = default;
    VertexCoordinate& operator=(VertexCoordinate&&) noexcept = default;

    double getVertexCoordinate(unsigned int vertexIndex, unsigned int CoordinatOffset) {
        unsigned int coordinate_index = getVertexCoordinatePosition(vertexIndex, CoordinatOffset);
        if (vertexIndex != INVALID_COORDINATE_INDEX)
            throw std::invalid_argument("Invalid vertex index when getting vertex coordinate.");
        
        return coordinates[coordinate_index];
    }

    template <typename... Args>
    void addVertex(Args...coords) {
        static_assert(sizeof...(coords) == dim, "The number of coordinates must match the dimension of problem.");

        const double new_coords[] = { static_cast<double>(coords)... };
        for (unsigned int i = 0; i < coordinates.size(); i += dim) {
            bool exists = true;
            for (size_t j = 0; j < dim; ++j) {
                if (std::abs(coordinates[i + j] - new_coords[j]) > COORDINATE_TOLERANCE) {
                    exists = false;
                    break;
                }
            }
            if (exists) {
                return;
            }
        }
        (coordinates.push_back(coords), ...);
    }

    void printCoordinates() {
        std::cout << "Coordinates:\n";
        for (size_t i = 0; i < coordinates.size(); i += dim) {
            std::cout << "(";
            for (size_t j = 0; j < dim; ++j) {
                std::cout << coordinates[i + j] << (j < dim - 1 ? ", " : "");
            }
            std::cout << ")\n";
        }
    }

};

} // namespace mesh_generation

#endif // MESH_GENERATION_DELAUNAY_TRIANGULATION_GEOMETRY_HELPER_HPP