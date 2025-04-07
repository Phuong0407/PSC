#ifndef GEOMETRY_POINT_HPP
#define GEOMETRY_POINT_HPP

#include <geometry/geometry.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>

namespace geometry {

template<unsigned int spdim>
class alignas(32) Point {
    static_assert(spdim == 2 || spdim == 3, "Point must be in 2D or 3D space.");

private:
    std::array<double, spdim> coords{};

public:
    Point() = default;
    ~Point() = default;

    constexpr Point(double x, double y) noexcept requires (spdim == 2)
        : coords{x, y} {}

    constexpr Point(double x, double y, double z) noexcept requires (spdim == 3)
        : coords{x, y, z} {}

    explicit constexpr Point(std::array<double, spdim>&& c) noexcept
        : coords(std::move(c)) {}

    explicit constexpr Point(const std::array<double, spdim>& c) noexcept
        : coords(c) {}

    Point(const Point&) = default;
    Point& operator=(const Point&) = default;

    Point(Point&&) noexcept = default;
    Point& operator=(Point&&) noexcept = default;

    [[nodiscard]] constexpr bool operator==(const Point<spdim>& other) const noexcept {
        for (std::size_t i = 0; i < spdim; ++i) {
            if (std::fabs(this->coords[i] - other.coordinate(i)) > std::numeric_limits<double>::epsilon()) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] constexpr bool operator!=(const Point<spdim>& other) const noexcept {
        return !(*this == other);
    }

    [[nodiscard]] constexpr const std::array<double, spdim>& coordinate() const noexcept {
        return coords;
    }

    [[nodiscard]] constexpr double coordinate(std::size_t i) const {
        if (i >= spdim)
            throw std::out_of_range("Index out of bounds.");
        return coords[i];
    }

    void set(std::size_t i, double value) {
        if (i >= spdim)
            throw std::out_of_range("Index out of bounds.");
        coords[i] = value;
    }

    void set(const std::array<double, spdim>& new_coords) noexcept {
        coords = std::move(new_coords);
    }

    [[nodiscard]] constexpr double distTo(const Point& other) const noexcept {
        double sum = 0.0;
        for (std::size_t i = 0; i < spdim; ++i) {
            sum += (coords[i] - other.coords[i]) * (coords[i] - other.coords[i]);
        }
        return std::sqrt(sum);
    }

    [[nodiscard]] constexpr double squaredDistTo(const Point& other) const noexcept {
        double sum = 0.0;
        for (std::size_t i = 0; i < spdim; ++i) {
            sum += (coords[i] - other.coords[i]) * (coords[i] - other.coords[i]);
        }
        return sum;
    }

    [[nodiscard]] constexpr double norm() const noexcept {
        double sum = 0.0;
        for (std::size_t i = 0; i < spdim; ++i) {
            sum += coords[i] * coords[i];
        }
        return sum;
    }

    void print(std::ostream& os = std::cout) const {
        os << "(";
        for (std::size_t i = 0; i < spdim; ++i) {
            os << std::fixed << std::setprecision(6) << coords[i];
            if (i < spdim - 1) os << ", ";
        }
        os << ")";
    }
};

} // namespace geometry

#endif // GEOMETRY_POINT_HPP