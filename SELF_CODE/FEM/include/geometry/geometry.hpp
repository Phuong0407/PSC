#ifndef GEOMETRY_GEOMETRY_HPP
#define GEOMETRY_GEOMETRY_HPP

#include <config.h>

namespace geometry
{
template<unsigned int spdim, typename real_t> class vector;
template<unsigned int spdim> class Point;
template<unsigned int spdim, typename real_t> class edge;
template<unsigned int spdim, typename real_t> class line;
template<unsigned int spdim, typename real_t> class polygon;
template<unsigned int spdim, typename real_t> class triangle;
template<typename real_t> class circle;

class Circle2D;
using Point2D = Point<2>;
}

#endif