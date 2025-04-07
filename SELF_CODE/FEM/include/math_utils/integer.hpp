#ifndef MATH_UTILS_INTEGER_HPP
#define MATH_UTILS_INTEGER_HPP

#include <cstdlib>
#include <limits>

namespace math_utils
{
    
template<typename index_t = std::size_t>
class integer {
public:

    static constexpr index_t INVALID_INDEX = std::numeric_limits<index_t>::max();

    static inline index_t mod(index_t i, index_t c) {
        return i >= c ? i % c : i;
    }
};
    
} // namespace math_utils


#endif