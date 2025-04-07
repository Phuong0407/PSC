#ifndef MATH_UTILS_FLOATING_POINT_HPP
#define MATH_UTILS_FLOATING_POINT_HPP

#include <vector>
#include <cmath>

namespace math_utils
{
    
template<typename index_t = std::size_t>
class floating_point {
public:

    static inline double kahan_sum(const std::vector<double>& x) {
        index_t N = x.size();
        double sum = x[0];
        double err = 0.0;

        for (index_t i = 1; i < N; ++i) {
            const double k = x[i];
            const double m = sum + k;
            err += std::fabs(sum) >= std::fabs(k) ? sum - m + k : k - m + sum;
            sum = m;
        }
        return sum + err;
    }
};


} // namespace math_utils


#endif