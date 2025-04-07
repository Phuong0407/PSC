#ifndef THOMAS_ALGORITHM_H
#define THOMAS_ALGORITHM_H

#include <cmath>
#include <vector>
#include <cstddef>
#include <stdexcept>

void thomas_algorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& f) {
    std::size_t n = d.size();
    if (n == 0) return;
    std::vector<double> c_star(n, 0.0);
    std::vector<double> d_star(n, 0.0);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (std::size_t i = 1; i < n; i++) {
        double m = b[i] - a[i] * c_star[i - 1];
        if (m == 0)
            throw std::runtime_error("DIVISION BY ZERO IN THOMAS ALGORITHM.");
        m = 1.0 / m;
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    f[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        f[i] = d_star[i] - c_star[i] * f[i + 1];
    }
}

#endif