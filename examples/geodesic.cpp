#include <cmath>
#include <iostream>
#include <numeric>

#include "marlin/solver.hpp"

#include "timer.hpp"

int main()
{
#if MARLIN_N_DIMS == 2
    constexpr std::array<std::pair<marlin::scalar_t, marlin::scalar_t>, marlin::dim>
        vertices = { { { -0.5, 0.5 }, { -0.5, 0.5 } } };

    constexpr std::array<int, marlin::dim> npts = { { 201, 201 } };

    auto pt = [&vertices, &npts](auto x) {
        std::array<marlin::scalar_t, marlin::dim> pt;
        for(auto i = 0ul; i < marlin::dim; ++i)
        {
            pt[i] =
                vertices[i].first +
                x[i] * (vertices[i].second - vertices[i].first) / (npts[i] - 1);
        }

        return pt;
    };

    auto grad = [](auto const& x) {
        return marlin::vector_t{ static_cast<marlin::scalar_t>(
                                  0.9 * 2 * M_PI * std::cos(2 * M_PI * x[0]) *
                                  std::sin(2 * M_PI * x[1])),
                              static_cast<marlin::scalar_t>(
                                  0.9 * 2 * M_PI * std::sin(2 * M_PI * x[0]) *
                                  std::cos(2 * M_PI * x[1])) };
    };

    auto speed = [&grad](auto const& x, auto omega) -> marlin::scalar_t {
        auto g = grad(x);
        auto s = std::sin(omega);
        auto ss = std::sin(2 * omega);
        auto c = std::cos(omega);
        auto num =
            1 + g[0] * g[0] * s * s + g[1] * g[1] * c * c - g[0] * g[1] * ss;
        auto den = 1 + g[0] * g[0] + g[1] * g[1];
        return std::sqrt(num / den);
    };

    auto h = [&pt, &speed](auto x, auto const& p) -> marlin::scalar_t {
        auto z = pt(x);
        auto norm = std::sqrt(
            std::inner_product(std::begin(p), std::end(p), std::begin(p), 0.0));
        auto omega = std::atan2(p[1], p[0]);
        return norm * speed(z, omega);
    };

    auto viscosity = [](auto x) { return marlin::vector_t{ 1.0, 1.0 }; };

    marlin::solver::params_t params;
    params.tolerance = 1.0e-4;
    params.maxval = 2.0;

    marlin::solver::solver_t s(
        "../data/geodesic.h5", h, vertices, viscosity, params);

    timer::timer_t t;

    s.solve();

    std::cout << "Took " << t.get_elapsed_sec<double>() << " seconds."
              << std::endl;

#endif
    return 0;
}
