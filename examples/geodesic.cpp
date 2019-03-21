#include <cmath>
#include <iostream>
#include <numeric>
#include "solver.hpp"

int main()
{
#if FSM_N_DIMS == 2
    std::array<std::pair<fsm::scalar_t, fsm::scalar_t>, fsm::dim> vertices = {
        { { -0.5, 0.5 }, { -0.5, 0.5 } }
    };

    std::array<int, fsm::dim> npts = { { 201, 201 } };

    auto pt = [&vertices, &npts ](auto x) -> auto
    {
        std::array<fsm::scalar_t, fsm::dim> pt;
        for(auto i = 0ul; i < fsm::dim; ++i)
        {
            pt[i] =
                vertices[i].first +
                x[i] * (vertices[i].second - vertices[i].first) / (npts[i] - 1);
        }

        return pt;
    };

    auto grad_x = [](auto const& x) -> fsm::scalar_t {
        return 0.9 * 2 * M_PI * std::cos(2 * M_PI * x[0]) *
               std::sin(2 * M_PI * x[1]);
    };

    auto grad_y = [](auto const& x) -> fsm::scalar_t {
        return 0.9 * 2 * M_PI * std::sin(2 * M_PI * x[0]) *
               std::cos(2 * M_PI * x[1]);
    };

    auto speed = [&grad_x, &grad_y](auto const& x,
                                    auto omega) -> fsm::scalar_t {
        auto gx = grad_x(x);
        auto gy = grad_y(x);
        auto s = std::sin(omega);
        auto ss = std::sin(2 * omega);
        auto c = std::cos(omega);
        auto num = 1 + gx * gx * s * s + gy * gy * c * c - gx * gy * ss;
        auto den = 1 + gx * gx + gy * gy;
        return std::sqrt(num / den);
    };

    auto h = [&pt, &speed](auto x, auto const& p) -> fsm::scalar_t {
        auto z = pt(x);
        auto norm = std::sqrt(
            std::inner_product(std::begin(p), std::end(p), std::begin(p), 0.0));
        auto omega = std::atan2(p[1], p[0]);
        return norm * speed(z, omega);
    };

    constexpr fsm::vector_t diss_coeffs = { { 1.0, 1.0 } };

    fsm::solver::solver_t s(
        "../data/geodesic.h5", h, vertices, diss_coeffs, 10.0, 1.0e-4);

    s.solve();
#endif
    return 0;
}
