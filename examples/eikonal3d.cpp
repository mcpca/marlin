#include <cmath>
#include <iostream>
#include <numeric>
#include "solver.hpp"

int main()
{
#if FSM_N_DIMS == 3
    auto h = [](auto x, auto&& p) -> double {
        // Norm of p
        return std::sqrt(
            std::inner_product(std::begin(p), std::end(p), std::begin(p), 0.0));
    };

    constexpr std::array<std::pair<fsm::scalar_t, fsm::scalar_t>, fsm::dim>
        vertices = { { { -1.0, 1.0 }, { -1.0, 1.0 }, { -1.0, 1.0 } } };

    constexpr fsm::vector_t diss_coeffs = { { 1.0, 1.0, 1.0 } };

    fsm::solver::solver_t s(
        "../data/eikonal3d.h5", h, vertices, diss_coeffs, 2.0, 1.0e-4);

    s.solve();
#endif
    return 0;
}
