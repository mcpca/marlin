#include <cmath>
#include <iostream>
#include <numeric>

#include "fsm/solver.hpp"

#include "timer.hpp"

int main()
{
#if FSM_N_DIMS == 2
    auto h = [](auto x, auto&& p) -> fsm::scalar_t {
        // Norm of p
        return std::sqrt(
            std::inner_product(std::begin(p), std::end(p), std::begin(p), 0.0));
    };

    constexpr std::array<std::pair<fsm::scalar_t, fsm::scalar_t>, fsm::dim>
        vertices = { { { -1.0, 1.0 }, { -1.0, 1.0 } } };

    auto viscosity = [](auto const& x) { return fsm::vector_t{ 1.0, 1.0 }; };

    fsm::solver::params_t params;
    params.tolerance = 1.0e-4;
    params.maxval = 2.0;

    fsm::solver::solver_t s(
        "../data/eikonal2d.h5", h, vertices, viscosity, params);

    timer::timer_t t;

    s.solve();

    std::cout << "Took " << t.get_elapsed_sec<double>() << " seconds."
              << std::endl;

#endif
    return 0;
}
