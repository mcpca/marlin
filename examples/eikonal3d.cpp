#include <cmath>
#include <iostream>
#include <numeric>

#include "marlin/solver.hpp"

#include "timer.hpp"

int main()
{
#if MARLIN_N_DIMS == 3
    auto h = [](auto x, auto&& p) -> marlin::scalar_t {
        // Norm of p
        return std::sqrt(
            std::inner_product(std::begin(p), std::end(p), std::begin(p), 0.0));
    };

    constexpr std::array<std::pair<marlin::scalar_t, marlin::scalar_t>, marlin::dim>
        vertices = { { { -1.0, 1.0 }, { -1.0, 1.0 }, { -1.0, 1.0 } } };

    auto viscosity = [](auto x){ return marlin::vector_t{1.0, 1.0, 1.0}; };

    marlin::solver::params_t params;
    params.tolerance = 1.0e-4;
    params.maxval = 2.0;

    marlin::solver::solver_t s(
        "../data/eikonal3d.h5", h, vertices, viscosity, params);

    timer::timer_t t;

    s.solve();

    std::cout << "Took " << t.get_elapsed_sec<double>() << " seconds."
              << std::endl;
#endif
    return 0;
}
