#include <cmath>
#include <iostream>
#include <numeric>

#include "marlin/solver.hpp"

#include "hdf5.hpp"
#include "timer.hpp"

int main()
{
#if MARLIN_N_DIMS == 3
    constexpr char const* filename = "../data/eikonal3d.h5";

    auto h = [](auto, auto&& p) -> marlin::scalar_t {
        // Norm of p
        return std::sqrt(
            std::inner_product(std::begin(p), std::end(p), std::begin(p), 0.0));
    };

    constexpr std::array<std::pair<marlin::scalar_t, marlin::scalar_t>,
                         marlin::dim>
        vertices = { { { -1.0, 1.0 }, { -1.0, 1.0 }, { -1.0, 1.0 } } };

    auto viscosity = [](auto) { return marlin::vector_t{ 1.0, 1.0, 1.0 }; };

    auto ds = h5io::read(filename, "cost_function");

    marlin::solver::params_t params;
    params.tolerance = 1.0e-4;
    params.maxval = 2.0;

    marlin::solver::solver_t s(std::move(ds.data), ds.size, vertices, params);

    timer::timer_t t;

    s.solve(h, viscosity);

    std::cout << "Took " << t.get_elapsed_sec<double>() << " seconds."
              << std::endl;

    h5io::write(filename, "value_function", s.steal(), ds.size);

#endif
    return 0;
}
