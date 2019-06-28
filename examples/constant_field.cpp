#include <cmath>
#include <cstdio>
#include <iostream>
#include <numeric>

#include "fsm/solver.hpp"

#include "timer.hpp"

int main(int argc, char** argv)
{
#if FSM_N_DIMS == 2

    std::array<std::pair<fsm::scalar_t, fsm::scalar_t>, fsm::dim> vertices = {
        { { -1.0, 1.0 }, { -1.0, 1.0 } }
    };

    fsm::scalar_t v[] = { -0.5, 0.5 };

    if(argc >= 3)
    {
        sscanf(argv[1], "%f", v);
        sscanf(argv[2], "%f", v + 1);
    }

    auto h = [&v](auto x, auto&& p) -> fsm::scalar_t {
        return std::sqrt(std::inner_product(
                   std::begin(p), std::end(p), std::begin(p), 0.0)) -
               p[0] * v[0] - p[1] * v[1];
    };

    auto viscosity = [&v](auto const& x) {
        return fsm::vector_t{ static_cast<fsm::scalar_t>(1.0 + std::abs(v[0])),
                              static_cast<fsm::scalar_t>(1.0 +
                                                         std::abs(v[1])) };
    };

    fsm::solver::params_t params;
    params.tolerance = 1.0e-4;
    params.maxval = 5.0;

    fsm::solver::solver_t s(
        "../data/constant_field.h5", h, vertices, viscosity, params);

    timer::timer_t t;

    s.solve();

    std::cout << "Took " << t.get_elapsed_sec<double>() << " seconds."
              << std::endl;

#endif
    return 0;
}
