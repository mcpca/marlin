#include <cmath>
#include <cstdio>
#include <iostream>
#include <numeric>
#include "solver.hpp"

int main(int argc, char** argv)
{
#if FSM_N_DIMS == 2

    std::array<std::pair<fsm::scalar_t, fsm::scalar_t>, fsm::dim> vertices = {
        { { -1.0, 1.0 }, { -1.0, 1.0 } }
    };

    std::array<int, fsm::dim> npts = { { 201, 201 } };

    auto pt = [&vertices, &npts](auto x) -> auto
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

    fsm::scalar_t v[] = { -0.5, 0.5 };

    if(argc >= 3)
    {
        sscanf(argv[1], "%f", v);
        sscanf(argv[2], "%f", v + 1);
    }

    std::cout << "v = {" << v[0] << ", " << v[1] << "}" << std::endl;

    auto h = [&v](auto x, auto&& p) -> double {
        return std::sqrt(std::inner_product(
                   std::begin(p), std::end(p), std::begin(p), 0.0)) -
               p[0] * v[0] - p[1] * v[1];
    };

    fsm::vector_t diss_coeffs;
    diss_coeffs[0] = 1.0 + std::fabs(v[0]);
    diss_coeffs[1] = 1.0 + std::fabs(v[1]);

    fsm::solver::solver_t s(
        "../data/constant_field.h5", h, vertices, diss_coeffs, 5.0, 1.0e-4);

    s.solve();
#endif
    return 0;
}
