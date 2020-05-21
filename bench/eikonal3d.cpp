#include "marlin/solver.hpp"

#if MARLIN_N_DIMS == 3
#    include <errno.h>
#    include <cmath>
#    include <cstdio>
#    include <cstdlib>
#    include <cstring>
#    include <vector>
#    include "bench.hpp"

using namespace marlin;

struct solver_startup_data
{
    std::vector<scalar_t> data;
    point_t size;
    std::array<std::pair<scalar_t, scalar_t>, dim> vertices;
    solver::params_t params;
};

static solver_startup_data get_solver_data(unsigned long domain_size);

static bench_stats run_once(solver_startup_data const& data);

struct arguments
{
    static constexpr unsigned c_max_domain_size = 1001;
    unsigned long domain_size = 101;
    unsigned long warmup_runs = 5;
    unsigned long min_compute_secs = 60;
};

static arguments handle_args(int argc, char** argv);

#endif

int main(int argc, char** argv)
{
#if MARLIN_N_DIMS == 3
    arguments args = handle_args(argc, argv);

    auto solver_data = get_solver_data(args.domain_size);

    for(int i = 0; (unsigned long)i < args.warmup_runs; ++i)
    {
        std::fprintf(
            stderr, "\rWarmup run %d/%lu        ", i + 1, args.warmup_runs);
        run_once(solver_data);
    }

    std::vector<bench_stats> stats;
    stats.reserve(1);

    std::fprintf(stderr,
                 "\rInitial measurement...                                ");
    stats.emplace_back(run_once(solver_data));

    int remaining = compute_n_runs(stats[0], args.min_compute_secs) - 1;
    stats.reserve(1 + remaining);

    for(; remaining > 0; remaining--)
    {
        report_current_avg(stats, remaining);
        stats.emplace_back(run_once(solver_data));
    }

    report_stats(stats);

#endif
    (void)argc;
    (void)argv;

    return 0;
}

#if MARLIN_N_DIMS == 3

solver_startup_data get_solver_data(unsigned long domain_size)
{
    solver_startup_data rv;

    rv.size = { domain_size, domain_size, domain_size };
    rv.vertices = { { { scalar_t{ -1.0 }, scalar_t{ 1.0 } },
                      { scalar_t{ -1.0 }, scalar_t{ 1.0 } },
                      { scalar_t{ -1.0 }, scalar_t{ 1.0 } } } };
    rv.params.maxval = 2.0;
    rv.params.tolerance = 1.0e-6;

    rv.data.resize(rv.size[0] * rv.size[1] * rv.size[2]);
    std::fill(std::begin(rv.data), std::end(rv.data), 1.0);

    marlin::point_t middle = { rv.size[0] / 2, rv.size[1] / 2, rv.size[2] / 2 };

    rv.data[middle[2] + rv.size[2] * (middle[1] + rv.size[1] * middle[0])] =
        -1.0;

    return rv;
}

bench_stats run_once(solver_startup_data const& sdata)
{
    bench_stats stats;

    auto hamiltonian = [](point_t const&, vector_t const& p) {
        return std::sqrt(std::inner_product(
            std::cbegin(p), std::cend(p), std::cbegin(p), 0.0));
    };

    auto viscosity = [](point_t const&) {
        return vector_t{ scalar_t{ 1.0 }, scalar_t{ 1.0 }, scalar_t{ 1.0 } };
    };

    chronometer c;

    solver::solver_t s(std::vector<scalar_t>(sdata.data),
                       sdata.size,
                       sdata.vertices,
                       sdata.params);

    stats.startup = c.get();

    c.reset();
    s.solve(hamiltonian, viscosity);
    stats.solve = c.get();

    return stats;
}

arguments handle_args(int argc, char** argv)
{
    arguments args;

    if(argc > 1 && ((std::strcmp(argv[1], "-h") == 0) ||
                    (std::strcmp(argv[1], "--help") == 0)))
    {
        std::fprintf(stderr,
                     "usage: %s DOMAIN_SIZE WARMUP_RUNS MIN_COMPUTE_SECS\n",
                     argv[0]);
        std::exit(0);
    }

    if(argc > 1)
    {
        errno = 0;
        args.domain_size = std::strtoul(argv[1], NULL, 10);

        if(errno != 0 || args.domain_size == 0)
            quit("Invalid argument for domain size.");

        if(args.domain_size > arguments::c_max_domain_size)
        {
            std::fprintf(
                stderr,
                "Warning: requested domain size is too large, clipping to %u",
                arguments::c_max_domain_size);
        }
    }

    if(argc > 2)
    {
        errno = 0;
        args.warmup_runs = std::strtoul(argv[2], NULL, 10);

        if(errno != 0)
            quit("Invalid argument for number of warmup runs");
    }

    if(argc > 3)
    {
        errno = 0;
        args.min_compute_secs = std::strtoul(argv[3], NULL, 10);

        if(errno != 0)
            quit("Invalid argument for minimum compute seconds");
    }

    return args;
}
#endif
