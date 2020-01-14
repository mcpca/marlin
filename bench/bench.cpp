#include <algorithm>
#include <cmath>
#include <cstdio>

#include "bench.hpp"

int compute_n_runs(bench_stats const& first_run, int min_time)
{
    int n_runs =
        static_cast<int>(std::floor((double)min_time / (first_run.solve)));

    // Round to next multiple of 10
    return n_runs + 10 - (n_runs % 10);
}

struct statistics
{
    double average = 0;
    double median = 0;
    double min = std::numeric_limits<double>::max();
    double max = 0;
    double std = 0;
};

void report_stats(std::vector<bench_stats> const& stats)
{
    statistics startup, compute;

    std::vector<double> sorted;
    sorted.reserve(stats.size());

#define COMPUTE_STATS(x, field)                                              \
    for(auto&& stat : stats)                                                 \
    {                                                                        \
        x.average += stat.field;                                             \
        x.min = std::min(stat.field, x.min);                                 \
        x.max = std::max(stat.field, x.max);                                 \
        sorted.insert(std::upper_bound(                                      \
                          std::begin(sorted), std::end(sorted), stat.field), \
                      stat.field);                                           \
    }                                                                        \
                                                                             \
    x.median = sorted[sorted.size() / 2];                                    \
    sorted.clear();                                                          \
                                                                             \
    x.average /= stats.size();                                               \
                                                                             \
    for(auto&& stat : stats)                                                 \
        x.std += std::pow(x.average - stat.field, 2);                        \
                                                                             \
    x.std = std::sqrt(x.std / (stats.size() - 1));

    COMPUTE_STATS(startup, startup);
    COMPUTE_STATS(compute, solve);

#undef COMPUTE_STATS

    std::fprintf(stderr,
                 "\rStatistics:                                                "
                 "        \n");

    std::fprintf(
        stderr,
        "\tStartup time:\tavg=%.2e, std=%.2e, range: %.2e .. %.2e .. %.2e\n",
        startup.average,
        startup.std,
        startup.min,
        startup.median,
        startup.max);

    std::fprintf(
        stderr,
        "\tCompute time:\tavg=%.2e, std=%.2e, range: %.2e .. %.2e .. %.2e\n",
        compute.average,
        compute.std,
        compute.min,
        compute.median,
        compute.max);
}

void report_current_avg(std::vector<bench_stats> const& stats, int remaining)
{
    double startup_avg = 0;
    double compute_avg = 0;

    for(auto&& stat : stats)
    {
        startup_avg += stat.startup;
        compute_avg += stat.solve;
    }

    startup_avg /= stats.size();
    compute_avg /= stats.size();

    double total = startup_avg + compute_avg;

    std::fprintf(stderr,
                 "\rCurrent estimate (%d left): %.2f + %.2f = %.2f "
                 "sec, ETA: %d sec        ",
                 remaining,
                 startup_avg,
                 compute_avg,
                 total,
                 (int)total * remaining);
}

void quit(const char* msg)
{
    std::fprintf(stderr, "%s\n", msg);
    std::exit(1);
}
