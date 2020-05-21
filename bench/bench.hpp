#pragma once

#include <chrono>
#include <vector>

class chronometer
{
    using clock_ = std::chrono::system_clock;
    clock_::time_point m_last;

  public:
    chronometer() : m_last(clock_::now()) {}

    double get() const
    {
        auto const now = clock_::now();

        std::chrono::duration<double> elapsed = now - m_last;

        return elapsed.count();
    }

    void reset() { m_last = clock_::now(); }
};

struct bench_stats
{
    double startup = 0;
    double solve = 0;
};

int compute_n_runs(bench_stats const& first_run, int min_time);

void report_current_avg(std::vector<bench_stats> const& stats, int remaining);

void report_stats(std::vector<bench_stats> const& stats);

void quit(const char* msg);
