////////////////////////////////////////////////////////////////////////////////
// MIT License
//
// Copyright (c) 2019 Miguel Aguiar
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
// https://github.com/mcpca/fsm

#include <algorithm>
#include <iostream>
#include <numeric>

#include "data.hpp"
#include "fsm/solver.hpp"
#include "grid.hpp"
#include "io.hpp"
#include "queue.hpp"
#include "solver_internals.hpp"

namespace fsm
{
    namespace solver
    {
        solver_t make_solver(
            std::string const& filename,
            hamiltonian_t const& hamiltonian,
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            std::function<vector_t(input_t const&)> const& viscosity,
            params_t const& params)
        {
            auto data = io::read(filename, std::string("cost_function"));

            if(io::dset_exists(filename, "value_function"))
            {
                throw std::runtime_error(
                    "A dataset named \'value_function\' already exists in the "
                    "given file.");
            }

            grid::grid_t g(vertices, data.size);
            return solver_t(filename,
                            hamiltonian,
                            std::move(data.data),
                            g,
                            viscosity,
                            params);
        }

        solver_t::solver_t(solver_t&&) noexcept = default;
        solver_t& solver_t::operator=(solver_t&&) noexcept = default;
        solver_t::~solver_t() = default;

        solver_t::solver_t(
            std::string const& filename,
            hamiltonian_t const& hamiltonian,
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            std::function<vector_t(input_t const&)> const& viscosity,
            params_t const& params)
            : solver_t(make_solver(filename,
                                   hamiltonian,
                                   vertices,
                                   viscosity,
                                   params))
        {}

        solver_t::solver_t(
            std::string const& filename,
            hamiltonian_t const& hamiltonian,
            data::data_t cost,
            grid::grid_t const& grid,
            std::function<vector_t(input_t const&)> const& viscosity,
            params_t const& params)
            : m_filename(filename),
              m_hamiltonian(hamiltonian),
              m_grid(std::make_unique<grid::grid_t>(grid)),
              m_soln(std::make_unique<data::data_t>(m_grid->npts(),
                                                    params.maxval)),
              m_cost(std::make_unique<data::data_t>(std::move(cost))),
              m_viscosity(viscosity),
              m_tolerance(params.tolerance),
              m_pool(std::make_unique<ThreadPool>(parallel::n_workers - 1)),
              m_queue(std::make_unique<queue::queue_t>())
        {
            m_worker.reserve(parallel::n_workers);

            for(unsigned i = 0; i < parallel::n_workers; ++i)
            {
                m_worker.emplace_back(std::make_unique<data::data_t>(
                    m_grid->npts(), params.maxval));

                if(i < parallel::n_workers - 1)
                {
                    m_queue->enqueue(m_worker.back().get());
                }
            }
        }

        void solver_t::solve()
        {
            std::cout << "Initializing problem instance..." << std::endl;
            initialize();
            std::cout << "Initialized. Solving (" << parallel::n_workers
                      << " threads)..." << std::endl;

#ifndef NDEBUG
            auto niter = 0;
            std::cerr << "Iteration " << niter++ << ":\n";
#endif

            while(!iterate())
            {
#ifndef NDEBUG
                std::cerr << "Iteration " << niter++ << ":\n";
#else
                ;
#endif
            }

            std::cout << "Done. Writing to " << m_filename << "..."
                      << std::endl;
            io::write(m_filename, "value_function", *m_soln, m_grid->size());
        }

        void solver_t::initialize()
        {
#ifndef NDEBUG
            auto ntargetpts = 0;
#endif
            for(index_t i = 0; i < m_cost->size(); ++i)
            {
                if(m_cost->at(i) < scalar_t{ 0.0 })
                {
#ifndef NDEBUG
                    ++ntargetpts;
#endif
                    m_soln->at(i) = -(m_cost->at(i) + scalar_t{ 1.0 });

                    for(unsigned j = 0; j < m_worker.size(); ++j)
                    {
                        m_worker[j]->at(i) = m_soln->at(i);
                    }
                }
            }
#ifndef NDEBUG
            std::cerr << "Found " << ntargetpts << " target grid points."
                      << '\n';
#endif
        }

        void sync(std::vector<std::future<void>>& results)
        {
            for(auto&& r : results)
            {
                r.get();
            }
        }

        bool solver_t::iterate()
        {
            std::vector<std::future<void>> results;
            results.reserve(n_sweeps + parallel::n_workers - 2);

            // Schedule a sweep in each worker.
            // The last sweep is done by the main thread.
            for(auto dir = 0; dir < n_sweeps - 1; ++dir)
            {
                results.emplace_back(m_pool->enqueue(detail::sweep_q,
                                                     dir,
                                                     m_queue.get(),
                                                     m_cost.get(),
                                                     m_grid.get(),
                                                     m_hamiltonian,
                                                     m_viscosity));
            }

            // Each worker enforces boundary conditions on one of the data
            // arrays.
            for(unsigned worker = 0; worker < parallel::n_workers - 1; ++worker)
            {
                results.emplace_back(m_pool->enqueue(detail::enforce_boundary_q,
                                                     m_queue.get(),
                                                     m_cost.get(),
                                                     m_grid.get()));
            }

            // Done with scheduling, main thread does work itself.
            detail::sweep(n_sweeps - 1,
                          m_worker.back().get(),
                          m_cost.get(),
                          m_grid.get(),
                          m_hamiltonian,
                          m_viscosity);

            detail::enforce_boundary(
                m_worker.back().get(), m_cost.get(), m_grid.get());

            sync(results);

            std::vector<std::future<scalar_t>> merged;
            merged.reserve(parallel::n_workers - 1);

            auto const points_per_worker = m_grid->npts() / parallel::n_workers;

            for(unsigned worker = 0; worker < parallel::n_workers - 1; ++worker)
            {
                auto const start = worker * points_per_worker;
                auto const end = start + points_per_worker;

                merged.emplace_back(
                    m_pool->enqueue(detail::merge,
                                    m_soln.get(),
                                    &m_worker,
                                    std::make_pair(start, end)));
            }

            auto diff = detail::merge(
                m_soln.get(),
                &m_worker,
                std::make_pair(
                    index_t{ (parallel::n_workers - 1) * points_per_worker },
                    index_t{ m_grid->npts() }));

            for(auto&& e : merged)
            {
                auto e_val = e.get();
                diff = std::max(diff, e_val);
            }

            return diff < m_tolerance;
        }
    }    // namespace solver
}    // namespace fsm
