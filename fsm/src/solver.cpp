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
              m_soln(),
              m_cost(std::make_unique<data::data_t>(std::move(cost))),
              m_viscosity(viscosity),
              m_tolerance(params.tolerance),
              m_pool(std::make_unique<ThreadPool>(parallel::n_workers))
        {
            m_soln.reserve(parallel::n_workers + 1);

            for(unsigned i = 0; i < parallel::n_workers + 1; ++i)
            {
                m_soln.emplace_back(std::make_unique<data::data_t>(
                    m_grid->npts(), params.maxval));
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
            io::write(
                m_filename, "value_function", *(m_soln[0]), m_grid->size());
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
                    for(unsigned j = 0; j < m_soln.size(); ++j)
                    {
                        m_soln[j]->at(i) = -(m_cost->at(i) + scalar_t{ 1.0 });
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
            for_each(std::begin(results), std::end(results), [](auto& a) {
                a.get();
            });
        }

        bool solver_t::iterate()
        {
            auto n_sweeps_done = 0;

            while(n_sweeps_done < n_sweeps)
            {
                auto n_workers_launched =
                    std::min(static_cast<unsigned>(parallel::n_workers),
                             static_cast<unsigned>(n_sweeps - n_sweeps_done)) -
                    1;

                std::vector<std::future<void>> results;
                results.reserve(n_workers_launched);

                unsigned worker = 0;

                for(; worker < n_workers_launched; ++worker)
                {
                    results.emplace_back(
                        m_pool->enqueue(detail::sweep,
                                        n_sweeps_done++,
                                        m_soln[1 + worker].get(),
                                        m_cost.get(),
                                        m_grid.get(),
                                        m_hamiltonian,
                                        m_viscosity));
                }

                detail::sweep(n_sweeps_done++,
                              m_soln[1 + worker].get(),
                              m_cost.get(),
                              m_grid.get(),
                              m_hamiltonian,
                              m_viscosity);

                sync(results);
            }

            auto n_boundaries_done = 0;

            while(n_boundaries_done < n_boundaries)
            {
                auto n_workers_launched =
                    std::min(static_cast<unsigned>(parallel::n_workers),
                             static_cast<unsigned>(n_boundaries -
                                                   n_boundaries_done)) -
                    1;

                std::vector<std::future<void>> results;
                results.reserve(n_workers_launched);

                unsigned worker = 0;

                for(; worker < n_workers_launched; ++worker)
                {
                    results.emplace_back(
                        m_pool->enqueue(detail::enforce_boundary,
                                        n_boundaries_done++,
                                        m_soln[1 + worker].get(),
                                        m_cost.get(),
                                        m_grid.get()));
                }

                detail::enforce_boundary(n_boundaries_done++,
                                         m_soln[1 + worker].get(),
                                         m_cost.get(),
                                         m_grid.get());
                sync(results);
            }

            return merge() < m_tolerance;
        }

        scalar_t solver_t::merge()
        {
            auto diff = scalar_t{ 0 };

            for(index_t i = 0; i < m_grid->npts(); ++i)
            {
                auto old = m_soln[0]->at(i);

                m_soln[0]->at(i) =
                    std::min_element(std::begin(m_soln),
                                     std::end(m_soln),
                                     [&](auto const& a, auto const& b) {
                                         return a->at(i) < b->at(i);
                                     })
                        ->get()
                        ->at(i);

                for(unsigned j = 1; j < m_soln.size(); ++j)
                {
                    m_soln[j]->at(i) = m_soln[0]->at(i);
                }

                diff = std::max(diff, old - m_soln[0]->at(i));
            }

            return diff;
        }
    }    // namespace solver
}    // namespace fsm
