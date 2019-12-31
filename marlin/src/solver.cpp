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
// https://github.com/mcpca/marlin

#if defined(NDEBUG) and not defined(MARLIN_PRINT_DEBUG_MSGS)
#    define MARLIN_DEBUG(x) ;
#else
#    define MARLIN_DEBUG(x) x
#endif

#include <algorithm>
#include <iostream>
#include <numeric>

#include "marlin/solver.hpp"

#include "io.hpp"
#include "levels.hpp"

#include "ThreadPool/ThreadPool.h"

namespace marlin
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
            auto data = io::read(filename, "cost_function");

            if(io::dset_exists(filename, "value_function"))
            {
                throw std::runtime_error(
                    "A dataset named \'value_function\' already exists in the "
                    "given file.");
            }

            grid_t g(vertices, data.size);
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
            data_t cost,
            grid_t const& grid,
            std::function<vector_t(input_t const&)> const& viscosity,
            params_t const& params)
            : m_filename(filename),
              m_hamiltonian(hamiltonian),
              m_grid(grid),
              m_soln(m_grid.npts(), params.maxval),
              m_cost(std::move(cost)),
              m_viscosity(viscosity),
              m_tolerance(params.tolerance),
              m_pool(std::make_unique<ThreadPool>(n_workers - 1))
        {
            for(index_t i = 0; i < m_grid.n_levels(); ++i)
            {
                level::level_t<dim> level(i, m_grid.size());

                std::vector<point_t> points;

                do
                {
                    point_t point;
                    level.get(point.data());

                    if(!m_grid.is_boundary(point))
                    {
                        points.push_back(point);
                    }
                } while(level.next());

                if(!points.empty())
                {
                    m_levels.emplace_back(std::move(points));
                }
            }
        }

        void solver_t::solve()
        {
            std::cout << "Initializing problem instance..." << std::endl;
            initialize();
            std::cout << "Initialized. Solving (" << n_workers << " threads)..."
                      << std::endl;

            MARLIN_DEBUG(auto niter = 0;
                         std::cerr << "Iteration " << niter++ << ":\n";)

            while(!iterate())
            {
                MARLIN_DEBUG(std::cerr << "Iteration " << niter++ << ":\n";)
            }

            std::cout << "Done. Writing to " << m_filename << "..."
                      << std::endl;
            io::write(m_filename, "value_function", m_soln, m_grid.size());
        }

        void solver_t::initialize()
        {
            MARLIN_DEBUG(auto ntargetpts = 0;)

            for(index_t i = 0; i < m_cost.size(); ++i)
            {
                if(m_cost.at(i) < scalar_t{ 0.0 })
                {
                    MARLIN_DEBUG(++ntargetpts;)

                    m_soln.at(i) = -(m_cost.at(i) + scalar_t{ 1.0 });
                }
            }

            MARLIN_DEBUG(std::cerr << "Found " << ntargetpts
                                   << " target grid points." << '\n';)
        }

        bool solver_t::iterate()
        {
            for(auto dir = 0; dir < n_sweeps; ++dir)
            {
                scalar_t diff = 0;

                diff = sweep(dir);
                MARLIN_DEBUG(std::cerr << "Sweep " << dir
                                       << ": delta = " << diff << '\n';)
                diff = std::max(diff, boundary());
                MARLIN_DEBUG(std::cerr << "Sweep " << dir
                                       << " (after boundary): delta = " << diff
                                       << '\n';)
                if(diff < m_tolerance)
                {
                    return true;
                }
            }

            return false;
        }
    }    // namespace solver
}    // namespace marlin
