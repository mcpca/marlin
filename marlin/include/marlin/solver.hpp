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

#pragma once

#include <array>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "defs.hpp"

#include "data.hpp"
#include "grid.hpp"

#include "ThreadPool/ThreadPool.h"

#if defined(NDEBUG) and not defined(MARLIN_PRINT_DEBUG_MSGS)
#    define MARLIN_DEBUG(x) ;
#else
#    define MARLIN_DEBUG(x) x
#endif

class ThreadPool;

namespace marlin
{
    using grid_t = grid::grid_t<dim, index_t, scalar_t>;
    using data_t = data::data_t<index_t, scalar_t>;

    namespace solver
    {
        //! @brief Numerical parameters supplied by the user.
        struct params_t
        {
            //! Maximum value of the solution over the whole domain.
            scalar_t maxval;
            //! Tolerance parameter for the convergence criterion.
            scalar_t tolerance;
        };

        //! @brief Solver API.
        //
        //! The main class, used for defining a problem, solving it and
        //! writing the solution to disk.
        class solver_t
        {
          public:
            //! Contructs a solver_t object given an HDF5 file and the problem
            //! data.
            //! Construction is delegated to the move constructor by calling a
            //! factory function, which makes it easier to read the data in the
            //! HDF5 file before constructing the solver_t object.
            //
            //! @param filename Path to a HDF5 containing the cost function.
            //! @param hamiltonian Callable for evaluating the Hamiltonian.
            //! @param vertices The limits of the grid.
            //! @param viscosity Callable for evaluating the viscosity
            //!                  coefficients.
            //! @param params Numerical parameters.
            solver_t(
                std::string const& filename,
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                params_t const& params);

            //! Compiler-generated move contructor.
            solver_t(solver_t&&) noexcept;
            //! Compiler-generated move assignment.
            solver_t& operator=(solver_t&&) noexcept;
            //! Compiler-generated destructor.
            ~solver_t();

            //! Initializes and solves the problem instance, and writes the
            //! solution to disk.
            template<typename Hamiltonian, typename Viscosity>
            void solve(Hamiltonian const& hamiltonian,
                       Viscosity const& viscosity);

          private:
            //! @brief Factory function.
            //
            //! Reads and processes the necessary data from the hdf5 file before
            //! calling a private constructor.
            friend solver_t make_solver(
                std::string const& filename,
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                params_t const& params);

            // Constructs a solver_t object after the relevant info has been
            // read from the HDF5 file.
            solver_t(std::string const& filename,
                     data_t cost,
                     grid_t const& grid,
                     params_t const& params);

            // Initialize the solution.
            void initialize();

            // Main loop.
            template<typename Hamiltonian, typename Viscosity>
            bool iterate(Hamiltonian const& hamiltonian,
                         Viscosity const& viscosity);

            void write() const;

            // Sweeps all gridpoints in the direction dir.
            template<typename Hamiltonian, typename Viscosity>
            scalar_t sweep(int dir,
                           Hamiltonian const& hamiltonian,
                           Viscosity const& viscosity) noexcept;

            // Updates the boundary points.
            scalar_t boundary() noexcept;

            // Updates a group of points.
            template<typename Hamiltonian, typename Viscosity>
            scalar_t update_points(int dir,
                                   int start,
                                   int end,
                                   std::vector<point_t> const& points,
                                   Hamiltonian const& hamiltonian,
                                   Viscosity const& viscosity) noexcept;

            // Path to a HDF5 file containing the cost function.
            std::string m_filename;

            grid_t m_grid;    // Grid.
            data_t m_soln;    // Solution.
            data_t m_cost;    // Cost function.

            // Tolerance for the convergence criterion.
            scalar_t m_tolerance;

            // Thread pool.
            std::unique_ptr<ThreadPool> m_pool;

            // Caches the sets of points which can be updated in parallel in
            // each sweep.
            std::vector<std::vector<point_t>> m_levels;
        };    // namespace detailclasssolver_t

        namespace detail
        {
            struct update_data_internal_t
            {
                vector_t p;
                vector_t avgs;
            };

            inline update_data_internal_t estimate_p(
                point_t const& point,
                data_t const& soln,
                grid_t const& grid) noexcept
            {
                assert(soln != nullptr);

                update_data_internal_t res;

                for(auto i = 0; i < dim; ++i)
                {
                    auto neighbor = point;

                    neighbor[i] += 1;
                    auto const right = soln.at(grid.index(neighbor));

                    neighbor[i] -= 2;
                    auto const left = soln.at(grid.index(neighbor));

                    res.p[i] = (right - left) / (scalar_t{ 2.0 } * grid.h(i));
                    res.avgs[i] =
                        (right + left) / (scalar_t{ 2.0 } * grid.h(i));
                }

                return res;
            }

            inline scalar_t update(scalar_t ham_value,
                                   scalar_t scale,
                                   scalar_t cost,
                                   vector_t const& avgs,
                                   vector_t const& viscosity) noexcept
            {
                assert(cost > 0);
                return scale * (cost - ham_value +
                                std::inner_product(std::begin(viscosity),
                                                   std::end(viscosity),
                                                   std::begin(avgs),
                                                   0.0));
            }

            inline scalar_t scale(vector_t const& viscosity,
                                         vector_t const& h) noexcept
            {
                return scalar_t{ 1.0 } /
                       std::inner_product(std::begin(viscosity),
                                          std::end(viscosity),
                                          std::begin(h),
                                          0.0,
                                          std::plus<>(),
                                          std::divides<>());
            }
        }    // namespace detail

        template<typename Hamiltonian, typename Viscosity>
        void solver_t::solve(Hamiltonian const& hamiltonian,
                             Viscosity const& viscosity)
        {
            std::cout << "Initializing problem instance..." << std::endl;
            initialize();
            std::cout << "Initialized. Solving (" << n_workers << " threads)..."
                      << std::endl;

            MARLIN_DEBUG(auto niter = 0;
                         std::cerr << "Iteration " << niter++ << ":\n";)

            while(!iterate(hamiltonian, viscosity))
            {
                MARLIN_DEBUG(std::cerr << "Iteration " << niter++ << ":\n";)
            }

            std::cout << "Done. Writing to " << m_filename << "..."
                      << std::endl;
            write();
        }

        template<typename Hamiltonian, typename Viscosity>
        bool solver_t::iterate(Hamiltonian const& hamiltonian,
                               Viscosity const& viscosity)
        {
            for(auto dir = 0; dir < n_sweeps; ++dir)
            {
                scalar_t diff = 0;

                diff = sweep(dir, hamiltonian, viscosity);
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

        template<typename Hamiltonian, typename Viscosity>
        scalar_t solver_t::sweep(int dir,
                                 Hamiltonian const& hamiltonian,
                                 Viscosity const& viscosity) noexcept
        {
            assert(dir >= 0);
            assert(dir < n_sweeps);

            scalar_t diff = 0;

            for(auto const& level : m_levels)
            {
                std::vector<std::future<scalar_t>> results;
                results.reserve(n_workers - 1);

                auto block_size = level.size() / n_workers;
                auto start = 0;

                for(auto i = 0; i < n_workers - 1; ++i)
                {
                    results.emplace_back(m_pool->enqueue([this,
                                                          dir,
                                                          start,
                                                          block_size,
                                                          &level,
                                                          &hamiltonian,
                                                          &viscosity] {
                        return update_points(dir,
                                             start,
                                             start + block_size,
                                             level,
                                             hamiltonian,
                                             viscosity);
                    }));
                    start += block_size;
                }

                // Update remaining points
                diff = std::max(diff,
                                update_points(dir,
                                              start,
                                              level.size(),
                                              level,
                                              hamiltonian,
                                              viscosity));

                for(auto& r : results)
                {
                    diff = std::max(diff, r.get());
                }
            }

            return diff;
        }

        template<typename Hamiltonian, typename Viscosity>
        scalar_t solver_t::update_points(int dir,
                                         int start,
                                         int end,
                                         std::vector<point_t> const& points,
                                         Hamiltonian const& hamiltonian,
                                         Viscosity const& viscosity) noexcept
        {
            assert(dir >= 0);
            assert(dir < n_sweeps);
            assert(start >= 0);
            assert(static_cast<size_t>(start) < points->size());

            assert(end >= 0);
            assert(static_cast<size_t>(end) <= points->size());
            assert(start <= end);

            scalar_t diff = 0;

            for(auto j = start; j < end; ++j)
            {
                auto point = m_grid.rotate_axes(points[j], dir);
                auto const index = m_grid.index(point);

                if(m_cost.at(index) < 0.0)
                {
                    continue;
                }

                auto const data = detail::estimate_p(point, m_soln, m_grid);
                auto const old = m_soln.at(index);

#ifdef MARLIN_USE_ROWMAJOR
                auto const sigma = viscosity(index);
                auto const scale_ = detail::scale(sigma, m_grid.h());

                m_soln.at(index) =
                    std::min(detail::update(hamiltonian(index, data.p),
                                            scale_,
                                            m_cost.at(index),
                                            data.avgs,
                                            sigma),
                             old);
#else
                auto const sigma = viscosity(point);
                auto const scale_ = detail::scale(sigma, m_grid.h());

                m_soln.at(index) =
                    std::min(detail::update(hamiltonian(point, data.p),
                                            scale_,
                                            m_cost.at(index),
                                            data.avgs,
                                            sigma),
                             old);
#endif

                diff = std::max(diff, old - m_soln.at(index));
            }

            return diff;
        }
    }    // namespace solver
}    // namespace marlin
