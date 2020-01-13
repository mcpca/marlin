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

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

#include "defs.hpp"

#include "data.hpp"
#include "grid.hpp"

#if defined(NDEBUG) and not defined(MARLIN_PRINT_DEBUG_MSGS)
#    define MARLIN_DEBUG(x) ;
#else
#    define MARLIN_DEBUG(x) x
#endif

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
            //! @param rhs Values of the right-hand side of the equation over
            //! the grid, stored in row-major order.
            //! @param dimensions Grid dimensions.
            //! @param vertices The limits of the grid.
            //! @param params Numerical parameters.
            solver_t(
                std::vector<scalar_t>&& rhs,
                point_t const& dimensions,
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                params_t const& params) noexcept;

            //! Compiler-generated move contructor.
            solver_t(solver_t&&) noexcept = default;
            //! Compiler-generated move assignment.
            solver_t& operator=(solver_t&&) noexcept = default;
            //! Compiler-generated destructor.
            ~solver_t() = default;

            //! Solves the problem instance.
            //
            //! @param hamiltonian Callable for evaluating the Hamiltonian.
            //! @param viscosity Callable for evaluating the viscosity
            //!                  coefficients.
            template<typename Hamiltonian, typename Viscosity>
            void solve(Hamiltonian const& hamiltonian,
                       Viscosity const& viscosity) noexcept;

            //! Extract the solution.
            //
            //! This moves the solution out of the solver object, so the solver
            //! object should not be used after calling this function.
            std::vector<scalar_t>&& steal() noexcept;

          private:
            // Initialize the solution.
            void initialize() noexcept;

            // Fill m_levels.
            void compute_levels() noexcept;

            // Fill m_bdry_idxs.
            void compute_bdry_idxs() noexcept;

            // Main loop.
            template<typename Hamiltonian, typename Viscosity>
            bool iterate(Hamiltonian const& hamiltonian,
                         Viscosity const& viscosity) noexcept;

            // Sweeps all gridpoints in the direction dir.
            template<typename Hamiltonian, typename Viscosity>
            scalar_t sweep(int dir,
                           Hamiltonian const& hamiltonian,
                           Viscosity const& viscosity) noexcept;

            // Updates the boundary points.
            scalar_t boundary() noexcept;

            // Updates a single boundary.
            scalar_t boundary_sweep(index_t boundary) noexcept;

            // Updates one point.
            template<typename Hamiltonian, typename Viscosity>
            scalar_t update_point(int dir,
                                  point_t point,
                                  Hamiltonian const& hamiltonian,
                                  Viscosity const& viscosity) noexcept;

            grid_t m_grid;    // Grid.
            data_t m_soln;    // Solution.
            data_t m_cost;    // Cost function.

            // Tolerance for the convergence criterion.
            scalar_t m_tolerance;

            // Caches the sets of points which can be updated in parallel in
            // each sweep.
            std::vector<std::vector<point_t>> m_levels;

            // Caches the indexes of the points in each boundary.
            std::array<std::vector<index_t>, n_boundaries> m_bdry_idxs;
        };

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
                             Viscosity const& viscosity) noexcept
        {
            std::cout << "Solving (" << omp_get_max_threads() << " threads)..."
                      << std::endl;

            MARLIN_DEBUG(auto niter = 0;
                         std::cerr << "Iteration " << niter++ << ":\n";)

            while(!iterate(hamiltonian, viscosity))
            {
                MARLIN_DEBUG(std::cerr << "Iteration " << niter++ << ":\n";)
            }

            std::cout << "Done." << std::endl;
        }

        template<typename Hamiltonian, typename Viscosity>
        bool solver_t::iterate(Hamiltonian const& hamiltonian,
                               Viscosity const& viscosity) noexcept
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
                std::vector<scalar_t> delta(level.size());

#pragma omp parallel default(none) \
    shared(delta, level, dir, hamiltonian, viscosity)
                {
#pragma omp for schedule(static) nowait
                    for(auto i = 0ul; i < level.size(); ++i)
                    {
                        delta[i] =
                            update_point(dir, level[i], hamiltonian, viscosity);
                    }
                }

                diff = std::max(
                    diff,
                    *std::max_element(std::cbegin(delta), std::cend(delta)));
            }

            return diff;
        }

        template<typename Hamiltonian, typename Viscosity>
        scalar_t solver_t::update_point(int dir,
                                        point_t point,
                                        Hamiltonian const& hamiltonian,
                                        Viscosity const& viscosity) noexcept
        {
            assert(dir >= 0);
            assert(dir < n_sweeps);

            point = m_grid.rotate_axes(point, dir);
            index_t const index = m_grid.index(point);
            scalar_t const cost = m_cost.at(index);

            if(cost < scalar_t{ 0.0 })
            {
                return scalar_t{ 0.0 };
            }

            auto const data = detail::estimate_p(point, m_soln, m_grid);
            scalar_t const old = m_soln.at(index);

#ifdef MARLIN_USE_ROWMAJOR
            vector_t const sigma = viscosity(index);
            scalar_t const hval = hamiltonian(index, data.p);
#else
            vector_t const sigma = viscosity(point);
            scalar_t const hval = hamiltonian(point, data.p);
#endif
            scalar_t const scale_ = detail::scale(sigma, m_grid.h());

            scalar_t const new_val = m_soln.at(index) = std::min(
                detail::update(hval, scale_, cost, data.avgs, sigma), old);

            return old - new_val;
        }
    }    // namespace solver
}    // namespace marlin
