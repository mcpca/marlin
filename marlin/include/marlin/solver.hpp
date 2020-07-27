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
#include <string>
#include <vector>

#include "defs.hpp"

#include "data.hpp"
#include "grid.hpp"

#if defined(NDEBUG) and not defined(MARLIN_PRINT_DEBUG_MSGS)
#    define MARLIN_DEBUG(x) ;
#else
#    include <iostream>
#    define MARLIN_DEBUG(x) x
#endif

namespace marlin
{
    using grid_t = grid::grid_t<dim, index_t, scalar_t>;
    using data_t = data::data_t<index_t, scalar_t>;

    namespace solver
    {
        //! @brief Implemented boundary conditions types.
        enum class boundary_condition_t
        {
            //! Extrapolate value of solution.
            extrapolate,

            //! Periodic boundary condition
            //! The two boundaries along the dimension are identified as the
            //! same point, and updated as if they were a regular gridpoint.
            periodic,

            //! Boundary is never updated.
            none,
        };

        //! @brief Numerical parameters supplied by the user.
        struct params_t
        {
            //! Maximum value of the solution over the whole domain.
            scalar_t maxval;
            //! Tolerance parameter for the convergence criterion.
            scalar_t tolerance;
        };

        namespace detail
        {
            template<index_t N, typename T>
            constexpr std::array<T, N> make_array(T init)
            {
                std::array<T, N> a;

                for(index_t i = 0; i < N; ++i)
                {
                    a[i] = init;
                }

                return a;
            }
        }    // namespace detail

        //! @brief Solver API.
        //
        //! The main class, used for defining and solving a problem.
        class solver_t
        {
          public:
            //! Constructs the solver given data for the right-hand-side of the
            //! equation (\p rhs), the grid sizes and bounds, the initial valu
            //! of the solution, the convergence parameter, and optionally the
            //! boundary condition types.
            //
            //! @param rhs Values of the right-hand side of the equation over
            //! the grid, stored in row-major order.
            //! @param dimensions Grid dimensions.
            //! @param vertices The limits of the grid.
            //! @param params Numerical parameters.
            //! @param boundary_conditions Boundary condition used for each
            //! piece of the computational boundary.
            solver_t(
                std::vector<scalar_t>&& rhs,
                point_t const& dimensions,
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                params_t const& params,
                std::array<boundary_condition_t, n_boundaries> const&
                    boundary_conditions = detail::make_array<n_boundaries>(
                        boundary_condition_t::extrapolate));

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
            void compute_bdry_idxs();

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
            template<typename Hamiltonian, typename Viscosity>
            scalar_t boundary(Hamiltonian const& hamiltonian,
                              Viscosity const& viscosity) noexcept;

            // Updates a single boundary via extrapolation.
            scalar_t boundary_sweep_extrapolate(index_t boundary) noexcept;

            // Updates a boundary with a periodic boundary condition.
            template<typename Hamiltonian, typename Viscosity>
            scalar_t boundary_sweep_periodic(
                index_t boundary,
                Hamiltonian const& hamiltonian,
                Viscosity const& viscosity) noexcept;

            // Updates one point.
            template<typename Hamiltonian,
                     typename Viscosity,
                     typename GradientEstimator>
            scalar_t update_point(int dir,
                                  point_t point,
                                  Hamiltonian const& hamiltonian,
                                  Viscosity const& viscosity,
                                  GradientEstimator const& ge) noexcept;

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

            // Boundary condition used in each piece of the grid boundary.
            std::array<boundary_condition_t, n_boundaries> m_bc;
        };

        namespace detail
        {
            struct update_data_internal_t
            {
                vector_t p;
                vector_t avgs;
            };

            inline update_data_internal_t estimate_p_interior(
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

            inline update_data_internal_t estimate_p_boundary(
                point_t const& point,
                data_t const& soln,
                grid_t const& grid) noexcept
            {
                update_data_internal_t res;

                for(auto i = 0; i < dim; ++i)
                {
                    auto neighbor = point;
                    auto const original_idx = point[i];

                    neighbor[i] += 1;

                    auto const right = soln.at(grid.index(neighbor));

                    if(original_idx > 0)
                        neighbor[i] = original_idx - 1;
                    else
                        neighbor[i] = grid.size(i) - 2;

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
            MARLIN_DEBUG(auto niter = 0;
                         std::cerr << "Iteration " << niter++ << ":\n";)

            while(!iterate(hamiltonian, viscosity))
            {
                MARLIN_DEBUG(std::cerr << "Iteration " << niter++ << ":\n";)
            }
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
                diff = std::max(diff, boundary(hamiltonian, viscosity));
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
                std::size_t const size = level.size();
                std::vector<scalar_t> delta(size);

#pragma omp parallel default(none) \
    shared(delta, level, dir, hamiltonian, viscosity, diff, size)
                {
#pragma omp for schedule(static) nowait
                    for(auto i = 0ul; i < size; ++i)
                    {
                        delta[i] = update_point(dir,
                                                level[i],
                                                hamiltonian,
                                                viscosity,
                                                detail::estimate_p_interior);
                    }

#pragma omp for schedule(static) reduction(max : diff) nowait
                    for(std::size_t i = 0; i < size; ++i)
                        diff = std::max(delta[i], diff);
                }
            }

            return diff;
        }

        template<typename Hamiltonian,
                 typename Viscosity,
                 typename GradientEstimator>
        scalar_t solver_t::update_point(int dir,
                                        point_t point,
                                        Hamiltonian const& hamiltonian,
                                        Viscosity const& viscosity,
                                        GradientEstimator const& ge) noexcept
        {
            point = m_grid.rotate_axes(point, dir);
            index_t const index = m_grid.index(point);
            scalar_t const cost = m_cost.at(index);

            if(cost < scalar_t{ 0.0 })
            {
                return scalar_t{ 0.0 };
            }

            auto const data = ge(point, m_soln, m_grid);

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

        template<typename Hamiltonian, typename Viscosity>
        scalar_t solver_t::boundary(Hamiltonian const& hamiltonian,
                                    Viscosity const& viscosity) noexcept
        {
            scalar_t diff = 0;

            for(index_t bdry = 0; bdry < n_boundaries; ++bdry)
            {
                switch(m_bc[bdry])
                {
                    case boundary_condition_t::extrapolate:
                        diff = std::max(diff, boundary_sweep_extrapolate(bdry));
                        break;
                    case boundary_condition_t::periodic:
                        diff = std::max(diff,
                                        boundary_sweep_periodic(
                                            bdry, hamiltonian, viscosity));
                        break;
                    case boundary_condition_t::none:
                        break;
                }
            }

            return diff;
        }

        // Updates a boundary with a periodic boundary condition.
        template<typename Hamiltonian, typename Viscosity>
        scalar_t solver_t::boundary_sweep_periodic(
            index_t boundary,
            Hamiltonian const& hamiltonian,
            Viscosity const& viscosity) noexcept
        {
            assert(boundary < n_boundaries);

            scalar_t delta = 0.0;

            for(index_t index : m_bdry_idxs[boundary])
            {
                auto const point = m_grid.point(index);

                if(boundary < dim)
                {
                    delta = std::max(delta,
                                     update_point(0,
                                                  point,
                                                  hamiltonian,
                                                  viscosity,
                                                  detail::estimate_p_boundary));
                }
                else
                {
                    auto opposing_point = point;
                    opposing_point[boundary - dim] = 0;
                    m_soln.at(index) = m_soln.at(m_grid.index(opposing_point));

                    // We don't need to keep track of delta in this case, since
                    // it will be the same as for the opposing boundary.
                }
            }

            return delta;
        }
    }    // namespace solver
}    // namespace marlin
