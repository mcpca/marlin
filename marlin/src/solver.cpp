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

#include <algorithm>
#include <numeric>

#include "marlin/solver.hpp"

#include "levels.hpp"

namespace marlin
{
    namespace solver
    {
        solver_t::solver_t(
            std::vector<scalar_t>&& rhs,
            point_t const& dimensions,
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            params_t const& params) noexcept
            : m_grid(vertices, dimensions),
              m_soln(m_grid.npts(), params.maxval),
              m_cost(std::move(rhs)),
              m_tolerance(params.tolerance)
        {
            compute_levels();
            compute_bdry_idxs();
            initialize();
        }

        std::vector<scalar_t>&& solver_t::steal() noexcept
        {
            return m_soln.steal();
        }

        void solver_t::compute_levels() noexcept
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

        void solver_t::compute_bdry_idxs() noexcept
        {
            for(int bd = 0; bd < n_boundaries; ++bd)
            {
                int const boundary_dim = bd >= dim ? bd - dim : bd;

                auto size = m_grid.size();
                size[boundary_dim] = index_t{ 1 };

                m_bdry_idxs[bd].reserve(std::accumulate(std::cbegin(size),
                                                        std::cend(size),
                                                        index_t{ 1 },
                                                        std::multiplies<>()));

                auto const npts = m_grid.npts();

                for(auto i = m_grid.next_in_boundary(npts, boundary_dim);
                    i != npts;
                    i = m_grid.next_in_boundary(i, boundary_dim))
                {
                    auto index = i;

                    if(bd >= dim)
                    {
                        point_t point = m_grid.point(i);
                        point[boundary_dim] = m_grid.size(boundary_dim) - 1;
                        index = m_grid.index(point);
                    }

                    if(m_cost.at(index) > scalar_t{ 0.0 })
                        m_bdry_idxs[bd].emplace_back(index);
                }
            }
        }

        void solver_t::initialize() noexcept
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

        static inline scalar_t update_boundary(index_t index,
                                               index_t boundary,
                                               data_t const& soln,
                                               grid_t const& grid) noexcept
        {
            assert(boundary < n_boundaries);
            assert(index < grid.npts());

            auto neighbor = grid.point(index);
            auto const boundary_dim =
                boundary >= dim ? boundary - dim : boundary;

            // Approximate based on the two points in a line orthogonal
            // to the boundary closest to the current point.
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto const outer = soln.at(grid.index(neighbor));
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto const inner = soln.at(grid.index(neighbor));

            return std::min(std::max(2 * outer - inner, inner), soln.at(index));
        }

        scalar_t solver_t::boundary_sweep(index_t boundary) noexcept
        {
            assert(boundary < n_boundaries);

            auto const& bdry_idxs = m_bdry_idxs[boundary];

            std::vector<scalar_t> deltas(bdry_idxs.size());

#pragma omp parallel default(none) shared(boundary, bdry_idxs, deltas)
#pragma omp for schedule(static) nowait
            for(size_t i = 0; i < bdry_idxs.size(); ++i)
            {
                index_t const index = bdry_idxs[i];
                scalar_t const old = m_soln.at(index);

                scalar_t const new_val = m_soln.at(index) = std::min(
                    update_boundary(index, boundary, m_soln, m_grid), old);

                deltas[i] = old - new_val;
            }

            return *std::max_element(std::cbegin(deltas), std::cend(deltas));
        }

        scalar_t solver_t::boundary() noexcept
        {
            scalar_t diff = 0;

            for(index_t bdry = 0; bdry < n_boundaries; ++bdry)
            {
                diff = std::max(diff, boundary_sweep(bdry));
            }

            return diff;
        }
    }    // namespace solver
}    // namespace marlin
