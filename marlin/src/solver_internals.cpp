////////////////////////////////////////////////////////////////////////////////
// MIT License
//
// Copyright (c) 2020 Miguel Aguiar
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

#include "data.hpp"
#include "grid.hpp"
#include "marlin/solver.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>

#include "ThreadPool/ThreadPool.h"

namespace marlin
{
    namespace solver
    {
        namespace
        {
            struct update_data_internal_t
            {
                vector_t p;
                vector_t avgs;
            };
        }    // namespace

        static inline update_data_internal_t estimate_p(
            point_t const& point,
            data::data_t const* soln,
            grid::grid_t const* grid)
        {
            assert(soln != nullptr);
            assert(grid != nullptr);

            update_data_internal_t res;

            for(auto i = 0; i < dim; ++i)
            {
                auto neighbor = point;

                neighbor[i] += 1;
                auto const right = soln->at(grid->index(neighbor));

                neighbor[i] -= 2;
                auto const left = soln->at(grid->index(neighbor));

                res.p[i] = (right - left) / (scalar_t{ 2.0 } * grid->h(i));
                res.avgs[i] = (right + left) / (scalar_t{ 2.0 } * grid->h(i));
            }

            return res;
        }

        static inline scalar_t update(scalar_t ham_value,
                                      scalar_t scale,
                                      scalar_t cost,
                                      vector_t const& avgs,
                                      vector_t const& viscosity)
        {
            assert(cost > 0);
            return scale * (cost - ham_value +
                            std::inner_product(std::begin(viscosity),
                                               std::end(viscosity),
                                               std::begin(avgs),
                                               0.0));
        }

        static inline scalar_t scale(vector_t const& viscosity,
                                     vector_t const& h)
        {
            return scalar_t{ 1.0 } / std::inner_product(std::begin(viscosity),
                                                        std::end(viscosity),
                                                        std::begin(h),
                                                        0.0,
                                                        std::plus<>(),
                                                        std::divides<>());
        }

        scalar_t solver_t::update_points(std::vector<point_t> const* points,
                                         int dir,
                                         int start,
                                         int end)
        {
            assert(points != nullptr);
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
                auto point = m_grid->rotate_axes((*points)[j], dir);
                auto const index = m_grid->index(point);

                if(m_cost->at(index) < 0.0)
                {
                    continue;
                }

                auto const data = estimate_p(point, m_soln.get(), m_grid.get());
                auto const old = m_soln->at(index);

#ifdef MARLIN_USE_ROWMAJOR
                auto const sigma = m_viscosity(index);
                auto const scale_ = scale(sigma, m_grid->h());

                m_soln->at(index) =
                    std::min(update(m_hamiltonian(index, data.p),
                                    scale_,
                                    m_cost->at(index),
                                    data.avgs,
                                    sigma),
                             old);
#else
                auto const sigma = m_viscosity(point);
                auto const scale_ = scale(sigma, m_grid->h());

                m_soln->at(index) =
                    std::min(update(m_hamiltonian(point, data.p),
                                    scale_,
                                    m_cost->at(index),
                                    data.avgs,
                                    sigma),
                             old);
#endif

                diff = std::max(diff, old - m_soln->at(index));
            }

            return diff;
        }

        scalar_t solver_t::sweep(int dir)
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
                    results.emplace_back(
                        m_pool->enqueue(&solver_t::update_points,
                                        this,
                                        &level,
                                        dir,
                                        start,
                                        start + block_size));
                    start += block_size;
                }

                // Update remaining points
                diff = std::max(
                    diff, update_points(&level, dir, start, level.size()));

                for(auto& r : results)
                {
                    diff = std::max(diff, r.get());
                }
            }

            return diff;
        }

        static inline scalar_t update_boundary(index_t index,
                                               index_t boundary,
                                               data::data_t const* soln,
                                               grid::grid_t const* grid)
        {
            assert(boundary < n_boundaries);
            assert(soln != nullptr);
            assert(grid != nullptr);
            assert(index < grid->npts());

            auto neighbor = grid->point(index);
            auto const boundary_dim =
                boundary >= dim ? boundary - dim : boundary;

            // Approximate based on the two points in a line orthogonal
            // to the boundary closest to the current point.
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto const outer = soln->at(grid->index(neighbor));
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto const inner = soln->at(grid->index(neighbor));

            return std::min(std::max(2 * outer - inner, inner),
                            soln->at(index));
        }

        scalar_t solver_t::boundary()
        {
            scalar_t diff = 0;
            auto const size = m_grid->npts();

            for(index_t bdry = 0; bdry < n_boundaries; ++bdry)
            {
                for(auto i = m_grid->next_in_boundary(size, bdry); i != size;
                    i = m_grid->next_in_boundary(i, bdry))
                {
                    if(m_cost->at(i) > scalar_t{ 0.0 })
                    {
                        auto old = m_soln->at(i);
                        m_soln->at(i) =
                            std::min(update_boundary(
                                         i, bdry, m_soln.get(), m_grid.get()),
                                     old);
                        diff = std::max(diff, old - m_soln->at(i));
                    }
                }
            }
            return diff;
        }
    }    // namespace solver
}    // namespace marlin
