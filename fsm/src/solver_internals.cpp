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

#include "solver_internals.hpp"
#include "data.hpp"
#include "grid.hpp"
#include "queue.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>

#include <iostream>

namespace fsm
{
    namespace solver
    {
        namespace detail
        {
            struct update_data_internal_t
            {
                vector_t p;
                vector_t avgs;
            };

            inline update_data_internal_t estimate_p(point_t const& point,
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
                    auto right = soln->at(grid->index(neighbor));

                    neighbor[i] -= 2;
                    auto left = soln->at(grid->index(neighbor));

                    res.p[i] = (right - left) / (scalar_t{ 2.0 } * grid->h(i));
                    res.avgs[i] =
                        (right + left) / (scalar_t{ 2.0 } * grid->h(i));
                }

                return res;
            }

            inline scalar_t update(scalar_t ham_value,
                                   scalar_t scale,
                                   scalar_t cost,
                                   vector_t const& avgs,
                                   vector_t const& viscosity)
            {
                return scale * (cost - ham_value +
                                std::inner_product(std::begin(viscosity),
                                                   std::end(viscosity),
                                                   std::begin(avgs),
                                                   0.0));
            }

            inline scalar_t scale(vector_t const& viscosity, vector_t const& h)
            {
                return scalar_t{ 1.0 } /
                       std::inner_product(std::begin(viscosity),
                                          std::end(viscosity),
                                          std::begin(h),
                                          0.0,
                                          std::plus<>(),
                                          std::divides<>());
            }

            void sweep(index_t dir,
                       queue::queue_t* queue,
                       data::data_t const* cost,
                       grid::grid_t const* grid,
                       hamiltonian_t const& hamiltonian,
                       std::function<vector_t(input_t const&)> const& viscosity)
            {
                assert(queue != nullptr);
                assert(cost != nullptr);
                assert(grid != nullptr);

                data::data_t* soln = queue->deque();

                for(auto const& i : grid::interior_visitor_t{ *grid, dir })
                {
                    if(cost->at(i) > scalar_t{ 0.0 })
                    {
                        auto point = grid->point(i);
                        auto data = estimate_p(point, soln, grid);
#ifdef FSM_USE_ROWMAJOR
                        auto sigma = viscosity(i);
                        auto scale_ = scale(sigma, grid->h());
                        soln->at(i) = std::min(update(hamiltonian(i, data.p),
                                                      scale_,
                                                      cost->at(i),
                                                      data.avgs,
                                                      sigma),
                                               soln->at(i));
#else
                        auto sigma = viscosity(point);
                        auto scale_ = scale(sigma, grid->h());
                        soln->at(i) =
                            std::min(update(hamiltonian(point, data.p),
                                            scale_,
                                            cost->at(i),
                                            data.avgs,
                                            sigma),
                                     soln->at(i));
#endif
                    }
                }

                queue->enqueue(soln);
            }

            inline scalar_t update_boundary(index_t index,
                                            index_t boundary,
                                            data::data_t const* soln,
                                            grid::grid_t const* grid)
            {
                auto neighbor = grid->point(index);
                auto boundary_dim = boundary >= dim ? boundary - dim : boundary;

                // Approximate based on the two points in a line orthogonal
                // to the boundary closest to the current point.
                neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
                auto outer = soln->at(grid->index(neighbor));
                neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
                auto inner = soln->at(grid->index(neighbor));

                return std::min(std::max(2 * outer - inner, inner),
                                soln->at(index));
            }

            void enforce_boundary(queue::queue_t* queue,
                                  data::data_t const* cost,
                                  grid::grid_t const* grid)
            {
                assert(queue != nullptr);
                assert(grid != nullptr);

                data::data_t* soln = queue->deque();

                for(index_t boundary = 0; boundary < n_boundaries; ++boundary)
                {
                    for(auto const& i :
                        grid::boundary_visitor_t{ *grid, boundary })
                    {
                        if(cost->at(i) > scalar_t{ 0.0 })
                        {
                            soln->at(i) = std::min(
                                update_boundary(i, boundary, soln, grid),
                                soln->at(i));
                        }
                    }
                }

                queue->enqueue(soln);
            }
        }    // namespace detail
    }        // namespace solver
}    // namespace fsm
