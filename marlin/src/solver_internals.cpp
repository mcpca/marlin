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

#include "marlin/solver.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>

namespace marlin
{
    namespace solver
    {
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

        scalar_t solver_t::boundary() noexcept
        {
            scalar_t diff = 0;
            auto const size = m_grid.npts();

            for(index_t bdry = 0; bdry < n_boundaries; ++bdry)
            {
                for(auto i = m_grid.next_in_boundary(size, bdry); i != size;
                    i = m_grid.next_in_boundary(i, bdry))
                {
                    if(m_cost.at(i) > scalar_t{ 0.0 })
                    {
                        auto old = m_soln.at(i);
                        m_soln.at(i) = std::min(
                            update_boundary(i, bdry, m_soln, m_grid), old);
                        diff = std::max(diff, old - m_soln.at(i));
                    }
                }
            }
            return diff;
        }
    }    // namespace solver
}    // namespace marlin
