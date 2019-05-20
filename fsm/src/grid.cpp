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
#include <cassert>
#include <numeric>

#include "grid.hpp"

namespace fsm
{
    namespace grid
    {
        grid_t::grid_t(
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            point_t const& size)
            : m_size(size),
              m_npts(std::accumulate(std::begin(m_size),
                                     std::end(m_size),
                                     1,
                                     std::multiplies<>())),
              m_h([&] {
                  vector_t h;

                  for(auto i = 0; i < dim; ++i)
                  {
                      h[i] = (vertices[i].second - vertices[i].first) /
                             (m_size[i] - 1);
                  }

                  return h;
              }())
        {
            assert([this] {
                for(auto i = 0; i < dim; ++i)
                {
                    if(m_h[i] <= 0)
                        return false;
                }
                return true;
            }());
        }

        index_t grid_t::size(index_t i) const
        {
            assert(i < dim);
            return m_size[i];
        }

        point_t const& grid_t::size() const { return m_size; }

        index_t grid_t::npts() const { return m_npts; }

        scalar_t grid_t::h(index_t i) const
        {
            assert(i < dim);
            return m_h[i];
        }

        vector_t const& grid_t::h() const { return m_h; }

        point_t grid_t::point(index_t index) const
        {
            assert(index < m_npts);

            point_t point;

            for(auto i = dim; i > 0; --i)
            {
                point[i - 1] = index % m_size[i - 1];
                index /= m_size[i - 1];
            }

            return point;
        }

        index_t grid_t::index(point_t const& point) const
        {
            assert([&] {
                for(auto i = 0; i < dim; ++i)
                    if(point[i] >= m_size[i])
                        return false;
                return true;
            }());

            index_t offset = 0;

            for(auto i = 0; i < dim - 1; ++i)
            {
                offset = m_size[i + 1] * (point[i] + offset);
            }

            offset += point[dim - 1];

            return offset;
        }

        inline bool backwards(index_t dir, index_t dimension)
        {
            assert(dir < n_sweeps);
            assert(dimension < dim);

            return static_cast<bool>(dir & (1 << dimension));
        }

        index_t grid_t::next(index_t idx, index_t dir) const
        {
            assert(idx <= m_npts);
            assert(dir < n_sweeps);

            if(idx == m_npts)
            {
                point_t p;

                for(auto i = 0; i < dim; ++i)
                {
                    p[i] = backwards(dir, i) ? m_size[i] - 2 : 1;
                }

                idx = index(p);
            }
            else
            {
                point_t p = point(idx);

                p[0] += backwards(dir, 0) ? -1 : +1;

                for(auto i = 0; i < dim - 1; ++i)
                {
                    if(p[i] != 0 && p[i] != m_size[i] - 1)
                    {
                        break;
                    }

                    p[i] = backwards(dir, i) ? m_size[i] - 2 : 1;
                    p[i + 1] += backwards(dir, i + 1) ? -1 : +1;
                }

                if(p[dim - 1] == 0 || p[dim - 1] == m_size[dim - 1] - 1)
                {
                    idx = m_npts;
                }
                else
                {
                    idx = index(p);
                }
            }

            return idx;
        }

        inline index_t rotate(index_t boundary_dim, index_t dimension)
        {
            assert(boundary_dim < dim);
            assert(dimension < dim - 1);

            auto idx = boundary_dim + dimension + 1;

            return idx >= dim ? idx - dim : idx;
        }

        index_t grid_t::next_in_boundary(index_t idx, index_t boundary) const
        {
            assert(idx <= m_npts);
            assert(boundary < 2 * dim);

            if(idx == m_npts)
            {
                point_t p;

                for(auto i = 0; i < dim; ++i)
                {
                    p[i] = 0;
                }

                if(boundary >= dim)
                {
                    p[boundary - dim] = m_size[boundary - dim] - 1;
                }

                idx = index(p);
            }
            else
            {
                point_t p = point(idx);

                index_t boundary_dim =
                    boundary >= dim ? boundary - dim : boundary;

                // Increment from the dimension after the boundary
                // dimension to the dimension right before.
                p[rotate(boundary_dim, 0)] += 1;

                for(auto i = 0; i < dim - 2; ++i)
                {
                    if(p[rotate(boundary_dim, i)] ==
                       m_size[rotate(boundary_dim, i)])
                    {
                        p[rotate(boundary_dim, i)] = 0;
                        p[rotate(boundary_dim, i + 1)] += 1;
                    }
                }

                if(p[rotate(boundary_dim, dim - 2)] ==
                   m_size[rotate(boundary_dim, dim - 2)])
                {
                    idx = m_npts;
                }
                else
                {
                    idx = index(p);
                }
            }

            return idx;
        }

    }    // namespace grid
}    // namespace fsm
