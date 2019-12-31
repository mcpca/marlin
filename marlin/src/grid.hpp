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
#include <cassert>
#include <numeric>

#include "marlin/defs.hpp"

namespace marlin
{
    //! @brief Functions and types related to \c grid_t.
    namespace grid
    {
        inline bool backwards(index_t dir, index_t dimension) noexcept
        {
            assert(dir < n_sweeps);
            assert(dimension < dim);

            return static_cast<bool>(dir & (1 << dimension));
        }

        inline index_t rotate(index_t boundary_dim, index_t dimension) noexcept
        {
            assert(boundary_dim < dim);
            assert(dimension < dim - 1);

            auto idx = boundary_dim + dimension + 1;

            return idx >= dim ? idx - dim : idx;
        }

        //! @brief Implements all grid-related operations.
        //
        //! Holds the grid parameters and implements the conversion
        //! between row major offsets and lists of indices.
        struct grid_t
        {
          public:
            //! @brief Constructor.
            //
            //! Construct a grid_t object given the limits of the grid and the
            //! number of points in each dimension.
            grid_t(
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                point_t const& size) noexcept
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
                  }()),
                  m_nlevels(std::accumulate(std::begin(m_size),
                                            std::end(m_size),
                                            0ul) -
                            dim + 1)
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

            //! @brief Compiler-generated copy constructor.
            grid_t(grid_t const&) noexcept = default;
            //! @brief Compiler-generated copy assignment.
            grid_t& operator=(grid_t const&) noexcept = default;

            //! @brief Compiler-generated move constructor.
            grid_t(grid_t&&) noexcept = default;
            //! @brief Compiler-generated move assignment.
            grid_t& operator=(grid_t&&) noexcept = default;

            //! @brief Compiler-generated destructor.
            ~grid_t() = default;

            //! @brief Size of the grid in a dimension.
            //
            //! @param i dimension.
            //! @return size along dimension \p i.
            index_t size(index_t i) const noexcept
            {
                assert(i < dim);
                return m_size[i];
            }

            //! @brief Size of the grid in each dimension.
            point_t const& size() const noexcept { return m_size; }
            //! @brief Total number of gridpoints.
            index_t npts() const noexcept { return m_npts; }

            //! @brief Grid resolution in a dimension.
            //
            //! @param i dimension.
            scalar_t h(index_t i) const noexcept
            {
                assert(i < dim);
                return m_h[i];
            }
            //! @brief Grid resolution in each dimension.
            vector_t const& h() const noexcept { return m_h; }

            //! @brief Number of levels.
            index_t n_levels() const noexcept { return m_nlevels; }

            //! @brief Convert from a row major offset to a list of indices.
            //
            //! @param index row major offset of a point.
            //! @return indices of the corresponding gridpoint.
            point_t point(index_t index) const noexcept
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

            //! @brief Convert from a list of indices to a row major offset.
            //
            //! @param point indices of a gridpoint.
            //! @return row major index of \p point.
            index_t index(point_t const& point) const noexcept
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

            //! @brief Rotate the axes to match a given sweep order.
            //
            //! Convert point from the standard ordering to the ordering given
            //! by a specified sweeping direction.
            //
            //! @param point indices of a gridpoint.
            //! @param dir sweeping direction.
            //! @return \p point after the axes rotation.
            point_t rotate_axes(point_t point, int dir) const noexcept
            {
                for(int i = 0; i < dim; ++i)
                {
                    if(backwards(dir, i))
                    {
                        point[i] = m_size[i] - point[i] - 1;
                    }
                }

                return point;
            }

            //! @brief Check if a point is in the computational boundary.
            //
            //! @param point indices of a gridpoint.
            //! @return true if \p point is in the computational boundary.
            bool is_boundary(point_t const& point) const noexcept
            {
                for(auto i = 0; i < dim; ++i)
                {
                    if(point[i] == 0 || point[i] == m_size[i] - 1)
                    {
                        return true;
                    }
                }

                return false;
            }

            //! @brief Get the next point according to sweeping order \p dir.
            //
            //! By calling this function repeatedly with the same value of
            //! \p dir, one can sweep the interior gridpoints in the direction
            //! given by \p dir.
            //
            //! @param idx row major index of current point.
            //! @param dir sweeping direction.
            //! @return index of the next point.
            index_t next(index_t idx, index_t dir) const noexcept
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

            //! @brief Get next point in boundary \p boundary.
            //
            //! By calling this function repeatedly with the same value of
            //! \p dir, one can iterate through the grid boundary given by
            //! \p boundary.
            //
            //! @param idx row major index of current point.
            //! @param boundary boundary identifier.
            //! @return index of the next point.
            index_t next_in_boundary(index_t idx, index_t boundary) const
                noexcept
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

          private:
            // Number of points in each dimension.
            point_t m_size;
            // Total number of gridpoints.
            index_t m_npts;
            // Diameters of the nodes.
            vector_t m_h;
            // Number of levels.
            index_t m_nlevels;
        };
    }    // namespace grid
}    // namespace marlin
