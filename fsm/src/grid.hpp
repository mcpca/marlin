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

#pragma once

#include <array>

#include "fsm/defs.hpp"

namespace fsm
{
    namespace grid
    {
        static_assert(dim > 1, "Number of dimensions must be greater than one");

        struct grid_t
        {
          public:
            grid_t(
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                point_t const& size);

            grid_t(grid_t const&) = default;
            grid_t& operator=(grid_t const&) = default;

            grid_t(grid_t&&) noexcept = default;
            grid_t& operator=(grid_t&&) noexcept = default;

            ~grid_t() = default;

            index_t size(index_t i) const;
            point_t const& size() const;
            index_t npts() const;
            scalar_t h(index_t i) const;
            vector_t const& h() const;

            //! Convert from a row major offset to a list of indices.
            point_t point(index_t index) const;

            //! Convert from a list of indices to a row major offset.
            index_t index(point_t const& point) const;

            //! Get next point in the sweep.
            index_t next(index_t idx, index_t dir) const;

            //! Get next boundary point.
            index_t next_in_boundary(index_t idx, index_t boundary) const;

          private:
            //! Number of points in each dimension.
            point_t m_size;
            //! Total number of gridpoints.
            index_t m_npts;
            //! Diameters of the nodes.
            vector_t m_h;
        };

        class interior_visitor_t
        {
          public:
            struct gridpoint_t
            {
                gridpoint_t(grid_t const* grid, index_t dir);
                gridpoint_t(index_t dir, index_t index);

                gridpoint_t& operator++();
                index_t operator*() const;

                index_t m_index;
                grid_t const* m_grid;
                index_t m_dir;
            };

            interior_visitor_t(grid_t const& grid, index_t dir);

            gridpoint_t const& begin();
            gridpoint_t end();

          private:
            gridpoint_t m_gridpoint;
        };

        bool operator==(interior_visitor_t::gridpoint_t const& a,
                        interior_visitor_t::gridpoint_t const& b);
        bool operator!=(interior_visitor_t::gridpoint_t const& a,
                        interior_visitor_t::gridpoint_t const& b);

        class boundary_visitor_t
        {
          public:
            struct gridpoint_t
            {
                gridpoint_t(grid_t const* grid, index_t boundary);
                gridpoint_t(index_t boundary, index_t index);

                gridpoint_t& operator++();
                index_t operator*() const;

                index_t m_index;
                grid_t const* m_grid;
                index_t m_boundary;
            };

            boundary_visitor_t(grid_t const& grid, index_t boundary);

            gridpoint_t const& begin();
            gridpoint_t end();

          private:
            gridpoint_t m_gridpoint;
        };

        bool operator==(interior_visitor_t::gridpoint_t const& a,
                        interior_visitor_t::gridpoint_t const& b);
        bool operator!=(boundary_visitor_t::gridpoint_t const& a,
                        boundary_visitor_t::gridpoint_t const& b);
    }    // namespace grid
}    // namespace fsm
