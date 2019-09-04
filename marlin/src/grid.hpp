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

#include "marlin/defs.hpp"

namespace marlin
{
    //! @brief Functions and types related to \c grid_t.
    namespace grid
    {
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
                point_t const& size);

            //! @brief Compiler-generated copy constructor.
            grid_t(grid_t const&) = default;
            //! @brief Compiler-generated copy assignment.
            grid_t& operator=(grid_t const&) = default;

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
            index_t size(index_t i) const;
            //! @brief Size of the grid in each dimension.
            point_t const& size() const;
            //! @brief Total number of gridpoints.
            index_t npts() const;

            //! @brief Grid resolution in a dimension.
            //
            //! @param i dimension.
            scalar_t h(index_t i) const;
            //! @brief Grid resolution in each dimension.
            vector_t const& h() const;

            //! @brief Number of levels.
            index_t n_levels() const;

            //! @brief Convert from a row major offset to a list of indices.
            //
            //! @param index row major offset of a point.
            //! @return indices of the corresponding gridpoint.
            point_t point(index_t index) const;

            //! @brief Convert from a list of indices to a row major offset.
            //
            //! @param point indices of a gridpoint.
            //! @return row major index of \p point.
            index_t index(point_t const& point) const;

            //! @brief Rotate the axes to match a given sweep order.
            //
            //! Convert point from the standard ordering to the ordering given
            //! by a specified sweeping direction.
            //
            //! @param point indices of a gridpoint.
            //! @param dir sweeping direction.
            //! @return \p point after the axes rotation.
            point_t rotate_axes(point_t point, int dir) const;

            //! @brief Check if a point is in the computational boundary.
            //
            //! @param point indices of a gridpoint.
            //! @return true if \p point is in the computational boundary.
            bool is_boundary(point_t const& point) const;

            //! @brief Get the next point according to sweeping order \p dir.
            //
            //! By calling this function repeatedly with the same value of
            //! \p dir, one can sweep the interior gridpoints in the direction
            //! given by \p dir.
            //
            //! @param idx row major index of current point.
            //! @param dir sweeping direction.
            //! @return index of the next point.
            index_t next(index_t idx, index_t dir) const;

            //! @brief Get next point in boundary \p boundary.
            //
            //! By calling this function repeatedly with the same value of
            //! \p dir, one can iterate through the grid boundary given by
            //! \p boundary.
            //
            //! @param idx row major index of current point.
            //! @param boundary boundary identifier.
            //! @return index of the next point.
            index_t next_in_boundary(index_t idx, index_t boundary) const;

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
