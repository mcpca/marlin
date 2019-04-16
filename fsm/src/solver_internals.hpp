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

#include "fsm/defs.hpp"

namespace grid
{
    struct grid_t;
}
namespace data
{
    class data_t;
}

namespace fsm
{
    namespace solver
    {
        namespace detail
        {
            void sweep(
                index_t dir,
                data::data_t* soln,
                data::data_t const* cost,
                grid::grid_t const* grid,
                hamiltonian_t const& hamiltonian,
                std::function<vector_t(input_t const&)> const& viscosity);

            void enforce_boundary(index_t boundary,
                                  data::data_t* soln,
                                  data::data_t const* cost,
                                  grid::grid_t const* grid);

        }    // namespace detail
    }        // namespace solver
}    // namespace fsm
