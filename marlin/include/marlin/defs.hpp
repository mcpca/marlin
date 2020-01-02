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
#include <functional>
#include <type_traits>

#ifndef MARLIN_N_DIMS
#    define MARLIN_N_DIMS 3
#endif

#ifndef MARLIN_N_WORKERS
#    include <thread>
#endif

namespace marlin
{
    //! Number of dimensions.
    constexpr auto dim = MARLIN_N_DIMS;
    static_assert(std::is_integral<decltype(dim)>::value,
                  "Number of dimensions must be a positive integer.");
    static_assert(dim > 1, "The number of dimensions must be at least two.");

    //! Number of sweeping directions.
    constexpr auto n_sweeps = 1 << dim;
    //! Number of boundaries.
    constexpr auto n_boundaries = 2 * dim;

#ifndef MARLIN_N_WORKERS
    //! Number of threads.
    //! If the number of threads is not given, try hardware_concurrency.
    static int const n_workers = [] {
        int const hc = std::thread::hardware_concurrency();
        return hc > 0 ? hc : 2;
    }();
#else
    constexpr auto n_workers = MARLIN_N_WORKERS;
    static_assert(std::is_integral<decltype(n_workers)>::value,
                  "Number of workers must be a positive integer.");
    static_assert(n_workers > 1, "The number of workers must be at least two.");
#endif

    //! Type of the integer indices for indexing into the grid and arrays.
    using index_t = size_t;
    //! Gridpoint representation as a list of indices.
    using point_t = std::array<index_t, dim>;

    //! Floating point type.
    using scalar_t = double;
    //! Representation of the value of a vector field at some gridpoint.
    using vector_t = std::array<scalar_t, dim>;
}    // namespace marlin
