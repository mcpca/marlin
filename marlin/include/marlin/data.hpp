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

#include <cassert>
#include <vector>

namespace marlin
{
    //! @brief Functions and types related to \c data_t.
    namespace data
    {
        //! @brief Range of scalar values.
        //
        //! Holds a heap-allocated array and its size.
        template<typename Index, typename Scalar>
        class data_t
        {
          public:
            using index_t = Index;
            using scalar_t = Scalar;

            //! @brief Construct from preexisting data.
            //
            //! @param values data array.
            explicit data_t(std::vector<scalar_t>&& values) noexcept
                : m_values(std::move(values)), m_memsize(m_values.size())
            {}

            //! @brief Allocate new array.
            //
            //! @param memsize number of datapoints.
            //! @param fill value the new array will be filled with.
            data_t(index_t memsize, scalar_t fill = scalar_t{ 0 }) noexcept
                : m_values(std::vector<scalar_t>(memsize)), m_memsize(memsize)
            {
                std::fill(std::begin(m_values), std::end(m_values), fill);
            }

            // Disable copying.
            data_t(data_t const&) = delete;
            data_t& operator=(data_t const&) = delete;

            //! @brief Compiler-generated move constructor.
            data_t(data_t&&) noexcept = default;
            //! @brief Compiler-generated move assignment.
            data_t& operator=(data_t&&) noexcept = default;

            //! @brief Compiler-generated destructor.
            ~data_t() = default;

            //! @brief Read value at an index.
            //
            //! Bounds checking is only performed in debug mode.
            //
            //! @param index index of the value to be read.
            scalar_t at(index_t index) const noexcept
            {
                assert(index < m_memsize);
                return m_values[index];
            }

            //! @brief Get reference to value at an index.
            //
            //! Bounds checking is only performed in debug mode.
            //
            //! @param index index of the value to be read or written to.
            scalar_t& at(index_t index) noexcept
            {
                assert(index < m_memsize);
                return m_values[index];
            }

            //! @brief Get raw pointer to first datapoint.
            scalar_t* get_values() const noexcept { return m_values.data(); }

            //! @brief Move data out of the class.
            std::vector<scalar_t>&& steal() noexcept
            {
                return std::move(m_values);
            }

            //! @brief Size of the underlying array.
            index_t size() const noexcept { return m_memsize; }

          private:
            // Data.
            std::vector<scalar_t> m_values;
            // Number of data points.
            index_t m_memsize;
        };
    }    // namespace data
}    // namespace marlin
