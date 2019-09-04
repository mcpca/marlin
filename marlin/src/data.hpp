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
#include <memory>

#include "marlin/defs.hpp"

namespace marlin
{
    //! @brief Functions and types related to \c data_t.
    namespace data
    {
        //! @brief Range of scalar values.
        //
        //! Holds a heap-allocated array and its size.
        class data_t
        {
          public:
            //! @brief Construct from preexisting data.
            //
            //! @param memsize number of datapoints.
            //! @param values unique_ptr holding the data.
            data_t(index_t memsize, std::unique_ptr<scalar_t[]> values);

            //! @brief Allocate new array.
            //
            //! @param memsize number of datapoints.
            //! @param fill value the new array will be filled with.
            explicit data_t(index_t memsize, scalar_t fill = scalar_t{ 0 });

            // Disable copying.
            data_t(data_t const&) = delete;
            data_t& operator=(data_t const&) = delete;

            //! @brief Compiler-generated move constructor.
            data_t(data_t&&) noexcept = default;
            //! @brief Compiler-generated move assignment.
            data_t& operator=(data_t&&) noexcept = default;

            //! @brief Compiler-generated destructor.
            ~data_t() = default;

            //! Iterator type.
            using iterator_t = scalar_t*;
            //! Constant iterator type.
            using const_iterator_t = scalar_t const*;

            //! @brief Read value at an index.
            //
            //! Bounds checking is only performed in debug mode.
            //
            //! @param index index of the value to be read.
            scalar_t at(index_t index) const
            {
                assert(index < m_memsize);
                return m_values[index];
            }

            //! @brief Get reference to value at an index.
            //
            //! Bounds checking is only performed in debug mode.
            //
            //! @param index index of the value to be read or written to.
            scalar_t& at(index_t index)
            {
                assert(index < m_memsize);
                return m_values[index];
            }

            //! @brief Get raw pointer to first datapoint.
            scalar_t* get_values() const;

            //! @brief Size of the underlying array.
            index_t size() const;

            //! @brief Iterator to the beggining of the data.
            iterator_t begin();
            //! @brief Constant iterator to the beggining of the data.
            const_iterator_t begin() const;

            //! @brief Iterator to the end of the data.
            iterator_t end();
            //! @brief Constant iterator to the end of the data.
            const_iterator_t end() const;

          private:
            // Size of the underlying array.
            index_t m_memsize;
            // Pointer to the beggining of the data.
            std::unique_ptr<scalar_t[]> m_values;
        };
    }    // namespace data
}    // namespace marlin
