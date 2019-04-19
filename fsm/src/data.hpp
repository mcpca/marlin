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
#include <cassert>
#include <memory>

#include "fsm/defs.hpp"

namespace fsm
{
    namespace data
    {
        //! A simple wrapper around a unique_ptr holding a dynamic array.
        //! @tparam scalar_t type of the array elements
        class data_t
        {
          public:
            //! Take over the data from an existing unique_ptr holding an array.
            //! Make sure that the allocated size of the array is at least
            //! memsize.
            //! @param[in] memsize size of the array
            //! @param[in] values unique_ptr to be moved from
            data_t(index_t memsize, std::unique_ptr<scalar_t[]> values);

            //! Creates a new unique_ptr and fills it with a given value.
            //! @param[in] memsize size of the array
            //! @param[in] fill value to fill the array with
            explicit data_t(index_t memsize, scalar_t fill = scalar_t{ 0 });

            // Disable copying
            data_t(data_t const&) = delete;
            data_t& operator=(data_t const&) = delete;

            data_t(data_t&&) noexcept = default;
            data_t& operator=(data_t&&) noexcept = default;

            ~data_t() = default;

            //! Iterator type.
            using iterator_t = scalar_t*;

            //! Constant iterator type.
            using const_iterator_t = scalar_t const*;

            //! Read value at an index.
            //! No bounds checking is performed.
            //! @param[in] index index of the value to be read in the array
            //! @return the value held by the array at the given index
            scalar_t at(index_t index) const
            {
                assert(index < m_memsize);
                return m_values[index];
            }

            //! Get reference to value at a gridpoint.
            //! No bounds checking is performed.
            //! @param[in] index index of the value to be read in the array
            //! @return reference to the corresponding position of the array
            scalar_t& at(index_t index)
            {
                assert(index < m_memsize);
                return m_values[index];
            }

            //! Get pointer to the first position of the array (read only)
            //! @return pointer to the first position of the array
            scalar_t* get_values() const;

            //! @return size of the underlying array
            index_t size() const;

            //! @return iterator to the beggining of the array
            iterator_t begin();

            //! @return iterator to the beggining of the array
            const_iterator_t begin() const;

            //! @return iterator to the element following the last element in
            //! the array
            iterator_t end();

            //! @return iterator to the element following the last element in
            //! the array
            const_iterator_t end() const;

          private:
            //! Number of gridpoints.
            index_t m_memsize;
            //! Array holding the value at each gridpoint.
            std::unique_ptr<scalar_t[]> m_values;
        };
    }    // namespace data
}    // namespace fsm
