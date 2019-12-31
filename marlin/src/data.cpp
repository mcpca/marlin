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

#include "data.hpp"

namespace marlin
{
    namespace data
    {
        data_t::data_t(index_t memsize,
                       std::unique_ptr<scalar_t[]> values) noexcept
            : m_memsize(memsize), m_values(std::move(values))
        {}

        data_t::data_t(index_t memsize, scalar_t fill) noexcept
            : m_memsize(memsize),
              m_values(std::make_unique<scalar_t[]>(m_memsize))
        {
            std::fill_n(this->begin(), m_memsize, fill);
        }

        scalar_t* data_t::get_values() const noexcept { return m_values.get(); }

        index_t data_t::size() const noexcept { return m_memsize; }

        typename data_t::iterator_t data_t::begin() noexcept
        {
            return &m_values[0];
        }

        typename data_t::const_iterator_t data_t::begin() const noexcept
        {
            return &m_values[0];
        }

        typename data_t::iterator_t data_t::end() noexcept
        {
            return &m_values[m_memsize - 1] + 1;
        }

        typename data_t::const_iterator_t data_t::end() const noexcept
        {
            return &m_values[m_memsize - 1] + 1;
        }
    }    // namespace data
}    // namespace marlin
