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

#include <numeric>

namespace fsm
{
    namespace level
    {
        class level_t
        {
          public:
            level_t(int sum, point_t const& limits);

            point_t const& get() const;
            void next_level();
            bool next();

          private:
            point_t m_range{ 0, 0, 0 };
            point_t m_limits;
            unsigned m_sum;
        };

        level_t::level_t(int sum, point_t const& limits)
            : m_limits(limits), m_sum(sum)
        {
            m_range[0] = m_sum;

            if(m_sum >= m_limits[0])
            {
                m_range[0] = m_limits[0] - 1;
                m_range[1] = m_sum - m_limits[0] + 1;
            }

            if(m_range[1] >= m_limits[1])
            {
                m_range[2] = m_range[1] - m_limits[1] + 1;
                m_range[1] = m_limits[1] - 1;
            }
        }

        bool level_t::next()
        {
            if((m_range[1] == 0) || (m_range[2] == m_limits[2] - 1))
            {
                if((m_range[0] == 0) || m_range[1] == m_limits[1] - 1)
                {
                    return false;
                }

                m_range[0]--;
                m_range[1] = m_sum - m_range[0];

                if(m_range[1] >= m_limits[1])
                {
                    m_range[2] = m_range[1] - m_limits[1] + 1;
                    m_range[1] = m_limits[1] - 1;
                }
                else
                {
                    m_range[2] = 0;
                }

                return true;
            }

            m_range[1]--;
            m_range[2]++;

            return true;
        }

        void level_t::next_level()
        {
            ++m_sum;

            m_range[0] = m_sum;

            if(m_sum >= m_limits[0])
            {
                m_range[0] = m_limits[0] - 1;
                m_range[1] = m_sum - m_limits[0] + 1;
            }

            if(m_range[1] >= m_limits[1])
            {
                m_range[2] = m_range[1] - m_limits[1] + 1;
                m_range[1] = m_limits[1] - 1;
            }
        }

        point_t const& level_t::get() const { return m_range; }
    }    // namespace level
}    // namespace fsm
