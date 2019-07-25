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

#include <cassert>
#include "fsm/defs.hpp"

#ifndef NDEBUG
#    include <numeric>
#endif

namespace fsm
{
    //! @brief Functions and types related to \c level_t.
    namespace level
    {
        //! @brief Generates the 'level' sets.
        //
        //! A 'level' is a set of gridpoints which can be updated in parallel
        //! (see Detrixhe et al. 2013).
        template<index_t N>
        class level_t
        {
            static_assert(N <= fsm::dim, "");
            friend class level_t<N + 1>;

          public:
            //! @brief Constructor.
            //
            //! @param sum sum of the coordinate indices of the the gridpoints
            //! in the set
            //! @param limits grid size along each dimension.
            level_t(index_t sum, point_t const& limits);

            //! @brief Get the indices of the current point.
            //
            //! @param range pointer to the first element of a range with N
            //! elements
            void get(index_t* range) const;

            //! @brief Generate next point.
            bool next();

          private:
            //! @brief Set the total sum to a new value.
            bool reset(index_t sum);
            index_t sum();

            //! @brief Upper limit for this element.
            index_t m_limit;
            //! @brief Current value of this element.
            index_t m_value;
            //! @brief Next element.
            level_t<N - 1> m_sublevel;
        };

        //! Class specialization for one dimension.
        template<>
        class level_t<1>
        {
            friend class level_t<2>;

          public:
            level_t(index_t sum, point_t const& limits);

            void get(index_t* range) const;
            bool next() const;

          private:
            bool reset(index_t sum);
            index_t sum();

            index_t m_value;
            index_t m_limit;
        };

        template<index_t N>
        level_t<N>::level_t(index_t sum, point_t const& limits)
            : m_limit(limits[dim - N]),
              m_value(std::min(sum, m_limit - 1)),
              m_sublevel(sum > m_value ? sum - m_value : 0, limits)
        {
            assert(std::accumulate(
                       std::begin(limits) + dim - N, std::end(limits), 0) -
                       N >=
                   sum);
        }

        level_t<1>::level_t(index_t sum, point_t const& limits)
            : m_value(sum), m_limit(limits.back())
        {
            assert(m_value < m_limit);
        }

        template<index_t N>
        void level_t<N>::get(index_t* range) const
        {
            assert(range != nullptr);

            *range = m_value;
            m_sublevel.get(range + 1);
        }

        void level_t<1>::get(index_t* range) const
        {
            assert(range != nullptr);
            *range = m_value;
        }

        template<index_t N>
        bool level_t<N>::next()
        {
            if(m_sublevel.next())
                return true;

            if(m_value == 0)
                return false;

            m_value--;
            return m_sublevel.reset(m_sublevel.sum() + 1);
        }

        bool level_t<1>::next() const { return false; }

        template<index_t N>
        bool level_t<N>::reset(index_t sum)
        {
            m_value = std::min(sum, m_limit - 1);
            return m_sublevel.reset(sum > m_value ? sum - m_value : 0);
        }

        bool level_t<1>::reset(index_t sum)
        {
            m_value = sum;
            return m_value < m_limit;
        }

        template<index_t N>
        index_t level_t<N>::sum()
        {
            return m_value + m_sublevel.sum();
        }

        index_t level_t<1>::sum() { return m_value; }

    }    // namespace level
}    // namespace fsm
