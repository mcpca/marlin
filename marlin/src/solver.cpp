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

#include <algorithm>
#include <numeric>

#include "marlin/solver.hpp"

#include "io.hpp"
#include "levels.hpp"

namespace marlin
{
    namespace solver
    {
        solver_t make_solver(
            std::string const& filename,
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            params_t const& params)
        {
            auto data = io::read(filename, "cost_function");

            if(io::dset_exists(filename, "value_function"))
            {
                throw std::runtime_error(
                    "A dataset named \'value_function\' already exists in the "
                    "given file.");
            }

            return solver_t(filename,
                            std::move(data.data),
                            grid_t(vertices, data.size),
                            params);
        }

        solver_t::solver_t(solver_t&&) noexcept = default;
        solver_t& solver_t::operator=(solver_t&&) noexcept = default;
        solver_t::~solver_t() = default;

        solver_t::solver_t(
            std::string const& filename,
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            params_t const& params)
            : solver_t(make_solver(filename, vertices, params))
        {}

        solver_t::solver_t(std::string filename,
                           data_t cost,
                           grid_t grid,
                           params_t const& params) noexcept
            : m_filename(std::move(filename)),
              m_grid(std::move(grid)),
              m_soln(m_grid.npts(), params.maxval),
              m_cost(std::move(cost)),
              m_tolerance(params.tolerance)
        {
            compute_levels();
            initialize();
        }

        void solver_t::compute_levels() noexcept
        {
            for(index_t i = 0; i < m_grid.n_levels(); ++i)
            {
                level::level_t<dim> level(i, m_grid.size());

                std::vector<point_t> points;

                do
                {
                    point_t point;
                    level.get(point.data());

                    if(!m_grid.is_boundary(point))
                    {
                        points.push_back(point);
                    }
                } while(level.next());

                if(!points.empty())
                {
                    m_levels.emplace_back(std::move(points));
                }
            }
        }

        void solver_t::initialize() noexcept
        {
            MARLIN_DEBUG(auto ntargetpts = 0;)

            for(index_t i = 0; i < m_cost.size(); ++i)
            {
                if(m_cost.at(i) < scalar_t{ 0.0 })
                {
                    MARLIN_DEBUG(++ntargetpts;)

                    m_soln.at(i) = -(m_cost.at(i) + scalar_t{ 1.0 });
                }
            }

            MARLIN_DEBUG(std::cerr << "Found " << ntargetpts
                                   << " target grid points." << '\n';)
        }

        void solver_t::write() const
        {
            std::cout << "Writing to " << m_filename << "..." << std::endl;
            io::write(m_filename, "value_function", m_soln, m_grid.size());
        }

        static inline scalar_t update_boundary(index_t index,
                                               index_t boundary,
                                               data_t const& soln,
                                               grid_t const& grid) noexcept
        {
            assert(boundary < n_boundaries);
            assert(index < grid.npts());

            auto neighbor = grid.point(index);
            auto const boundary_dim =
                boundary >= dim ? boundary - dim : boundary;

            // Approximate based on the two points in a line orthogonal
            // to the boundary closest to the current point.
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto const outer = soln.at(grid.index(neighbor));
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto const inner = soln.at(grid.index(neighbor));

            return std::min(std::max(2 * outer - inner, inner), soln.at(index));
        }

        scalar_t solver_t::boundary_sweep(index_t boundary) noexcept
        {
            assert(boundary < n_boundaries);

            auto const size = m_grid.npts();

            auto diff = scalar_t{ 0 };

            for(auto i = m_grid.next_in_boundary(size, boundary); i != size;
                i = m_grid.next_in_boundary(i, boundary))
            {
                if(m_cost.at(i) > scalar_t{ 0.0 })
                {
                    auto old = m_soln.at(i);
                    m_soln.at(i) = std::min(
                        update_boundary(i, boundary, m_soln, m_grid), old);

                    diff = std::max(diff, old - m_soln.at(i));
                }
            }

            return diff;
        }

        scalar_t solver_t::boundary() noexcept
        {
            scalar_t diff = 0;

            for(index_t bdry = 0; bdry < dim; ++bdry)
            {
                scalar_t diffs[2];

#pragma omp parallel shared(diffs)
                {
#pragma omp sections
                    {
#pragma omp section
                        diffs[0] = boundary_sweep(bdry);
#pragma omp section
                        diffs[1] = boundary_sweep(bdry + dim);
                    }
                }

                diff = std::max(diff, std::max(diffs[0], diffs[1]));
            }

            return diff;
        }
    }    // namespace solver
}    // namespace marlin
