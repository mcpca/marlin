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

#include <iostream>
#include <numeric>

#include "data.hpp"
#include "grid.hpp"
#include "io.hpp"
#include "fsm/solver.hpp"

namespace fsm
{
    namespace solver
    {
        struct solver_t::update_data_internal_t
        {
            vector_t p;
            vector_t avgs;
        };

        solver_t make_solver(
            std::string const& filename,
            hamiltonian_t const& hamiltonian,
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            vector_t const& viscosity,
            params_t const& params)
        {
            auto data = io::read(filename, std::string("cost_function"));

            if(io::dset_exists(filename, "value_function"))
            {
                throw std::runtime_error(
                    "A dataset named \'value_function\' already exists in the "
                    "given file.");
            }

            grid::grid_t g(vertices, data.size);
            return solver_t(filename,
                            hamiltonian,
                            std::move(data.data),
                            g,
                            viscosity,
                            params);
        }

        solver_t::solver_t(solver_t&&) noexcept = default;
        solver_t& solver_t::operator=(solver_t&&) noexcept = default;
        solver_t::~solver_t() = default;

        solver_t::solver_t(
            std::string const& filename,
            hamiltonian_t const& hamiltonian,
            std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
            vector_t const& viscosity,
            params_t const& params)
            : solver_t(make_solver(filename,
                                   hamiltonian,
                                   vertices,
                                   viscosity,
                                   params))
        {}

        solver_t::solver_t(std::string const& filename,
                           hamiltonian_t const& hamiltonian,
                           data::data_t cost,
                           grid::grid_t const& grid,
                           vector_t const& viscosity,
                           params_t const& params)
            : m_filename(filename),
              m_hamiltonian(hamiltonian),
              m_grid(std::make_unique<grid::grid_t>(grid)),
              m_soln(std::make_unique<data::data_t>(m_grid->npts(), params.maxval)),
              m_cost(std::make_unique<data::data_t>(std::move(cost))),
              m_viscosity(viscosity),
              m_c(scalar_t{ 1.0 } /
                  std::inner_product(std::begin(m_viscosity),
                                     std::end(m_viscosity),
                                     std::begin(m_grid->h()),
                                     0.0,
                                     std::plus<>(),
                                     std::divides<>())),
              m_tolerance(params.tolerance)
        {}

        void solver_t::solve()
        {
            std::cout << "Initializing problem instance..." << std::endl;
            initialize();
            std::cout << "Initialized. Solving..." << std::endl;

#ifndef NDEBUG
            auto niter = 0;
            std::cerr << "Iteration " << niter++ << ":\n";
#endif

            while(!iterate())
            {
#ifndef NDEBUG
                std::cerr << "Iteration " << niter++ << ":\n";
#else
                ;
#endif
            }

            std::cout << "Done. Writing to " << m_filename << "..."
                      << std::endl;
            io::write(m_filename, "value_function", *m_soln, m_grid->size());
        }

        void solver_t::initialize()
        {
#ifndef NDEBUG
            auto ntargetpts = 0;
#endif
            for(index_t i = 0; i < m_cost->size(); ++i)
            {
                if(m_cost->at(i) < scalar_t{ 0.0 })
                {
#ifndef NDEBUG
                    ++ntargetpts;
#endif
                    m_soln->at(i) = -(m_cost->at(i) + scalar_t{ 1.0 });
                }
            }
#ifndef NDEBUG
            std::cerr << "Found " << ntargetpts << " target grid points."
                      << '\n';
#endif
        }

        bool solver_t::iterate()
        {
            for(index_t dir = 0; dir < n_sweeps; ++dir)
            {
                auto diff = sweep(dir);

#ifndef NDEBUG
                std::cerr << "Sweep " << dir << ": "
                          << "delta = " << diff << '\n';
#endif

                for(index_t boundary = 0; boundary < n_boundaries; ++boundary)
                {
                    diff = std::max(diff, enforce_boundary(boundary));
                }

#ifndef NDEBUG
                std::cerr << "Sweep " << dir << " (after boundary): "
                          << "delta = " << diff << '\n';
#endif

                // If the difference from the last iteration is less than the
                // specified tolerance value, stop.
                if(diff < m_tolerance)
                {
                    return true;
                }
            }

            return false;
        }

        scalar_t solver_t::sweep(index_t dir)
        {
            auto diff = scalar_t{ 0.0 };
            for(auto i = m_grid->next(m_grid->npts(), dir); i != m_grid->npts();
                i = m_grid->next(i, dir))
            {
                if(m_cost->at(i) > scalar_t{ 0.0 })
                {
                    scalar_t old = m_soln->at(i);
                    m_soln->at(i) = std::min(old, update(i));

                    diff = std::max(diff, old - m_soln->at(i));
                }
            }

            return diff;
        }

        scalar_t solver_t::enforce_boundary(index_t boundary)
        {
            auto diff = scalar_t{ 0.0 };

            for(auto i = m_grid->next_in_boundary(m_grid->npts(), boundary);
                i != m_grid->npts();
                i = m_grid->next_in_boundary(i, boundary))
            {
                if(m_cost->at(i) > scalar_t{ 0.0 })
                {
                    auto old = m_soln->at(i);
                    auto res = update_boundary(i, boundary);
                    m_soln->at(i) = std::min(old, res);
                    diff = std::max(diff, old - m_soln->at(i));
                }
            }

            return diff;
        }

        inline scalar_t solver_t::update(index_t index) const
        {
            auto data = estimate_p(index);

            return m_c * (m_cost->at(index) -
                          m_hamiltonian(m_grid->point(index), data.p) +
                          std::inner_product(std::begin(m_viscosity),
                                             std::end(m_viscosity),
                                             std::begin(data.avgs),
                                             0.0));
        }

        inline scalar_t solver_t::update_boundary(index_t index,
                                                  index_t boundary) const
        {
            auto neighbor = m_grid->point(index);
            auto boundary_dim = boundary >= dim ? boundary - dim : boundary;

            // Approximate based on the two points in a line perpendicular to
            // the boundary closest to the current point.
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto outer = m_soln->at(m_grid->index(neighbor));
            neighbor[boundary_dim] += boundary >= dim ? -1 : +1;
            auto inner = m_soln->at(m_grid->index(neighbor));

            return std::min(std::max(2 * outer - inner, inner),
                            m_soln->at(index));
        }

        inline typename solver_t::update_data_internal_t solver_t::estimate_p(
            index_t index) const
        {
            // Centered finite-difference estimate of the gradient.
            vector_t p;
            // Scaled averages of the neighboring values along each dimension
            // for the dissipation terms.
            vector_t avgs;

            auto point = m_grid->point(index);

            for(index_t i = 0; i < dim; ++i)
            {
                auto neighbor = point;
                neighbor[i] += 1;
                auto right = m_soln->at(m_grid->index(neighbor));
                neighbor[i] -= 2;
                auto left = m_soln->at(m_grid->index(neighbor));

                p[i] = (right - left) / (scalar_t{ 2.0 } * m_grid->h(i));
                avgs[i] = (right + left) / (scalar_t{ 2.0 } * m_grid->h(i));
            }

            return { p, avgs };
        }
    }    // namespace solver
}    // namespace fsm
