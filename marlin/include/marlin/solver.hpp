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
#include <memory>
#include <string>
#include <vector>

#include "defs.hpp"

class ThreadPool;

namespace marlin
{
    namespace grid
    {
        struct grid_t;
    }
    namespace data
    {
        class data_t;
    }

    namespace solver
    {
        //! @brief Numerical parameters supplied by the user.
        struct params_t
        {
            //! Maximum value of the solution over the whole domain.
            scalar_t maxval;
            //! Tolerance parameter for the convergence criterion.
            scalar_t tolerance;
        };

        //! @brief Solver API.
        //
        //! The main class, used for defining a problem, solving it and writing
        //! the solution to disk.
        class solver_t
        {
          public:
            //! Contructs a solver_t object given an HDF5 file and the problem
            //! data.
            //! Construction is delegated to the move constructor by calling a
            //! factory function, which makes it easier to read the data in the
            //! HDF5 file before constructing the solver_t object.
            //
            //! @param filename Path to a HDF5 containing the cost function.
            //! @param hamiltonian Callable for evaluating the Hamiltonian.
            //! @param vertices The limits of the grid.
            //! @param viscosity Callable for evaluating the viscosity
            //!                  coefficients.
            //! @param params Numerical parameters.
            solver_t(
                std::string const& filename,
                hamiltonian_t const& hamiltonian,
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                std::function<vector_t(input_t const&)> const& viscosity,
                params_t const& params);

            //! Compiler-generated move contructor.
            solver_t(solver_t&&) noexcept;
            //! Compiler-generated move assignment.
            solver_t& operator=(solver_t&&) noexcept;
            //! Compiler-generated destructor.
            ~solver_t();

            //! Initializes and solves the problem instance, and writes the
            //! solution to disk.
            void solve();

          private:
            //! @brief Factory function.
            //
            //! Reads and processes the necessary data from the hdf5 file before
            //! calling a private constructor.
            friend solver_t make_solver(
                std::string const& filename,
                hamiltonian_t const& hamiltonian,
                std::array<std::pair<scalar_t, scalar_t>, dim> const& vertices,
                std::function<vector_t(input_t const&)> const& viscosity,
                params_t const& params);

            // Constructs a solver_t object after the relevant info has been
            // read from the HDF5 file.
            solver_t(std::string const& filename,
                     hamiltonian_t const& hamiltonian,
                     data::data_t cost,
                     grid::grid_t const& grid,
                     std::function<vector_t(input_t const&)> const& viscosity,
                     params_t const& params);

            // Initialize the solution.
            void initialize();
            // Main loop.
            bool iterate();

            // Sweeps all gridpoints in the direction dir.
            scalar_t sweep(int dir);
            // Updates the boundary points.
            scalar_t boundary();

            // Updates a group of points.
            scalar_t update_points(std::vector<point_t> const* points,
                                   int dir,
                                   int start,
                                   int end);

            // Path to a HDF5 file containing the cost function.
            std::string m_filename;
            // Callable for evaluating the Hamiltonian.
            hamiltonian_t m_hamiltonian;

            // Privately implemented members.
            std::unique_ptr<grid::grid_t> m_grid;    // Grid.
            std::unique_ptr<data::data_t> m_soln;    // Solution.
            std::unique_ptr<data::data_t> m_cost;    // Cost function.

            // Callable for evalutating the viscosity coefficients.
            std::function<vector_t(input_t const&)> m_viscosity;
            // Tolerance for the convergence criterion.
            scalar_t m_tolerance;

            // Thread pool.
            std::unique_ptr<ThreadPool> m_pool;

            // Caches the sets of points which can be updated in parallel in
            // each sweep.
            std::vector<std::vector<point_t>> m_levels;
        };
    }    // namespace solver
}    // namespace marlin
