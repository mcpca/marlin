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

#include "data.hpp"
#include "defs.hpp"

namespace fsm
{
    //! Data input/output
    namespace io
    {
        struct solver_data_t
        {
            data::data_t data;
            point_t size;
        };

        //! Read a dataset from a hdf5 file.
        //! @tparam dim number of spatial dimensions
        //! @tparam scalar_t element type of the container which will hold the
        //! data
        //! @param[in] filename path to the hdf5 file containing the data
        //! @param[in] dsetname name of the hdf5 dataset containing the data
        //! @param[in, out] size records the size along each dimension of the
        //! data
        //! @return gds::Function holding the data
        solver_data_t read(std::string const& filename,
                           std::string const& dsetname);

        //! Write a dataset to a hdf5 file.
        //! @tparam dim number of spatial dimensions
        //! @tparam scalar_t element type of the container holding the data
        //! @param[in] filename path to the hdf5 file to which the data will be
        //! written
        //! @param[in] dsetname name of the hdf5 dataset to which the data will
        //! be written
        //! @param[in] data gds::Function containing the data
        void write(std::string const& filename,
                   std::string const& dsetname,
                   data::data_t const& data,
                   point_t const& size);

        bool dset_exists(std::string const& filename,
                         std::string const& dsetname);

    }    // namespace io
}    // namespace fsm

