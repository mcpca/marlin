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

#include "data.hpp"
#include "marlin/defs.hpp"

namespace marlin
{
    //! @brief Functions for data input/output with hdf5 files.
    namespace io
    {
        //! @brief Holds a dataset and its size along each dimension.
        struct solver_data_t
        {
            data::data_t data;
            point_t size;
        };

        //! @brief Read a dataset from a hdf5 file.
        //
        //! @param filename path to the file.
        //! @param dsetname path to the dataset in the file.
        //! @return the data and its dimensions.
        solver_data_t read(std::string const& filename,
                           std::string const& dsetname);

        //! @brief Write data to a hdf5 file.
        //
        //! Creates a new dataset and writes to it.
        //
        //! @param filename path to the file.
        //! @param dsetname path of the dataset to be created.
        //! @param data the data to be written.
        //! @param size the dimensions of the data.
        void write(std::string const& filename,
                   std::string const& dsetname,
                   data::data_t const& data,
                   point_t const& size);

        //! @brief Check whether a dataset exists in a hdf5 file.
        //
        //! @param filename path to the file.
        //! @param dsetname path to the dataset being queried.
        //! @return true if the dataset exists.
        bool dset_exists(std::string const& filename,
                         std::string const& dsetname);

    }    // namespace io
}    // namespace marlin
