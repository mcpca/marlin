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

#include <iostream>

#include <h5cpp/hdf5.hpp>
#include <h5cpp/utilities/array_adapter.hpp>
#include "io.hpp"

namespace marlin
{
    namespace io
    {
        solver_data_t read(std::string const& filename,
                           std::string const& dsetname)
        {
            auto f =
                hdf5::file::open(filename, hdf5::file::AccessFlags::READONLY);

            auto dataset = f.root().get_dataset(dsetname);
            auto dataspace = hdf5::dataspace::Simple(dataset.dataspace());
            auto dimensions = dataspace.current_dimensions();

            if(dimensions.size() != dim)
            {
                throw std::runtime_error(
                    "Dataset has an invalid number of dimensions.");
            }

            std::string msg("Found dataset \'" + dsetname + "\' with size (");

            point_t size;

            for(auto i = 0ul; i < dimensions.size(); ++i)
            {
                size[i] = dimensions[i];

                if(i != dimensions.size() - 1)
                {
                    msg += std::to_string(dimensions[i]) + ", ";
                }
            }

            msg += std::to_string(dimensions.back()) + ") (total " +
                   std::to_string(dataspace.size()) + ")";

            std::cout << msg << '\n';

            auto raw_data = std::make_unique<scalar_t[]>(dataspace.size());

            hdf5::ArrayAdapter<scalar_t> adapter(raw_data.get(),
                                                 dataspace.size());
            dataset.read(adapter);

            auto data = data_t(dataspace.size(), std::move(raw_data));

            std::cout << "Successfully read dataset \'" + dsetname + "\'."
                      << '\n';

            return { std::move(data), size };
        }

        bool dset_exists(hdf5::file::File const& f, std::string const& dsetname)
        {
            return f.root().has_dataset(dsetname);
        }

        bool dset_exists(std::string const& filename,
                         std::string const& dsetname)
        {
            return dset_exists(
                hdf5::file::open(filename, hdf5::file::AccessFlags::READONLY),
                dsetname);
        }

        void write(std::string const& filename,
                   std::string const& dsetname,
                   data_t const& data,
                   point_t const& size)
        {
            auto f =
                hdf5::file::open(filename, hdf5::file::AccessFlags::READWRITE);

            if(dset_exists(f, dsetname))
            {
                throw std::runtime_error("A dataset named \'" + dsetname +
                                         "\' already exists!");
            }

            auto* data_values = data.get_values();

            auto dims = hdf5::Dimensions(size.size());
            std::copy(std::begin(size), std::end(size), std::begin(dims));

            hdf5::ArrayAdapter<scalar_t> adapter(data_values, data.size());
            auto dataset =
                hdf5::node::Dataset(f.root(),
                                    dsetname,
                                    hdf5::datatype::create<scalar_t>(),
                                    hdf5::dataspace::Simple(dims));

            dataset.write(adapter);
        }

    }    // namespace io
}    // namespace marlin
