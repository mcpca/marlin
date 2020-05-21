#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "hdf5.hpp"
#include "marlin/defs.hpp"

namespace h5io
{
    // Safe conversion from type to H5T.
    template<typename T>
    constexpr H5::PredType get_type();

    template<>
    H5::PredType get_type<float>()
    {
        return H5::PredType::NATIVE_FLOAT;
    }

    template<>
    H5::PredType get_type<double>()
    {
        return H5::PredType::NATIVE_DOUBLE;
    }

    template<>
    H5::PredType get_type<long double>()
    {
        return H5::PredType::NATIVE_LDOUBLE;
    }

    dataset_t read(std::string const& filename, std::string const& dsetname)
    {
        auto f = H5::H5File(filename, H5F_ACC_RDONLY);

        auto dataset = f.openDataSet(dsetname);
        auto dataspace = dataset.getSpace();

        auto const n_dims = dataspace.getSimpleExtentNdims();
        std::vector<hsize_t> dimensions(n_dims);
        dataspace.getSimpleExtentDims(dimensions.data());

        if(n_dims != marlin::dim)
        {
            throw std::runtime_error(
                "Dataset has an invalid number of dimensions.");
        }

        std::string msg("Found dataset \'" + dsetname + "\' with size (");

        dataset_t rv;

        for(auto i = 0ul; i < static_cast<unsigned long>(n_dims); ++i)
        {
            rv.size[i] = dimensions[i];

            if(i != dimensions.size() - 1)
            {
                msg += std::to_string(dimensions[i]) + ", ";
            }
        }

        auto const n_points = dataspace.getSimpleExtentNpoints();

        msg += std::to_string(dimensions.back()) + ") (total " +
               std::to_string(n_points) + ")";

        std::cerr << msg << '\n';

        rv.data = std::vector<marlin::scalar_t>(n_points);
        dataset.read(rv.data.data(), get_type<marlin::scalar_t>());

        std::cerr << "Successfully read dataset \'" + dsetname + "\'." << '\n';

        return rv;
    }

    static bool dset_exists(H5::H5File const& f,
                            std::string const& dsetname)
    {
        if (!f.nameExists(dsetname))
            return false;
        try {
            f.openDataSet(dsetname);
            return true;
        } catch (...) {
            return false;
        }
    }

    void write(std::string const& filename,
               std::string const& dsetname,
               std::vector<marlin::scalar_t> const& data,
               marlin::point_t const& size)
    {
        auto f = H5::H5File(filename, H5F_ACC_RDWR);

        if(dset_exists(f, dsetname))
        {
            throw std::runtime_error("A dataset named \'" + dsetname +
                                     "\' already exists!");
        }

        std::vector<hsize_t> dims(size.size());
        std::copy(std::begin(size), std::end(size), std::begin(dims));
        auto const space = H5::DataSpace(dims.size(), dims.data());
        auto const type = get_type<marlin::scalar_t>();

        auto dataset = f.createDataSet(dsetname, type, space);
        dataset.write(data.data(), type, space);
    }
}    // namespace h5io
