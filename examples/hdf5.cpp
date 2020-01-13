#include <array>
#include <string>
#include <vector>

#include <h5cpp/hdf5.hpp>

#include "hdf5.hpp"
#include "marlin/defs.hpp"

namespace h5io
{
    dataset_t read(std::string const& filename, std::string const& dsetname)
    {
        auto f = hdf5::file::open(filename, hdf5::file::AccessFlags::READONLY);

        auto dataset = f.root().get_dataset(dsetname);
        auto dataspace = hdf5::dataspace::Simple(dataset.dataspace());
        auto dimensions = dataspace.current_dimensions();

        if(dimensions.size() != marlin::dim)
        {
            throw std::runtime_error(
                "Dataset has an invalid number of dimensions.");
        }

        std::string msg("Found dataset \'" + dsetname + "\' with size (");

        dataset_t rv;

        for(auto i = 0ul; i < dimensions.size(); ++i)
        {
            rv.size[i] = dimensions[i];

            if(i != dimensions.size() - 1)
            {
                msg += std::to_string(dimensions[i]) + ", ";
            }
        }

        msg += std::to_string(dimensions.back()) + ") (total " +
               std::to_string(dataspace.size()) + ")";

        std::cout << msg << '\n';

        rv.data = std::vector<marlin::scalar_t>(dataspace.size());
        dataset.read(rv.data);

        std::cout << "Successfully read dataset \'" + dsetname + "\'." << '\n';

        return rv;
    }

    static bool dset_exists(hdf5::file::File const& f,
                            std::string const& dsetname)
    {
        return f.root().has_dataset(dsetname);
    }

    void write(std::string const& filename,
               std::string const& dsetname,
               std::vector<marlin::scalar_t> const& data,
               marlin::point_t const& size)
    {
        auto f = hdf5::file::open(filename, hdf5::file::AccessFlags::READWRITE);

        if(dset_exists(f, dsetname))
        {
            throw std::runtime_error("A dataset named \'" + dsetname +
                                     "\' already exists!");
        }

        auto dims = hdf5::Dimensions(size.size());
        std::copy(std::begin(size), std::end(size), std::begin(dims));

        auto dataset =
            hdf5::node::Dataset(f.root(),
                                dsetname,
                                hdf5::datatype::create<marlin::scalar_t>(),
                                hdf5::dataspace::Simple(dims));

        dataset.write(data);
    }
}    // namespace h5io

