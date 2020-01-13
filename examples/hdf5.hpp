#pragma once

#include <array>
#include <string>
#include <vector>

#include "marlin/defs.hpp"

namespace h5io
{
    struct dataset_t
    {
        std::vector<marlin::scalar_t> data;
        std::array<marlin::index_t, marlin::dim> size;
    };

    dataset_t read(std::string const& filename, std::string const& dsetname);

    void write(std::string const& filename,
               std::string const& dsetname,
               std::vector<marlin::scalar_t> const& data,
               marlin::point_t const& size);

}    // namespace h5io
