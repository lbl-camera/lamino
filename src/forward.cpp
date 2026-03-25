#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include <toml++/toml.hpp>

#include "array_ops.h"
#include "tomocam.h"

constexpr double PADDING =
    1.41421356237; // sqrt(2) to avoid cropping corners of the sample

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " config.toml\n";
        return 1;
    }
    // get input data
    toml::table table;
    try {
        table = toml::parse_file(argv[1]);
    } catch (const toml::parse_error &err) {
        std::cerr << std::format("Parsing failed:\n{}\n", err.description());
        return 1;
    }

    // get parameters from config file
    auto paths = table["input"];
    auto basedir = paths["data_path"].value<std::string>().value();
    auto comp1 = paths["component_1"].value<std::string>().value();
    auto comp2 = paths["component_2"].value<std::string>().value();
    auto comp3 = paths["component_3"].value<std::string>().value();

    auto angles_file = table["angles"]["filename"].value<std::string>().value();

    auto gamma = table["orientation"]["gamma"].value<float>().value();
    gamma = gamma * M_PI / 180.0f; // convert to radians
    auto output_file = table["output"]["filename"].value<std::string>().value();
    auto output_basedir = table["output"]["basedir"].value<std::string>().value();
    std::vector<int> output_dims;
    auto dim_arr = table["output"]["dims"].as_array();
    if (dim_arr && dim_arr->size() == 2) {
        for (const auto &elem : *dim_arr) {
            if (auto val = elem.value<int>()) {
                output_dims.push_back(*val);
            } else {
                std::cerr << "Invalid output dimension value in config file\n";
                return 1;
            }
        }
    } else {
        std::cerr << "Output dimensions must be an array of two integers\n";
        return 1;
    }

    // sanity checks
    std::filesystem::path out_basedir(output_basedir);
    // check if data path exists
    if (!std::filesystem::exists(out_basedir)) {
        std::cerr << "Data path does not exist: " << output_basedir << "\n";
        return 1;
    }
    // set output path
    auto output = (out_basedir / output_file).string();

    // check if components exist
    if (!std::filesystem::exists(std::filesystem::path(basedir) / comp1) ||
        !std::filesystem::exists(std::filesystem::path(basedir) / comp2) ||
        !std::filesystem::exists(std::filesystem::path(basedir) / comp3)) {
        std::cerr << "One or more components do not exist in data path: " << basedir
                  << "\n";
        return 1;
    }

    // read angles from text file
    std::ifstream angles_stream(angles_file);
    if (!angles_stream.is_open()) {
        std::cerr << "Could not open angles file: " << angles_file << "\n";
        return 1;
    }
    std::vector<float> angles;
    float angle;
    while (angles_stream >> angle) { angles.push_back(angle); }
    if (angles.empty()) {
        std::cerr << "No angles found in file: " << angles_file << "\n";
        return 1;
    }

    auto minangle = *std::min_element(angles.begin(), angles.end());
    auto maxangle = *std::max_element(angles.begin(), angles.end());
    if (std::abs(maxangle) > 2 * M_PI) {
        for (auto &a : angles) { a = a * M_PI / 180.0f; }
    }

    // print parameters
    std::cerr << "----------------------------------------\n";
    std::cerr << "Data path: " << basedir << "\n";
    std::cerr << "Component 1: " << comp1 << "\n";
    std::cerr << "Component 2: " << comp2 << "\n";
    std::cerr << "Component 3: " << comp3 << "\n";
    std::cerr << "Angles: [" << minangle << ", " << maxangle << "] with "
              << angles.size() << " steps\n";
    std::cerr << "Gamma: " << gamma << "\n";
    std::cerr << "Output: " << output << "\n";
    std::cerr << "----------------------------------------\n";

    // load data
    auto base_path = std::filesystem::path(basedir);
    std::array<std::string, 3> components = {comp1, comp2, comp3};
    std::array<tomocam::Array<float>, 3> m_data;
    tomocam::Timer t0;
    t0.start();
    for (int i = 0; i < 3; ++i) {
        auto filename = (base_path / components[i]).string();
        m_data[i] = tomocam::tiff::read(filename);
    }

    t0.stop();
    std::cerr << "Time to read data: " << t0.seconds() << "(s)\n";
    std::cerr << "Data dimensions: [" << m_data[0].nslices() << ", "
              << m_data[0].nrows() << ", " << m_data[0].ncols() << "]\n";

    // transpose data [n1, n2, n3] -> [n3, n1, n2]
    t0.start();
    for (int i = 0; i < 3; ++i) {
        m_data[i] = tomocam::array::transpose(m_data[i], {2, 0, 1});
    }
    t0.stop();
    std::cerr << "Time to transpose data: " << t0.seconds() << "(s)\n";

    // pad the sample
    t0.start();
    for (int i = 0; i < 3; ++i) {
        m_data[i] = tomocam::pad3d<float>(m_data[i], PADDING - 1,
                                          tomocam::PadType::SYMMETRIC);
    }
    t0.stop();
    std::cerr << "Time to pad data: " << t0.seconds() << "(s)\n";
    std::cerr << "Padded data dimensions: [" << m_data[0].nslices() << ", "
              << m_data[0].nrows() << ", " << m_data[0].ncols() << "]\n";

    auto padded = [&](int dim) {
        int p = (dim * (PADDING - 1)) / 2;
        return static_cast<size_t>(dim + 2 * p);
    };

    // create a polar grid
    t0.start();
    // oversample polar-grid
    size_t nrows = static_cast<size_t>(padded(output_dims[0]));
    size_t ncols = static_cast<size_t>(padded(output_dims[1]));
    tomocam::PolarGrid<float> grid(angles, nrows, ncols, gamma);
    t0.stop();
    std::cerr << "Time to build a polar grid: " << t0.seconds() << "(s)\n";
    std::cerr << "Polar grid dimensions: [" << grid.dims().n1 << ", "
              << grid.dims().n2 << ", " << grid.dims().n3 << "]\n";

    // do the forward projection
    t0.start();
    auto proj = tomocam::forward(m_data, grid, gamma);
    t0.stop();
    std::cerr << "time to do forward projection: " << t0.seconds() << "(s)\n";

    // crop the projection to original size
    tomocam::dims_t crop_dims = {angles.size(), static_cast<size_t>(output_dims[0]),
                                 static_cast<size_t>(output_dims[1])};
    t0.start();
    proj = tomocam::crop2d<float>(proj, crop_dims, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << "Time to crop data: " << t0.seconds() << "(s)\n";

    //  save data to tiff-stack
    tomocam::tiff::write(output, proj);

    return 0;
}
