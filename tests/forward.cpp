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

#include "tomocam.h"

constexpr double PADDING = 1.4142;

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
    auto output = table["output"]["filename"].value<std::string>().value();
    // sanity checks
    // check if data path exists
    if (!std::filesystem::exists(basedir)) {
        std::cerr << "Data path does not exist: " << basedir << "\n";
        return 1;
    }
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

    // check if output path is valid
    auto output_path = std::filesystem::path(output);
    if (output_path.has_parent_path() &&
        !std::filesystem::exists(output_path.parent_path())) {
        std::cerr << "Output path does not exist: " << output_path.parent_path()
                  << "\n";
        return 1;
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

    // save the original size
    tomocam::dims_t dims = m_data[0].dims();

    // pad the sample
    t0.start();
    for (int i = 0; i < 3; ++i) {
        m_data[i] =
            tomocam::pad3d<float>(m_data[i], PADDING, tomocam::PadType::SYMMETRIC);
    }
    t0.stop();
    std::cerr << "Time to pad data: " << t0.seconds() << "(s)\n";

    // create a polar grid
    t0.start();
    // oversample polar-grid
    auto nrows = m_data[0].nrows();
    auto ncols = m_data[0].ncols();
    tomocam::PolarGrid<float> grid(angles, nrows, ncols, gamma);
    t0.stop();
    std::cerr << "Time to build a polar grid: " << t0.seconds() << "(s)\n";

    // do the forward projection
    t0.start();
    auto proj = tomocam::forward(m_data, grid, gamma);
    t0.stop();
    std::cerr << "time to do forward projection: " << t0.seconds() << "(s)\n";

    // crop the projection to original size
    tomocam::dims_t crop_dims = {angles.size(), dims.n2, dims.n3};
    t0.start();
    proj = tomocam::crop2d<float>(proj, crop_dims, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << "Time to crop data: " << t0.seconds() << "(s)\n";
    // save data to tiff-stack
    tomocam::tiff::write(output, proj);

    return 0;
}
