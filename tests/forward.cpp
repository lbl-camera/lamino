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

    auto angles = table["angles"];
    auto minangle = angles["min"].value<float>().value();
    auto maxangle = angles["max"].value<float>().value();
    auto nangles = angles["steps"].value<int>().value();

    auto gamma = table["orientation"]["gamma"].value<float>().value();
    gamma = gamma * M_PI / 180.0f; // convert to radians
    auto output = table["output"]["output"].value<std::string>().value();

    // print parameters
    std::cerr << "----------------------------------------\n";
    std::cerr << "Data path: " << basedir << "\n";
    std::cerr << "Component 1: " << comp1 << "\n";
    std::cerr << "Component 2: " << comp2 << "\n";
    std::cerr << "Component 3: " << comp3 << "\n";
    std::cerr << "Angles: [" << minangle << ", " << maxangle << "] with " << nangles
              << " steps\n";
    std::cerr << "Gamma: " << gamma << "\n";
    std::cerr << "Output: " << output << "\n";
    std::cerr << "----------------------------------------\n";

    // load data
    auto base_path = std::filesystem::path(basedir);
    std::array<tomocam::Array<float>, 3> m_data;
    tomocam::Timer t0;
    t0.start();
    for (int i = 0; i < 3; ++i) {
        auto filename = (base_path / comp1).string();
        m_data[i] = tomocam::tiff::read(filename);
    }
    t0.stop();
    std::cerr << "Time to read data: " << t0.seconds() << "(s)\n";

    // read angles
    std::vector<float> theta(nangles);
    float dtheta = (maxangle - minangle) / static_cast<float>(nangles - 1);
    for (int i = 0; i < nangles; ++i) { theta[i] = minangle + i * dtheta; }
    // convert angles to radians
    for (auto &a : theta) { a = a * M_PI / 180.0f; }

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
    tomocam::PolarGrid<float> grid(theta, nrows, ncols, gamma);
    t0.stop();
    std::cerr << "Time to build a polar grid: " << t0.seconds() << "(s)\n";

    // do the forward projection
    t0.start();
    auto proj = tomocam::forward(m_data, grid, gamma);
    t0.stop();
    std::cerr << "time to do forward projection: " << t0.seconds() << "(s)\n";

    // crop the projection to original size
    tomocam::dims_t crop_dims = {theta.size(), dims.n2, dims.n3};
    t0.start();
    proj = tomocam::crop2d<float>(proj, crop_dims, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << "Time to crop data: " << t0.seconds() << "(s)\n";
    // save data to tiff-stack
    tomocam::tiff::write(output, proj);

    return 0;
}
