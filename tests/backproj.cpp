#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <toml++/toml.h>

#include "padding.h"
#include "polar_grid.h"
#include "tiff.h"
#include "timer.h"
#include "tomocam.h"

constexpr double PADDING = 1.4142;

void usage(char **argv) {
    std::cout << "Usage: " << argv[0] << " TOML input configuration" << std::endl;
    exit(0);
}

int main(int argc, char **argv) {

    // parse TOML input
    if (argc < 2) usage(argv);
    std::string config_file = argv[1];
    toml::table config;
    try {
        config = toml::parse_file(config_file);
    } catch (const toml::parse_error &err) {
        std::cerr << "Error parsing TOML file:\n" << err << "\n";
        return 1;
    }

    std::string filename = config["input"]["filename"].value_or("");
    if (!std::filesystem::exists(filename)) {
        std::cerr << std::format("Error: projection file '{}' does not exist\n",
                                 filename);
        return 1;
    }
    std::string angles = config["input"]["angles"].value_or("");
    if (!std::filesystem::exists(angles)) {
        std::cerr << std::format("Error: angles file '{}' does not exist\n", angles);
        return 1;
    }
    std::string output = config["output"]["filename"].value_or("output.tiff");
    size_t thickness = config["recon_params"]["thickness"].value_or(21);

    tomocam::Timer t0;
    t0.start();
    auto projs = tomocam::tiff::read(filename);
    t0.stop();

    size_t nproj = projs.nslices();
    size_t nrows = projs.nrows();
    size_t ncols = projs.ncols();
    tomocam::dims_t recon_dims = {thickness, nrows, ncols};

    std::cerr << std::format("Time to read data: {}(s)\n", t0.seconds());
    std::cerr << std::format("Projections size: ({}, {}, {})\n", nproj, nrows,
                             ncols);
    std::cerr << std::format("Reconstruction size: ({}, {}, {})\n", recon_dims.n1,
                             recon_dims.n2, recon_dims.n3);

    std::vector<float> theta;
    std::ifstream angles_file(angles);
    if (!angles_file.is_open()) {
        std::cerr << std::format("Could not open angles file: {}\n", angles);
        return 1;
    }
    float angle;
    while (angles_file >> angle) { theta.push_back(angle); }
    angles_file.close();
    // convert degrees to radians
    for (auto &a : theta) { a = a * M_PI / 180.0f; }
    if (theta.size() != nproj) {
        std::cerr << std::format(
            "Number of angles ({}) does not match number of projections ({})\n",
            theta.size(), nproj);
        return 1;
    }

    // pad projections

    t0.start();
    auto projs2 = tomocam::pad2d<float>(projs, PADDING, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << std::format("Time to pad projections: {}(s)\n", t0.seconds());
    std::cerr << std::format("Padded projections size: ({}, {})\n", projs2.ncols(),
                             projs2.nrows());

    // create a polar grid
    t0.start();
    tomocam::PolarGrid<float> pgrid(theta, projs2.nrows(), projs2.ncols());
    t0.stop();
    std::cerr << std::format("Polar grid size: ({}, {})\n", projs2.nrows(),
                             projs2.ncols());
    std::cerr << std::format("Time to create polar grid: {}(s)\n", t0.seconds());

    float gamma = 0.0f;
    t0.start();
    auto img = tomocam::backproj(projs2, pgrid, recon_dims, gamma);
    t0.stop();
    std::cerr << std::format("Time to backproject: {}(s)\n", t0.seconds());

    // crop the image to original size

    tomocam::dims_t crop_dims = {thickness, projs.nrows(), projs.ncols()};
    t0.start();
    img = tomocam::crop2d<float>(img, crop_dims, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << std::format("Cropped image size: ({}, {}, {})\n", img.nslices(),
                             img.nrows(), img.ncols());
    std::cerr << std::format("Time to crop image: {}(s)\n", t0.seconds());

    // save data to tiff-stack
    tomocam::tiff::write(output, img);

    return 0;
}
