#include <cstdint>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <toml++/toml.hpp>

#include "padding.h"
#include "polar_grid.h"
#include "recon_params.h"
#include "tiff.h"
#include "timer.h"
#include "tomocam.h"

constexpr double PADDING = 1.4142;

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) {
        std::cerr << std::format("Usage: {} TOML input configuration\n", argv[0]);
        return 1;
    }

    // read data
    std::ifstream toml_file(argv[1]);
    if (!toml_file.is_open()) {
        std::cerr << std::format("Usage: {} TOML input configuration\n", argv[0]);
        return 1;
    }
    auto params = ReconParams::from_toml(toml_file);
    toml_file.close();

    params.print_toml(std::cout);
    std::exit(0);

    tomocam::Timer t0;
    t0.start();
    auto projs = tomocam::tiff::read(filename);
    t0.stop();

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

    // out dims
    tomocam::dims_t out_dims = {thickness, projs.nrows(), projs.ncols()};

    // padd projections
    t0.start();
    auto projs2 = tomocam::pad2d<float>(projs, PADDING, tomocam::PadType::SYMMETRIC);
    t0.stop();

    // create a polar grid
    t0.start();
    auto nrows = projs2.nrows();
    auto ncols = projs2.ncols();
    tomocam::PolarGrid<float> pgrid(theta, nrows, ncols);
    t0.stop();

    // do the forward projection
    tomocam::dims_t dims = {thickness, nrows, ncols};
    t0.start();
    auto img = tomocam::adjoint(projs2, pgrid, dims, gamma);
    t0.stop();

    // crop to original size
    t0.start();
    std::array<tomocam::Array<float>, 3> crop_imgs;

    for (size_t i = 0; i < 3; ++i) {
        crop_imgs[i] =
            tomocam::crop3d<float>(img[i], out_dims, tomocam::PadType::SYMMETRIC);
    }
    t0.stop();

    // write output
    t0.start();
    tomocam::tiff::write3<float>(output, crop_imgs);
    t0.stop();

    return 0;
}
