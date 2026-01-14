#include <cstdint>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <toml++/toml.hpp>

#include "tomocam.h"

constexpr double PADDING = 1.4142;

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) {
        std::cerr << std::format("Usage: {} TOML input configuration\n", argv[0]);
        return 1;
    }

    auto params = tomocam::ReconParams::from_toml(argv[1]);
    params.print_toml(std::cout);

    tomocam::Timer t0;
    t0.start();
    auto projs = tomocam::tiff::read(params.input_file_path);
    t0.stop();

    std::vector<float> theta;
    std::ifstream angles_file(params.angles_file_path);
    if (!angles_file.is_open()) {
        std::cerr << std::format("Could not open angles file: {}\n",
                                 params.angles_file_path);
        return 1;
    }
    float angle;
    while (angles_file >> angle) { theta.push_back(angle); }
    angles_file.close();
    // convert degrees to radians
    for (auto &a : theta) { a = a * M_PI / 180.0f; }

    // out dims
    auto thickness = params.thickness;
    auto gamma = params.orientation;
    tomocam::dims_t out_dims = {thickness, projs.nrows(), projs.ncols()};

    // padd projections
    t0.start();
    auto projs2 = tomocam::pad2d<float>(projs, PADDING, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cout << std::format("Padded projections to size: {} x {}\n", projs2.nrows(),
                             projs2.ncols());
    std::cout << std::format("Time elapsed for padding: {:.3f} s\n", t0.seconds());

    // create a polar grid
    t0.start();
    auto nrows = projs2.nrows();
    auto ncols = projs2.ncols();
    tomocam::PolarGrid<float> pgrid(theta, nrows, ncols, gamma);
    t0.stop();
    std::cout << std::format("Time elapsed for creating polar grid: {:.3f} s\n",
                             t0.seconds());

    // do the forward projection
    tomocam::dims_t dims = {thickness, nrows, ncols};
    t0.start();
    auto img = tomocam::adjoint(projs2, pgrid, dims, gamma);
    t0.stop();
    std::cout << std::format("Time elapsed for backprojection: {:.3f} s\n",
                             t0.seconds());

    // crop to original size
    t0.start();
    std::array<tomocam::Array<float>, 3> crop_imgs;

    for (size_t i = 0; i < 3; ++i) {
        crop_imgs[i] =
            tomocam::crop3d<float>(img[i], out_dims, tomocam::PadType::SYMMETRIC);
    }
    t0.stop();
    std::cout << std::format("Time elapsed for cropping: {:.3f} s\n", t0.seconds());

    // write output
    auto output = params.output_file_path;
    t0.start();
    tomocam::tiff::write3<float>(output, crop_imgs);
    t0.stop();
    std::cout << std::format("Time elapsed for writing output: {:.3f} s\n",
                             t0.seconds());

    return 0;
}
