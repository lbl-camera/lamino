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

    auto config = toml::parse_file(argv[1]);
    auto params = tomocam::ReconParams(config);
    auto output = tomocam::OutputParams(config);
    params.print(std::cout);

    // read projections
    tomocam::Timer t0;
    t0.start();
    auto dset = tomocam::parse_input_datasets<float>(config);
    t0.stop();
    std::cout << std::format("Time elapsed for reading dataset: {:.3f} s\n",
                             t0.seconds());
    // out dims
    auto &[projs, theta, gamma_ref] = dset[0];
    auto gamma = gamma_ref;
    tomocam::dims_t out_dims = params.recon_dims;

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
    t0.start();
    auto img = tomocam::adjoint(projs2, pgrid, out_dims, gamma);
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
    t0.start();
    tomocam::tiff::write3<float>(output.filepath, crop_imgs);
    t0.stop();
    std::cout << std::format("Time elapsed for writing output: {:.3f} s\n",
                             t0.seconds());

    return 0;
}
