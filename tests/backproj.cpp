#include <cstdint>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <toml++/toml.hpp>

#include "array_ops.h"
#include "tomocam.h"

constexpr double PADDING = 1.4142;

inline size_t pad(size_t n) {
    size_t npad = 2 * (size_t(std::ceil(n * (PADDING - 1) / 2)));
    return n + npad;
}

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
    auto out_dims = params.recon_dims;

    // shift dims so that radial direction is the last dimension
    tomocam::dims_t recon_dims = {out_dims[2], out_dims[0], out_dims[1]};
    // pad output dimensions
    recon_dims.n1 = pad(recon_dims.n1);
    recon_dims.n2 = pad(recon_dims.n2);
    recon_dims.n3 = pad(recon_dims.n3);

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
    auto img = tomocam::adjoint(projs2, pgrid, recon_dims, gamma);
    t0.stop();
    std::cout << std::format("Time elapsed for backprojection: {:.3f} s\n",
                             t0.seconds());

    t0.start();
    // transpose to (n2, n3, n1)
    for (size_t i = 0; i < 3; ++i) {
        img[i] = tomocam::array::transpose<float>(img[i], {1, 2, 0});
    }
    t0.stop();
    std::cout << std::format("Time elapsed for transposing: {:.3f} s\n",
                             t0.seconds());

    // crop to original size
    t0.start();
    std::array<tomocam::Array<float>, 3> crop_imgs;
    for (size_t i = 0; i < 3; ++i) {
        crop_imgs[i] = tomocam::crop3d<float>(img[i], params.recon_dims,
                                              tomocam::PadType::SYMMETRIC);
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
