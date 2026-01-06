/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals from the
 * U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 * IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 * the U.S. Government has been granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 * to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#include <cstdint>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

#include "tomocam.h"

constexpr double PADDING = 1.4142;

int main(int argc, char **argv) {

    // read input arguments
    if (argc < 2) {
        std::cerr << std::format("Usage: {} TOML input configuration\n", argv[0]);
        return 1;
    }
    auto params = tomocam::ReconParams::from_toml(argv[1]);

    tomocam::Timer t0;
    t0.start();
    auto projs = tomocam::tiff::read(params.input_file_path);
    t0.stop();
    std::cerr << std::format("Projections size: ({}, {})\n", projs.nrows(),
                             projs.ncols());
    std::cerr << std::format("No. of projections: {}\n", projs.nslices());

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
    std::cerr << std::format("Number of angles: {}\n", theta.size());

    // out dims
    auto thickness = params.thickness;
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
    std::cerr << std::format("Polar grid size: ({}, {})\n", nrows, ncols);

    // do the forward projection
    tomocam::dims_t dims = {thickness, nrows, ncols};
    std::cerr << std::format("Backprojecting to size: ({}, {}, {})\n", dims.n1,
                             dims.n2, dims.n3);
    t0.start();
    auto img = tomocam::adjoint<float>(projs2, pgrid, dims, params.orientation);
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
    std::string output = params.output_file_path;
    t0.start();
    tomocam::tiff::write3<float>(output, crop_imgs);
    t0.stop();

    return 0;
}
