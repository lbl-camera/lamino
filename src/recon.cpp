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

#include <format>
#include <fstream>
#include <iostream>
#include <string>

#include "tomocam.h"

int main(int argc, char **argv) {

    // paser JSON input
    if (argc < 2) {
        std::cerr << std::format("Usage: {} <input.toml>\n", argv[0]);
        std::cerr << "Please see config_template.toml for an example input file.\n";
        tomocam::ReconParams::dump_config();
        return 1;
    }

    std::string input_file = argv[1];
    std::ifstream ifs(input_file);
    if (!ifs.is_open()) {
        std::cerr << "Error: Could not open input file " << input_file << "\n";
        return 1;
    }
    auto params = tomocam::ReconParams::from_toml(input_file);
    auto projs = tomocam::tiff::read(params.input_file_path);

    // read projection angles
    std::ifstream fp(params.angles_file_path);
    if (!fp.is_open()) {
        std::cerr << std::format("Error: Could not open angles file {}\n",
                                 params.angles_file_path);
        return 1;
    }
    std::vector<float> angles;
    float angle;
    while (fp >> angle) { angles.push_back(angle); }
    fp.close();

    // find keys in json, else set default values
    float gamma = params.orientation;
    size_t max_iter = params.maxIters;
    float sigma = params.sigma;
    float p = 1.2;
    float tol = params.tol;
    float xtol = params.xtol;
    size_t thickness = params.thickness;

    // log parameters
    params.print_toml(std::cout);
    exit(1);

    // ensure angles are in radians
    auto max_angle = *std::max_element(angles.begin(), angles.end());
    if (std::abs(max_angle) > 2 * M_PI) {
        for (auto &angle : angles) { angle = angle * M_PI / 180.0f; }
    }

    // set reconstruction dimensions
    tomocam::dims_t rec_dims = {thickness, projs.nrows(), projs.ncols()};

    tomocam::Timer t0;
    t0.start();
    auto recon = tomocam::MBIR<float>(projs, angles, gamma, rec_dims, max_iter,
                                      sigma, p, tol, xtol);
    t0.stop();
    std::cout << std::format("Reconstruction completed in {:.2f} seconds.\n",
                             t0.seconds());

    // save result to tiff
    std::array<std::string, 3> tifnames = {"recon_x.tif", "recon_y.tif",
                                           "recon_z.tif"};
    for (size_t i = 0; i < 3; i++) { tomocam::tiff::write(tifnames[i], recon[i]); }
    tomocam::vti::write_vectors("recon.vti", recon);
    return 0;
}
