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

#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "array_ops.h"
#include "tomocam.h"

int main(int argc, char **argv) {

    // check for input file (toml)
    if (argc < 2) {
        std::cerr << std::format("Usage: {} <input.toml>\n", argv[0]);
        std::cerr << "Please see config_template.toml for an example input file.\n";
        tomocam::dump_config();
        return 1;
    }

    // read parameters from toml file
    auto config = tomocam::read_toml_file(argv[1]);
    auto datasets = tomocam::parse_input_datasets<float>(config);
    auto params = tomocam::ReconParams(config);
    auto output = tomocam::OutputParams(config);

    size_t max_iter = params.maxIters;
    float sigma = params.sigma;
    float p = 1.2;
    float tol = params.tol;
    float xtol = params.xtol;
    auto dims = params.recon_dims;

    // print parameters
    params.print(std::cout);

    // Reorder reconstruction dimensions so that the polar components of the
    // cylindrical coordinates vary fastest in memory, improving Fourier transform
    // performance.
    tomocam::dims_t recon_dims = {dims[2], dims[0], dims[1]};

    tomocam::Timer t0;
    t0.start();
    std::array<tomocam::Array<float>, 3> recon;
    switch (params.regularizer) {
        case tomocam::Regularizer::SPLIT_BREGMAN:
            recon = tomocam::MBIR2<float>(datasets, recon_dims, params);
            break;
        case tomocam::Regularizer::qGGMRF:
            recon = tomocam::MBIR3<float>(datasets, recon_dims, params);
            break;
        default: recon = tomocam::MBIR1<float>(datasets, recon_dims, params);
    }

    // transpose to match original input dimensions
    for (size_t i = 0; i < 3; ++i) {
        recon[i] = tomocam::array::transpose(recon[i], {1, 2, 0});
    }

    t0.stop();
    double elapsed = t0.seconds();
    if (elapsed > 3600) {
        int hours = static_cast<int>(elapsed / 3600);
        int minutes = static_cast<int>((elapsed - hours * 3600) / 60);
        double seconds = elapsed - hours * 3600 - minutes * 60;
        std::cout << std::format("Reconstruction completed in {} hours, {} minutes, "
                                 "and {:.2f} seconds.\n",
                                 hours, minutes, seconds);
    } else if (elapsed > 60) {
        int minutes = static_cast<int>(elapsed / 60);
        double seconds = elapsed - minutes * 60;
        std::cout << std::format(
            "Reconstruction completed in {} minutes and {:.2f} seconds.\n", minutes,
            seconds);
    } else {
        std::cout << std::format("Reconstruction completed in {:.2f} seconds.\n",
                                 elapsed);
    }

    // save result to tiff
    auto base_dir = std::filesystem::path(output.filepath).parent_path();
    if (!std::filesystem::exists(base_dir)) {
        std::filesystem::create_directories(base_dir);
    }

    if (output.has_format("tiff")) { tomocam::tiff::write3(output.filepath, recon); }
    if (output.has_format("vti")) {
        tomocam::vti::write_vectors(output.filepath, recon);
    }
    return 0;
}
