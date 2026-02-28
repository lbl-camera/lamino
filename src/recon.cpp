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
#ifdef DEBUG
    params.print(std::cout);
#endif

    // set reconstruction dimensions
    tomocam::dims_t recon_dims = {dims[2], dims[0], dims[1]};

    tomocam::Timer t0;
    t0.start();
    auto recon = tomocam::MBIR2<float>(datasets, recon_dims, params);
    t0.stop();
    std::cout << std::format("Reconstruction completed in {:.2f} seconds.\n",
                             t0.seconds());

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
