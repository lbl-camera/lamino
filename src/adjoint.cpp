#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <toml++/toml.hpp>

#include "array_ops.h"
#include "config.h"
#include "projection.h"
#include "tomocam.h"

constexpr double PADDING =
    1.41421356237; // sqrt(2) to avoid cropping corners of the sample

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " config.toml\n";
        return 1;
    }

    // read parameters from toml file
    auto config = tomocam::read_toml_file(argv[1]);
    auto datasets = tomocam::parse_input_datasets<float>(config);
    auto params = tomocam::ReconParams(config);
    auto output = tomocam::OutputParams(config);

    // print input parameters
    params.print(std::cerr);

    auto &[proj, angles, gamma_ref] = datasets[0];
    float gamma = gamma_ref;

    // record start time
    // record unpadded projection spatial dimensions
    int proj_nrows = static_cast<int>(proj.nrows());
    int proj_ncols = static_cast<int>(proj.ncols());

    // pad projections in 2D — mirrors the pad3d applied to magnetization in
    // forward.cpp so that the polar grid covers the same frequency extent
    tomocam::Timer t0;
    t0.start();
    proj = tomocam::pad2d<float>(proj, static_cast<float>(PADDING - 1),
                                 tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << "Time to pad projections: " << t0.seconds() << "(s)\n";
    std::cerr << "Padded projection dims: [" << proj.nslices() << ", "
              << proj.nrows() << ", " << proj.ncols() << "]\n";

    // compute padded sizes consistent with forward.cpp
    auto padded = [&](int dim) {
        int p = (dim * (PADDING - 1)) / 2;
        return static_cast<size_t>(dim + 2 * p);
    };

    // create polar grid with padded projection dimensions
    t0.start();
    size_t nrows = padded(proj_nrows);
    size_t ncols = padded(proj_ncols);
    tomocam::PolarGrid<float> grid(angles, nrows, ncols, gamma);
    t0.stop();
    std::cerr << "Time to build polar grid: " << t0.seconds() << "(s)\n";
    std::cerr << "Polar grid dimensions: [" << grid.dims().n1 << ", "
              << grid.dims().n2 << ", " << grid.dims().n3 << "]\n";

    tomocam::dims_t recon_dims = params.recon_dims;

    // run adjoint projection
    t0.start();
    auto m_data = tomocam::adjoint(proj, grid, recon_dims, gamma);
    t0.stop();
    std::cerr << "Time to run adjoint: " << t0.seconds() << "(s)\n";

    t0.stop();
    std::cerr << "Time to transpose: " << t0.seconds() << "(s)\n";

    // save results as three separate TIFF files
    auto output_dir = std::filesystem::path(output.filepath).parent_path();
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }
    tomocam::tiff::write3(output.filepath, m_data);
    std::cerr << "Results saved to " << output.filepath << "\n";

    return 0;
}
