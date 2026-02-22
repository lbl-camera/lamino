#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <toml++/toml.h>

#include "padding.h"
#include "tiff.h"
#include "timer.h"
#include "tomocam.h"

constexpr double PADDING = 1.4142;

void usage() {
    std::cerr << std::format("Usage: <config.toml>\n");
    std::cerr << "Example config.toml:\n";
    std::cerr << R"(
    filename = "/path/do/volume/data.tiff"
    angles = "/path/to/angles/angles.txt"
    output = "/path/to/output/proj.tiff"
    gamma = 1.0
)";
    exit(1);
}

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) { usage(); }
    const std::string toml_file = argv[1];
    // read config file
    std::ifstream ifs(toml_file);
    if (!ifs.is_open()) {
        std::cerr << std::format("error: can't open config file: {}\n", toml_file);
        return 1;
    }
    auto config = toml::parse(ifs);

    // get input data
    std::string filename = config["filename"].value_or("");
    if (filename.empty()) {
        std::cerr << "error: filename is not specified in config file\n";
        return 1;
    }
    if (!std::filesystem::exists(filename)) {
        std::cerr << std::format("error: file does not exist: {}\n", filename);
        return 1;
    }

    std::string angles_file = config["angles"].value_or("");
    if (angles_file.empty()) {
        std::cerr << "error: angles file is not specified in config file\n";
        return 1;
    }
    if (!std::filesystem::exists(angles_file)) {
        std::cerr << std::format("error: angles file does not exist: {}\n",
                                 angles_file);
        return 1;
    }

    std::string output = config["output"].value_or("");
    if (output.empty()) {
        output = std::filesystem::path(filename).stem().string() + "_proj.tiff";
    }
    float gamma = config["gamma"].value_or(0.f);

    tomocam::Timer t0;
    t0.start();
    auto sample = tomocam::tiff::read(filename);
    t0.stop();
    std::cerr << std::format("Time to read data: {} (s)\n", t0.seconds());

    // read angles
    std::vector<float> theta;
    std::ifstream fp(angles_file);
    if (!fp.is_open()) {
        std::cerr << std::format("error: can't open angles file: {}\n", angles_file);
        return 1;
    }
    float angle;
    while (fp >> angle) { theta.push_back(angle); }
    fp.close();
    // convert angles to radians
    for (auto &a : theta) { a = a * M_PI / 180.0f; }

    // pad the sample
    t0.start();
    auto sample2 =
        tomocam::pad3d<float>(sample, PADDING, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << std::format("Time to pad data: {} (s)\n", t0.seconds());

    // create a polar grid
    t0.start();
    // oversample polar-grid
    auto nrows = sample2.nrows();
    auto ncols = sample2.ncols();
    tomocam::PolarGrid<float> pgrid(theta, nrows, ncols);
    t0.stop();
    std::cerr << std::format("Time to build a polar grid: {} (s)\n", t0.seconds());

    // do the forward projection
    t0.start();
    auto proj = tomocam::forward<float>(sample2, pgrid, gamma);
    t0.stop();
    std::cerr << std::format("Time to do forward projection: {} (s)\n",
                             t0.seconds());

    // crop the projection to original size
    tomocam::dims_t crop_dims = {theta.size(), sample.nrows(), sample.ncols()};
    t0.start();
    auto proj2 =
        tomocam::crop2d<float>(proj, crop_dims, tomocam::PadType::SYMMETRIC);
    t0.stop();
    std::cerr << std::format("Time to crop data: {} (s)\n", t0.seconds());
    // save data to tiff-stack
    tomocam::tiff::write(output, proj2);

    return 0;
}
