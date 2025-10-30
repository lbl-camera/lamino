#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <ostream>
#include <string>
using json = nlohmann::json;

#include "tiff.h"
#include "timer.h"
#include "tomocam.h"

constexpr double PADDING = 1.4142;

void usage(char **argv) {
    std::cout << "Usage: " << argv[0] << " JSON input configuration\n";
    exit(0);
}

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) { usage(argv); }

    // read data
    std::ifstream json_file(argv[1]);
    if (!json_file.is_open()) {
        std::cerr << "error: can't open JSON file: " << argv[1] << '\n';
        return 1;
    }

    // get input data
    json config = json::parse(json_file);
    std::string filename = config["filename"];
    std::string output = config["output"];
    float gamma = config["gamma"];

    tomocam::Timer t0;
    t0.start();
    auto sample = tomocam::tiff::read(filename);
    t0.stop();
    std::cerr << "Time to read data: " << t0.seconds() << "(s)\n";

    std::vector<float> theta(141, 0);
    for (int i = 0; i < theta.size(); i++) {
        int deg = (-1 * (theta.size() / 2)) + i + 90;
        theta[i] = deg * M_PI / 180.F;
    }

    // create a polar grid
    t0.start();
    // oversample polar-grid
    auto nrows = 2 * (static_cast<uint64_t>(PADDING * sample.nrows()) / 2);
    auto ncols = 2 * (static_cast<uint64_t>(PADDING * sample.ncols()) / 2);
    tomocam::PolarGrid<float> pgrid(theta, nrows, ncols);
    t0.stop();
    std::cerr << "Time to build a polar grid: " << t0.seconds() << "(s)\n";

    // do the forward projection
    t0.start();
    auto proj = tomocam::forward(sample, pgrid, gamma);
    t0.stop();
    std::cerr << "time to do forward projection: " << t0.seconds() << "(s)\n";

    // save data to tiff-stack
    tomocam::tiff::write(output, proj);

    return 0;
}
