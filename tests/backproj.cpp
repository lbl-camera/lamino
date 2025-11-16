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
    std::cout << "Usage: " << argv[0] << " JSON input configuration" << std::endl;
    exit(0);
}

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) usage(argv);

    // read data
    std::ifstream json_file(argv[1]);
    if (!json_file.is_open()) {
        std::cerr << "error: can't open JSON file: " << argv[1] << std::endl;
        return 1;
    }

    // get input data
    json config = json::parse(json_file);
    std::string filename = config["filename"];
    std::string output = config["output"];
    float gamma = config["gamma"];

    tomocam::Timer t0;
    t0.start();
    auto projs = tomocam::tiff::read(filename);
    t0.stop();
    std::cerr << "Time to read data: " << t0.seconds() << "(s)\n";

    std::vector<float> theta(141, 0);
    for (int i = 0; i < theta.size(); i++) {
        int deg = -1 * (theta.size() / 2) + i + 90;
        theta[i] = deg * M_PI / 180.f;
    }

    // create a polar grid
    t0.start();
    auto nrows = 2 * (static_cast<uint64_t>(PADDING * projs.nrows()) / 2);
    auto ncols = 2 * (static_cast<uint64_t>(PADDING * projs.nrows()) / 2);
    tomocam::PolarGrid<float> pgrid(theta, nrows, ncols);
    t0.stop();
    std::cerr << "Time to build a polar grid: " << t0.seconds() << "(s)\n";

    // do the forward projection
    tomocam::dims_t dims = {20, projs.nrows(), projs.ncols()};
    t0.start();
    auto img = tomocam::backward(projs, pgrid, dims, gamma, false);
    t0.stop();
    std::cerr << "time to do back-projection: " << t0.seconds() << "(s)\n";

    // save data to tiff-stack
    tomocam::tiff::write(output, img);

    return 0;
}
