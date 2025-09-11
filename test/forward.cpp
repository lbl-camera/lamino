#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <ostream>
#include <string>
using json = nlohmann::json;

#include "tiff.h"
#include "timer.h"
#include "tomocam.h"

void usage(char **argv) {
    std::cout << "Usage: " << argv[0] << " JSON input configuration" << std::endl;
    exit(0);
}

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2)
        usage(argv);

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

    tomocam::Timer t0;
    t0.start();
    auto sample = tomocam::tiff::read(filename);
    t0.stop();
    std::cerr << "Time to read data: " << t0.ms() << std::endl;

    std::vector<float> theta(141, 0);
    for (int i = 0; i < theta.size(); i++) {
        int deg = -1 * (theta.size() / 2) + i;
        theta[i] = deg * M_PI / 180.f;
    }

    // create a polar grid
    t0.start();
    tomocam::PolarGrid<float> pgrid(theta, sample.nrows(), sample.ncols());
    t0.stop();
    std::cerr << "Time to build a polar grid: " << t0.ms() << std::endl;

    // do the forward projection
    t0.start();
    auto proj = tomocam::forward(sample, pgrid);
    t0.stop();
    std::cerr << "time to do forward projection: " << t0.ms() << std::endl;

    // save data to tiff-stack
    tomocam::tiff::write(output, proj);

    return 0;
}
