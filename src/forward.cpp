#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include <toml++/toml.hpp>

#include "array_ops.h"
#include "tomocam.h"

constexpr double PADDING =
    1.41421356237; // sqrt(2) to avoid cropping corners of the sample

int main(int argc, char **argv) {

    // sanity check
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " config.toml\n";
        return 1;
    }
    // get input data
    toml::table table;
    try {
        table = toml::parse_file(argv[1]);
    } catch (const toml::parse_error &err) {
        std::cerr << std::format("Parsing failed:\n{}\n", err.description());
        return 1;
    }

    // get parameters from config file
    auto paths = table["input"];
    auto basedir = paths["data_path"].value<std::string>().value();
    auto comp1 = paths["component_1"].value<std::string>().value();
    auto comp2 = paths["component_2"].value<std::string>().value();
    auto comp3 = paths["component_3"].value<std::string>().value();

    auto angles_file = table["angles"]["filename"].value<std::string>().value();

    auto output_basedir = table["output"]["basedir"].value<std::string>().value();
    std::vector<int> output_dims;
    auto dim_arr = table["output"]["dims"].as_array();
    if (dim_arr && dim_arr->size() == 2) {
        for (const auto &elem : *dim_arr) {
            if (auto val = elem.value<int>()) {
                output_dims.push_back(*val);
            } else {
                std::cerr << "Invalid output dimension value in config file\n";
                return 1;
            }
        }
    } else {
        std::cerr << "Output dimensions must be an array of two integers\n";
        return 1;
    }

    // parse [[projections]] array: each entry has gamma (degrees) and filename
    auto proj_array = table["projections"].as_array();
    if (!proj_array || proj_array->empty()) {
        std::cerr << "Missing or empty [[projections]] array in config file\n";
        return 1;
    }
    struct ProjEntry {
        float gamma_rad;
        std::string output_path;
    };
    std::vector<ProjEntry> projections;
    for (const auto &elem : *proj_array) {
        auto entry = elem.as_table();
        if (!entry) {
            std::cerr << "Invalid [[projections]] entry\n";
            return 1;
        }
        auto gamma_opt = (*entry)["gamma"].value<float>();
        auto filename_opt = (*entry)["filename"].value<std::string>();
        if (!gamma_opt) {
            std::cerr << "[[projections]] entry missing 'gamma' field\n";
            return 1;
        }
        if (!filename_opt) {
            std::cerr << "[[projections]] entry missing 'filename' field\n";
            return 1;
        }
        float gamma_rad = *gamma_opt * static_cast<float>(M_PI) / 180.0f;
        auto output_path =
            (std::filesystem::path(output_basedir) / *filename_opt).string();
        projections.push_back({gamma_rad, output_path});
    }

    // sanity checks
    std::filesystem::path out_basedir(output_basedir);
    // create output directory if it does not exist
    if (!std::filesystem::exists(out_basedir)) {
        std::filesystem::create_directories(out_basedir);
    }

    // check if components exist
    if (!std::filesystem::exists(std::filesystem::path(basedir) / comp1) ||
        !std::filesystem::exists(std::filesystem::path(basedir) / comp2) ||
        !std::filesystem::exists(std::filesystem::path(basedir) / comp3)) {
        std::cerr << "One or more components do not exist in data path: " << basedir
                  << "\n";
        return 1;
    }

    // read angles from text file
    std::ifstream angles_stream(angles_file);
    if (!angles_stream.is_open()) {
        std::cerr << "Could not open angles file: " << angles_file << "\n";
        return 1;
    }
    std::vector<float> angles;
    float angle;
    while (angles_stream >> angle) { angles.push_back(angle); }
    if (angles.empty()) {
        std::cerr << "No angles found in file: " << angles_file << "\n";
        return 1;
    }

    auto minangle = *std::min_element(angles.begin(), angles.end());
    auto maxangle = *std::max_element(angles.begin(), angles.end());
    if (std::abs(minangle) > M_PI || std::abs(maxangle) > M_PI) {
        for (auto &a : angles) { a = a * M_PI / 180.0f; }
    }

    // print parameters
    std::cerr << "----------------------------------------\n";
    std::cerr << "Data path: " << basedir << "\n";
    std::cerr << "Component 1: " << comp1 << "\n";
    std::cerr << "Component 2: " << comp2 << "\n";
    std::cerr << "Component 3: " << comp3 << "\n";
    std::cerr << "Angles: [" << minangle << ", " << maxangle << "] with "
              << angles.size() << " steps\n";
    std::cerr << "Output basedir: " << output_basedir << "\n";
    std::cerr << "Projections:\n";
    for (const auto &p : projections) {
        std::cerr << "  gamma=" << (p.gamma_rad * 180.0f / static_cast<float>(M_PI))
                  << " deg -> " << p.output_path << "\n";
    }
    std::cerr << "----------------------------------------\n";

    // load data
    auto base_path = std::filesystem::path(basedir);
    std::array<std::string, 3> components = {comp1, comp2, comp3};
    std::array<tomocam::Array<float>, 3> m_data;
    tomocam::Timer t0;
    t0.start();
    for (int i = 0; i < 3; ++i) {
        auto filename = (base_path / components[i]).string();
        m_data[i] = tomocam::tiff::read(filename);
    }

    t0.stop();
    std::cerr << "Time to read data: " << t0.seconds() << "(s)\n";
    std::cerr << "Data dimensions: [" << m_data[0].nslices() << ", "
              << m_data[0].nrows() << ", " << m_data[0].ncols() << "]\n";

    // pad the sample
    t0.start();
    for (int i = 0; i < 3; ++i) {
        m_data[i] = tomocam::pad3d<float>(m_data[i], PADDING - 1,
                                          tomocam::PadType::SYMMETRIC);
    }
    t0.stop();
    std::cerr << "Time to pad data: " << t0.seconds() << "(s)\n";
    std::cerr << "Padded data dimensions: [" << m_data[0].nslices() << ", "
              << m_data[0].nrows() << ", " << m_data[0].ncols() << "]\n";

    auto padded = [&](int dim) {
        int p = (dim * (PADDING - 1)) / 2;
        return static_cast<size_t>(dim + 2 * p);
    };

    size_t nrows = static_cast<size_t>(padded(output_dims[0]));
    size_t ncols = static_cast<size_t>(padded(output_dims[1]));
    tomocam::dims_t crop_dims = {angles.size(), static_cast<size_t>(output_dims[0]),
                                 static_cast<size_t>(output_dims[1])};

    // loop over projections
    for (size_t k = 0; k < projections.size(); ++k) {
        const auto &[gamma_rad, output_path] = projections[k];
        std::cerr << std::format("\n[{}/{}] gamma = {:.1f} deg\n", k + 1,
                                 projections.size(),
                                 gamma_rad * 180.0f / static_cast<float>(M_PI));

        // build polar grid for this gamma
        t0.start();
        tomocam::PolarGrid<float> grid(angles, nrows, ncols, gamma_rad);
        t0.stop();
        std::cerr << "Time to build polar grid: " << t0.seconds() << "(s)\n";

        // do the forward projection
        t0.start();
        auto proj = tomocam::forward(m_data, grid, gamma_rad);
        t0.stop();
        std::cerr << "Time to do forward projection: " << t0.seconds() << "(s)\n";

        // crop the projection to original size
        t0.start();
        proj = tomocam::crop2d<float>(proj, crop_dims, tomocam::PadType::SYMMETRIC);
        t0.stop();
        std::cerr << "Time to crop data: " << t0.seconds() << "(s)\n";

        //  save data to tiff-stack
        tomocam::tiff::write(output_path, proj);
        std::cerr << "Written: " << output_path << "\n";
    }

    return 0;
}
