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
#include "config.h"
#include "ovf.h"
#include "padding.h"
#include "polar_grid.h"
#include "projection.h"
#include "tiff.h"
#include "timer.h"
#include "tomocam.h"

constexpr double PADDING = tomocam::DEFAULT_PAD_FACTOR;

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
        double gamma_rad;
        std::string output_path;
    };
    std::vector<ProjEntry> projections;
    for (const auto &elem : *proj_array) {
        auto entry = elem.as_table();
        if (!entry) {
            std::cerr << "Invalid [[projections]] entry\n";
            return 1;
        }
        auto gamma_opt = (*entry)["gamma"].value<double>();
        auto filename_opt = (*entry)["filename"].value<std::string>();
        if (!gamma_opt) {
            std::cerr << "[[projections]] entry missing 'gamma' field\n";
            return 1;
        }
        if (!filename_opt) {
            std::cerr << "[[projections]] entry missing 'filename' field\n";
            return 1;
        }
        double gamma_rad = *gamma_opt * M_PI / 180.0;
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

    // if extension is tiff, check all 3 files exist
    bool tiff_format = false;
    if (std::filesystem::path(comp1).extension() == ".tiff" ||
        std::filesystem::path(comp1).extension() == ".tif") {
        tiff_format = true;
        if (!std::filesystem::exists(std::filesystem::path(basedir) / comp1) ||
            !std::filesystem::exists(std::filesystem::path(basedir) / comp2) ||
            !std::filesystem::exists(std::filesystem::path(basedir) / comp3)) {
            std::cerr << "One or more components do not exist in data path: "
                      << basedir << "\n";
            return 1;
        }
    } else if (std::filesystem::path(comp1).extension() == ".ovf") {
        // check that the ovf file exists
        if (!std::filesystem::exists(std::filesystem::path(basedir) / comp1)) {
            std::cerr << "OVF file does not exist in data path: " << basedir << "\n";
            return 1;
        }
    } else {
        std::cerr << "Unsupported component file extension: "
                  << std::filesystem::path(comp1).extension() << "\n";
        return 1;
    }

    // read angles — either from a file or generated from {start, end, num_projs}
    std::vector<double> angles;
    auto angles_node = table["angles"];
    if (auto filename_opt = angles_node["filename"].value<std::string>()) {
        try {
            // same reader (and degree detection) recon uses
            angles = tomocam::read_angles_file<double>(*filename_opt);
        } catch (const std::runtime_error &err) {
            std::cerr << err.what() << "\n";
            return 1;
        }
    } else if (auto start_opt = angles_node["begin"].value<double>()) {
        auto end_opt = angles_node["end"].value<double>();
        auto num_opt = angles_node["num_projs"].value<int>();
        if (!end_opt || !num_opt || *num_opt < 2) {
            std::cerr
                << "angles struct requires 'begin', 'end', and 'num_projs' (>= 2)\n";
            return 1;
        }
        // generate angles linearly spaced between start and end [begin, end)
        int n = *num_opt;
        double start = *start_opt, end = *end_opt;
        angles.resize(n);
        for (int i = 0; i < n; ++i) angles[i] = start + i * (end - start) / n;
        tomocam::to_radians_if_degrees(angles);
    } else {
        std::cerr
            << "angles must specify either 'filename' or {start, end, num_projs}\n";
        return 1;
    }

    auto minangle = *std::min_element(angles.begin(), angles.end());
    auto maxangle = *std::max_element(angles.begin(), angles.end());

    // print parameters
    std::cerr << "----------------------------------------\n";
    std::cerr << "Data path: " << basedir << "\n";
    std::cerr << "Component 1: " << comp1 << "\n";
    std::cerr << "Component 2: " << comp2 << "\n";
    std::cerr << "Component 3: " << comp3 << "\n";
    std::cerr << std::format("Angles: [{:.2f}, {:.2f}] deg with {} steps\n",
                             minangle * 180.0 / M_PI, maxangle * 180.0 / M_PI,
                             angles.size());
    std::cerr << "Output basedir: " << output_basedir << "\n";
    std::cerr << "Projections:\n";
    for (const auto &p : projections) {
        std::cerr << "  gamma=" << (p.gamma_rad * 180.0 / M_PI)
                  << " deg -> " << p.output_path << "\n";
    }
    std::cerr << "----------------------------------------\n";

    // load data
    auto base_path = std::filesystem::path(basedir);
    std::array<std::string, 3> components = {comp1, comp2, comp3};
    std::array<tomocam::Array<float>, 3> m_float;
    tomocam::Timer t0;
    t0.start();
    if (tiff_format) {
        for (int i = 0; i < 3; ++i) {
            auto filename = (base_path / components[i]).string();
            m_float[i] = tomocam::tiff::read(filename);
        }
    } else {
        auto ovf_filename = (base_path / comp1).string();
        m_float = tomocam::ovf::read<float>(ovf_filename);
    }

    t0.stop();
    std::cerr << "Time to read data: " << t0.seconds() << "(s)\n";
    std::cerr << "Data dimensions: [" << m_float[0].nslices() << ", "
              << m_float[0].nrows() << ", " << m_float[0].ncols() << "]\n";

    // convert to double and pad the sample
    t0.start();
    std::array<tomocam::Array<double>, 3> m_data;
    for (int i = 0; i < 3; ++i) {
        m_data[i] = tomocam::pad3d<double>(
            tomocam::array::cast<float, double>(m_float[i]), PADDING,
            tomocam::PadType::SYMMETRIC);
    }
    t0.stop();
    std::cerr << "Time to pad data: " << t0.seconds() << "(s)\n";
    std::cerr << "Padded data dimensions: [" << m_data[0].nslices() << ", "
              << m_data[0].nrows() << ", " << m_data[0].ncols() << "]\n";

    // same (odd) padded size as the volume and as recon's padded projections
    size_t nrows = tomocam::padded_dim(static_cast<size_t>(output_dims[0]), PADDING);
    size_t ncols = tomocam::padded_dim(static_cast<size_t>(output_dims[1]), PADDING);
    if (nrows != m_data[0].nrows() || ncols != m_data[0].ncols()) {
        std::cerr << "Warning: output dims differ from the volume's in-plane "
                     "dims; projections will be offset from the volume center\n";
    }
    tomocam::dims_t crop_dims = {angles.size(), static_cast<size_t>(output_dims[0]),
                                 static_cast<size_t>(output_dims[1])};

    // loop over projections
    for (size_t k = 0; k < projections.size(); ++k) {
        const auto &[gamma_rad, output_path] = projections[k];
        std::cerr << std::format("\n[{}/{}] gamma = {:.1f} deg\n", k + 1,
                                 projections.size(),
                                 gamma_rad * 180.0 / M_PI);

        // build polar grid for this gamma
        t0.start();
        tomocam::PolarGrid<double> grid(angles, nrows, ncols, gamma_rad);
        t0.stop();
        std::cerr << "Time to build polar grid: " << t0.seconds() << "(s)\n";

        // do the forward projection
        t0.start();
        auto proj = tomocam::forward(m_data, grid, gamma_rad, 0.0);
        t0.stop();
        std::cerr << "Time to do forward projection: " << t0.seconds() << "(s)\n";

        // crop the projection to original size
        t0.start();
        proj = tomocam::crop2d<double>(proj, crop_dims, tomocam::PadType::SYMMETRIC);
        t0.stop();
        std::cerr << "Time to crop data: " << t0.seconds() << "(s)\n";

        //  save data to tiff-stack
        tomocam::tiff::write(output_path,
                             tomocam::array::cast<double, float>(proj));
        std::cerr << "Written: " << output_path << "\n";
    }
    // save angles to a text file
    auto angles_path =
        (std::filesystem::path(output_basedir) / "angles.txt").string();
    std::ofstream angles_file(angles_path);
    if (!angles_file.is_open()) {
        std::cerr << "Could not open angles output file: " << angles_path << "\n";
        return 1;
    }
    for (const auto &angle : angles) {
        double deg_angle = angle * 180.0 / M_PI;
        angles_file << std::format("{:.6f}\n", deg_angle);
    }
    std::cout << std::format("Written angles to: {}\n", angles_path);

    return 0;
}
