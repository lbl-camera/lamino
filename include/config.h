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

#ifndef CONFIG_H
#define CONFIG_H

#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <toml++/toml.h>
#include <tuple>
#include <vector>

#include "array.h"
#include "tiff.h"

namespace tomocam {

    // Function to read and parse a TOML file
    inline toml::table read_toml_file(const std::string &filepath) {
        if (!std::filesystem::exists(filepath)) {
            throw std::runtime_error(
                std::format("TOML file does not exist: {}", filepath));
        }

        try {
            return toml::parse_file(filepath);
        } catch (const toml::parse_error &err) {
            throw std::runtime_error(std::format(
                "Failed to parse TOML file '{}': {}", filepath, err.description()));
        }
    }

    // Dataset type: (projections, angles, gamma)
    template <typename T>
    using Dataset_t = std::tuple<Array<T>, std::vector<T>, T>;

    // Function to read angles from a text file
    template <typename T>
    inline std::vector<T> read_angles_file(const std::string &filepath) {
        std::ifstream fp(filepath);
        if (!fp.is_open()) {
            throw std::runtime_error(
                std::format("Could not open angles file: {}", filepath));
        }

        std::vector<T> angles;
        T angle;
        while (fp >> angle) { angles.push_back(angle); }
        fp.close();

        // Convert to radians if necessary
        auto max_angle = *std::max_element(angles.begin(), angles.end());
        if (std::abs(max_angle) > 2 * M_PI) {
            for (auto &a : angles) { a = a * M_PI / (T)180.0; }
        }
        return angles;
    }

    // Function to parse input datasets from TOML config
    template <typename T>
    [[nodiscard]] std::vector<Dataset_t<T>>
    parse_input_datasets(const toml::table &config) {

        auto input_array = config["input"].as_array();
        if (!input_array) {
            throw std::runtime_error("Missing [[input]] array in TOML file");
        }

        std::vector<Dataset_t<T>> datasets;

        for (auto &elem : *input_array) {
            auto input_table = elem.as_table();
            if (!input_table) {
                throw std::runtime_error("Invalid [[input]] entry");
            }

            auto filename = (*input_table)["filename"].value<std::string>();
            auto angles_file = (*input_table)["angles"].value<std::string>();
            auto gamma = (*input_table)["gamma"].value<T>();

            if (!filename || !angles_file || !gamma) {
                throw std::runtime_error("[[input]] entry must have 'filename', "
                                         "'angles', and 'gamma' fields");
            }

            auto projs = tomocam::tiff::read(*filename);
            auto angles = read_angles_file<T>(*angles_file);
            T gamm_radians = (*gamma) * M_PI / (T)180.0;
            datasets.push_back(
                std::make_tuple(std::move(projs), std::move(angles), gamm_radians));
        }
        return datasets;
    }

    // Reconstruction parameters
    struct ReconParams {
        std::array<size_t, 3> recon_dims = {51, 511,
                                            511}; // Reconstruction dimensions
        size_t maxIters = 50;                     // Maximum number of iterations
        float sigma = 1000.0f;                    // Regularization parameter
        float tol = 1e-5f;                        // Tolerance for convergence
        float xtol = 1e-5f;                       // Tolerance for solution change

        ReconParams() = default;
        ReconParams(const toml::table &config) {

            // Read [recon_params] section
            auto recon = config["recon_params"];
            if (!recon) {
                throw std::runtime_error(
                    "Missing [recon_params] section in config file");
            }
            maxIters = recon["max_iters"].value_or<size_t>(50);
            const auto *dims = recon["recon_dims"].as_array();
            if (dims && dims->size() == 3) {
                for (size_t i = 0; i < 3; ++i) {
                    recon_dims[i] = (*dims)[i].value_or<size_t>(0);
                }
            } else {
                throw std::runtime_error("[recon_params] 'recon_dims' must be an "
                                         "array of three integers");
            }
            // if reocn dims[0]/dims[2] > 0.15, remind user that since this is
            // laminography, the thickness is expected to be much smaller than the
            // in-plane dimensions, and this may lead to increased memory usage
            float ratio = static_cast<float>(recon_dims[0]) /
                          static_cast<float>(recon_dims[2]);
            if (ratio > 0.15f) {
                std::cerr << std::format(
                    "Warning: recon_dims[0] / recon_dims[2] = {:.2f} > 0.15. "
                    "In laminography, the thickness (recon_dims[0]) "
                    "is expected to be much smaller than the in-plane "
                    "dimensions (recon_dims[1], recon_dims[2]). "
                    "This may lead to increased memory usage.\n",
                    ratio);
            }
            sigma = recon["sigma"].value_or<float>(1000.0f);
            tol = recon["tol"].value_or<float>(1e-5f);
            xtol = recon["xtol"].value_or<float>(1e-5f);
        }

        void print(std::ostream &os) const {

            std::string dims = "[";
            for (const auto &d : recon_dims) { dims += std::to_string(d) + ", "; }
            dims = dims.substr(0, dims.size() - 2) + "]";
            auto config = toml::table{
                {"recon_params",
                 toml::table{{"max_iters", static_cast<int64_t>(maxIters)},
                             {"sigma", sigma},
                             {"recon_dims", dims},
                             {"tol", tol},
                             {"xtol", xtol}}}};
            os << config << std::endl;
        }
    };

    // Output parameters
    struct OutputParams {
        std::string filepath;

        OutputParams() = default;
        OutputParams(const toml::table &config) {

            // Read [output] section
            auto output = config["output"];
            if (!output) {
                throw std::runtime_error("Missing [output] section in config file");
            }
            filepath = output["filename"].value_or<std::string>("./recon.tiff");
        }
    };

    // Function to dump an example configuration file
    inline void dump_config(const std::string &filepath = "config.toml") {
        std::ofstream outfile(filepath);
        if (!outfile.is_open()) {
            throw std::runtime_error(
                std::format("Could not open file for writing: {}", filepath));
        }

        outfile << "[[input]]\n";
        outfile << "filename = \"/path/to/gamma0_stack.tiff\"\n";
        outfile << "angles = \"/path/to/gamma0_angles.txt\"\n";
        outfile << "gamma = 0\n";
        outfile << "\n";
        outfile << "[[input]]\n";
        outfile << "filename = \"/path/to/gamma45_stack.tiff\"\n";
        outfile << "angles = \"/path/to/gamma45_angles.txt\"\n";
        outfile << "gamma = 45\n";
        outfile << "\n";
        outfile << "[output]\n";
        outfile << "filename = \"output.tiff\"\n";
        outfile << "\n";
        outfile << "[recon_params]\n";
        outfile << "max_iters = 50\n";
        outfile << "sigma = 1000.0\n";
        outfile << "tol = 1e-5\n";
        outfile << "xtol = 1e-5\n";
        outfile << "recon_dims = [51, 511, 511]\n";
        outfile << "\n";

        outfile.close();
    }
} // namespace tomocam
#endif // CONFIG_H
