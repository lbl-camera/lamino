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

#ifndef RECON_PARAMS_H
#define RECON_PARAMS_H

#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <toml++/toml.h>

namespace tomocam {
    struct ReconParams {
        std::string input_file_path;  // Path to the input data file
        std::string angles_file_path; // Path to the angles file
        std::string output_file_path; // Path to the output file
        float orientation = 0.0f;     // Orientation angle
        size_t thickness = 41;
        size_t maxIters = 50;  // Maximum number of iterations
        float sigma = 1000.0f; // Regularization parameter
        float tol = 1e-5f;     // Tolerance for convergence
        float xtol = 1e-5f;    // Tolerance for solution change

        static ReconParams from_toml(const std::string &input_deck) {
            ReconParams params;

            toml::table config;
            try {
                config = toml::parse_file(input_deck);
            } catch (const toml::parse_error &err) {
                throw std::runtime_error(
                    std::format("Could not open or parse input file {}, error: {}\n",
                                input_deck, err.description()));
            }

            // Read [input] section
            if (auto input = config["input"].as_table()) {
                params.input_file_path =
                    (*input)["filename"].value_or<std::string>("");
                params.angles_file_path =
                    (*input)["angles"].value_or<std::string>("");
                params.orientation = (*input)["orientation"].value_or<float>(0.0f);

                // Validate input paths
                if (params.input_file_path.empty()) {
                    throw std::runtime_error("Input filename is required");
                }
                if (!std::filesystem::exists(params.input_file_path)) {
                    throw std::runtime_error("Input file does not exist: " +
                                             params.input_file_path);
                }

                if (params.angles_file_path.empty()) {
                    throw std::runtime_error("Angles filename is required");
                }
                if (!std::filesystem::exists(params.angles_file_path)) {
                    throw std::runtime_error("Angles file does not exist: " +
                                             params.angles_file_path);
                }
            }

            // Read [output] section
            if (auto output = config["output"].as_table()) {
                params.output_file_path =
                    (*output)["filename"].value_or<std::string>("");

                // Validate output path
                if (params.output_file_path.empty()) {
                    throw std::runtime_error("Output filename is required");
                }
                std::filesystem::path outPath(params.output_file_path);
                std::filesystem::path outDir = outPath.parent_path();
                if (!outDir.empty() && !std::filesystem::exists(outDir)) {
                    throw std::runtime_error("Output directory does not exist: " +
                                             outDir.string());
                }
            }

            // Read [recon_params] section
            if (auto recon = config["recon_params"].as_table()) {
                params.maxIters = (*recon)["max_iters"].value_or<size_t>(50);
                params.thickness = (*recon)["thickness"].value_or<size_t>(41);
                params.sigma = (*recon)["sigma"].value_or<float>(1000.0f);
                params.tol = (*recon)["tol"].value_or<float>(1e-5f);
                params.xtol = (*recon)["xtol"].value_or<float>(1e-5f);
            }

            return params;
        }

        static void dump_config() {

            auto config = toml::table{
                {"input", toml::table{{"filename", "path/to/input/file.tiff"},
                                      {"angles", "path/to/angles/file.txt"},
                                      {"orientation", 0.0f}}},
                {"output", toml::table{{"filename", "path/to/output/directory"}}},
                {"recon_params", toml::table{{"max_iters", static_cast<int64_t>(50)},
                                             {"sigma", 1000.0f},
                                             {"thickness", static_cast<int64_t>(41)},
                                             {"tol", 1e-5f},
                                             {"xtol", 1e-5f}}}};
            std::ofstream os("config.toml");
            os << config << std::endl;
        }
        void print_toml(std::ostream &os) const {

            auto config = toml::table{
                {"input", toml::table{{"filename", input_file_path},
                                      {"angles", angles_file_path},
                                      {"orientation", orientation}}},
                {"output", toml::table{{"filename", output_file_path}}},
                {"recon_params",
                 toml::table{{"max_iters", static_cast<int64_t>(maxIters)},
                             {"sigma", sigma},
                             {"thickness", static_cast<int64_t>(thickness)},
                             {"tol", tol},
                             {"xtol", xtol}}}};
            os << config << std::endl;
        }
    };
} // namespace tomocam
#endif // RECON_PARAMS_H
