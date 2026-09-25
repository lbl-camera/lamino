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

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <toml++/toml.h>
#include <tuple>
#include <vector>

#include "array.h"
#include "mask.h"
#include "recon_params.h"
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

    // Angles are taken to be degrees if any |angle| > 2*pi, otherwise radians.
    template <typename T>
    inline void to_radians_if_degrees(std::vector<T> &angles) {
        T max_abs = T(0);
        for (auto a : angles) max_abs = std::max(max_abs, std::abs(a));
        if (max_abs > T(2 * M_PI)) {
            for (auto &a : angles) { a = a * T(M_PI) / T(180); }
        }
    }

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
        if (angles.empty()) {
            throw std::runtime_error(
                std::format("No angles found in file: {}", filepath));
        }

        to_radians_if_degrees(angles);
        return angles;
    }

    // Read per-projection COR shifts from a two-column text file (dx dy per line).
    template <typename T>
    inline std::vector<std::array<T, 2>>
    read_shifts_file(const std::string &filepath) {
        std::ifstream fp(filepath);
        if (!fp.is_open()) {
            throw std::runtime_error(
                std::format("Could not open shifts file: {}", filepath));
        }
        std::vector<std::array<T, 2>> shifts;
        T dx, dy;
        while (fp >> dx >> dy) { shifts.push_back({dx, dy}); }
        if (shifts.empty()) {
            throw std::runtime_error(
                std::format("No shifts found in file: {}", filepath));
        }
        return shifts;
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

            if (!input_table->contains("filename") ||
                !input_table->contains("angles") ||
                !input_table->contains("gamma")) {
                throw std::runtime_error("[[input]] entry must have 'filename', "
                                         "'angles', and 'gamma' fields");
            }
            auto filename = (*input_table)["filename"].value<std::string>();
            if (!filename.has_value()) {
                throw std::runtime_error("[[input]] 'filename' must be a string");
            }
            if (!std::filesystem::exists(*filename)) {
                throw std::runtime_error(
                    std::format("Projection file does not exist: {}", *filename));
            }
            auto angles_file = (*input_table)["angles"].value<std::string>();
            if (!angles_file.has_value()) {
                throw std::runtime_error("[[input]] 'angles' must be a string");
            }
            if (!std::filesystem::exists(*angles_file)) {
                throw std::runtime_error(
                    std::format("Angles file does not exist: {}", *angles_file));
            }
            auto gamma = (*input_table)["gamma"].value<T>();
            if (!gamma.has_value()) {
                throw std::runtime_error("[[input]] 'gamma' field must be a number");
            }

            // gamma: stored as positive radians (rotation of beam azimuth)
            T gamma_rad = *gamma * T(M_PI) / T(180);
            // beta: out-of-plane tilt, defaults to 0
            T beta_rad = (*input_table)["beta"].value_or<T>(0) * T(M_PI) / T(180);

            // COR shifts: 'shifts' file takes priority over 'cor-offset' scalar.
            std::vector<std::array<T, 2>> per_proj_shifts;
            if (input_table->contains("shifts")) {
                auto shifts_path = (*input_table)["shifts"].value<std::string>();
                if (!shifts_path.has_value())
                    throw std::runtime_error(
                        "[[input]] 'shifts' must be a string path");
                if (!std::filesystem::exists(*shifts_path))
                    throw std::runtime_error(
                        std::format("Shifts file does not exist: {}", *shifts_path));
                per_proj_shifts = read_shifts_file<T>(*shifts_path);
            } else if (input_table->contains("cor-offset")) {
                auto offsets_array = (*input_table)["cor-offset"].as_array();
                if (!offsets_array || offsets_array->size() != 2) {
                    throw std::runtime_error(
                        "[[input]] 'cor-offset' must be in form [dx, dy]");
                }
                std::array<T, 2> buf;
                for (size_t i = 0; i < 2; ++i) {
                    auto val = (*offsets_array)[i].value<T>();
                    if (!val.has_value())
                        throw std::runtime_error(
                            "[[input]] 'cor-offset' must be an array of numbers");
                    buf[i] = *val;
                }
                per_proj_shifts = {buf}; // sentinel: one element = broadcast
            }

            auto projs = tomocam::tiff::read(*filename);
            projs = tomocam::mask_infs_nans(projs);
            auto angles = read_angles_file<T>(*angles_file);

            // Broadcast scalar cor-offset to all projections now that N is known.
            if (per_proj_shifts.size() == 1) {
                per_proj_shifts.assign(angles.size(), per_proj_shifts[0]);
            } else if (per_proj_shifts.empty()) {
                per_proj_shifts.assign(angles.size(), std::array<T, 2>{T(0), T(0)});
            } else if (per_proj_shifts.size() != angles.size()) {
                throw std::runtime_error(std::format(
                    "shifts file has {} entries but {} projections were loaded",
                    per_proj_shifts.size(), angles.size()));
            }

            // Override angles, gamma, beta, and shifts from external alignment TOML.
            if (input_table->contains("alignment")) {
                auto align_path = (*input_table)["alignment"].value<std::string>();
                if (!align_path.has_value())
                    throw std::runtime_error(
                        "[[input]] 'alignment' must be a string path");
                if (!std::filesystem::exists(*align_path))
                    throw std::runtime_error(std::format(
                        "Alignment file does not exist: {}", *align_path));
                auto align_tbl = read_toml_file(*align_path);

                if (auto gv = align_tbl["gamma_deg"].value<T>())
                    gamma_rad = *gv * T(M_PI) / T(180);
                if (auto bv = align_tbl["beta_deg"].value<T>())
                    beta_rad = *bv * T(M_PI) / T(180);

                auto *ang_arr = align_tbl["corrected_angles_deg"].as_array();
                if (!ang_arr)
                    throw std::runtime_error(
                        "Alignment TOML missing 'corrected_angles_deg'");
                angles.clear();
                for (auto &a : *ang_arr)
                    angles.push_back(a.value<T>().value() * T(M_PI) / T(180));

                auto *spx = align_tbl["shifts_px"].as_array();
                if (spx) {
                    per_proj_shifts.clear();
                    for (auto &e : *spx) {
                        auto *pair = e.as_array();
                        if (!pair || pair->size() < 2)
                            throw std::runtime_error("Alignment TOML 'shifts_px' "
                                                     "entries must be [dx, dy]");
                        T dx = (*pair)[0].value<T>().value();
                        T dy = (*pair)[1].value<T>().value();
                        per_proj_shifts.push_back({dx, dy});
                    }
                    if (per_proj_shifts.size() != angles.size())
                        throw std::runtime_error(std::format(
                            "Alignment TOML: shifts_px has {} entries but "
                            "corrected_angles_deg has {}",
                            per_proj_shifts.size(), angles.size()));
                }
            }

            datasets.push_back({std::move(projs), std::move(angles), gamma_rad,
                                beta_rad, std::move(per_proj_shifts)});
        }
        return datasets;
    }

    inline ReconParams parse_recon_params(const toml::table &config) {

        ReconParams p;
        // Read [recon_params] section
        auto recon = config["recon_params"];
        if (!recon) {
            throw std::runtime_error(
                "Missing [recon_params] section in config file");
        }
        // read max_outer_iters
        p.maxIters = recon["max_iters"].value_or<size_t>(100);

        // read recon_dims
        std::array<size_t, 3> recon_dims;
        const auto *dims = recon["recon_dims"].as_array();
        if (dims && dims->size() == 3) {
            for (size_t i = 0; i < 3; ++i) {
                size_t temp = (*dims)[i].value_or<size_t>(0);
                if (temp == 0) {
                    throw std::runtime_error(
                        std::format("[recon_params] 'recon_dims[{}]' must be a "
                                    "positive integer",
                                    i));
                }
                if (temp % 2 == 0) {
                    temp -= 1; // make sure it's odd
                }
                recon_dims[i] = temp;
            }
        } else {
            throw std::runtime_error("[recon_params] 'recon_dims' must be an "
                                     "array of three integers");
        }
        p.recon_dims = recon_dims;
        // if reocn dims[0]/dims[2] > 0.15, remind user that since this is
        // laminography, the thickness is expected to be much smaller than the
        // in-plane dimensions, and this may lead to increased memory usage
        float ratio =
            static_cast<float>(recon_dims[0]) / static_cast<float>(recon_dims[2]);
        if (ratio > 0.15f) {
            std::cerr << std::format(
                "Warning: recon_dims[0] / recon_dims[2] = {:.2f} > 0.15. "
                "In laminography, the thickness (recon_dims[0]) "
                "is expected to be much smaller than the in-plane "
                "dimensions (recon_dims[1], recon_dims[2]). "
                "This may lead to increased memory usage.\n",
                ratio);
        }
        p.tol = recon["tol"].value_or<float>(1e-5f);
        p.xtol = recon["xtol"].value_or<float>(1e-5f);
        // read regularizer type
        auto reg = recon["regularizer"].as_table();
        if (!reg) {
            p.regularizer = Regularizer::UNCONSTRAINED;
        } else {
            auto reg_str = (*reg)["method"].value_or<std::string>("split_bregman");
            if (reg_str == "split_bregman") {
                p.regularizer = Regularizer::SPLIT_BREGMAN;
                auto params = (*reg)["split_bregman"].as_table();
                if (!params) {
                    throw std::runtime_error(
                        "Missing [recon_params.regularizer.split_bregman] section "
                        "in config file");
                }
                p.innerIters = (*params)["inner_iters"].value_or<size_t>(1);
                p.lambda = (*params)["lambda"].value_or<float>(0.1f);
                p.mu = (*params)["mu"].value_or<float>(10.0f);
            } else {
                throw std::runtime_error("[recon_params] 'regularizer' must be "
                                         "either 'qGGMRF' or 'split_bregman'");
            }
        }
        return p;
    }

    // Output parameters
    struct OutputParams {
        std::string filepath;
        std::vector<std::string> formats;

        OutputParams() = default;
        OutputParams(const toml::table &config) {

            // Read [output] section
            auto output = config["output"];
            if (!output) {
                throw std::runtime_error("Missing [output] section in config file");
            }
            filepath = output["filename"].value_or<std::string>("./recon.tiff");

            // Read formats array
            auto formats_array = output["formats"].as_array();
            if (formats_array) {
                for (auto &elem : *formats_array) {
                    auto fmt = elem.value<std::string>();
                    if (!fmt.has_value()) {
                        throw std::runtime_error(
                            "[output] 'formats' array must contain strings");
                    }
                    if (*fmt != "tiff" && *fmt != "vti") {
                        throw std::runtime_error(
                            std::format("[output] invalid format '{}'. Must be "
                                        "'tiff' or 'vti'",
                                        *fmt));
                    }
                    formats.push_back(*fmt);
                }
            } else {
                // Default to both formats if not specified
                formats = {"tiff", "vti"};
            }

            // Remove duplicates while preserving order
            std::vector<std::string> unique_formats;
            for (const auto &fmt : formats) {
                if (std::find(unique_formats.begin(), unique_formats.end(), fmt) ==
                    unique_formats.end()) {
                    unique_formats.push_back(fmt);
                }
            }
            formats = std::move(unique_formats);
        }

        bool has_format(const std::string &fmt) const {
            return std::find(formats.begin(), formats.end(), fmt) != formats.end();
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
        outfile << "beta = 0        # optional: out-of-plane tilt in degrees\n";
        outfile << "# cor-offset = [0.0, 0.0]  # optional: [dx, dy] COR shift in "
                   "pixels\n";
        outfile << "# shifts = \"/path/to/gamma0_shifts.txt\"  # optional: "
                   "per-projection shifts\n";
        outfile << "# alignment = \"/path/to/alignment.toml\"  # optional: override "
                   "angles/shifts\n";
        outfile << "\n";
        outfile << "[[input]]\n";
        outfile << "filename = \"/path/to/gamma45_stack.tiff\"\n";
        outfile << "angles = \"/path/to/gamma45_angles.txt\"\n";
        outfile << "gamma = 45\n";
        outfile << "beta = 0\n";
        outfile << "\n";
        outfile << "[output]\n";
        outfile << "filename = \"output.tiff\"\n";
        outfile << "formats = [\"tiff\", \"vti\"]  # available: \"tiff\", \"vti\"\n";
        outfile << "\n";
        outfile << "[recon_params]\n";
        outfile << "max_outer_iters = 50\n";
        outfile << "tol = 1e-5\n";
        outfile << "xtol = 1e-5\n";
        outfile << "recon_dims = [51, 511, 511]\n";
        outfile << "\n";
        outfile << "[recon_params.regularizer]\n";
        outfile << "method = \"split_bregman\"\n";
        outfile << "\n";
        outfile << "[recon_params.regularizer.split_bregman]\n";
        outfile << "lambda = 0.1\n";
        outfile << "mu = 10.0\n";
        outfile << "\n";
        outfile << "# alternatively, for qGGMRF regularization:\n";
        outfile << "# [recon_params.regularizer]\n";
        outfile << "# method = \"qGGMRF\"\n";
        outfile << "# [recon_params.regularizer.qGGMRF]\n";
        outfile << "# sigma = 1000.0\n";
        outfile << "# p = 1.2\n";
        outfile.close();
    }
} // namespace tomocam
#endif // CONFIG_H
