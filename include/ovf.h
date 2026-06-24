#pragma once

#include "array.h"
#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace tomocam::ovf {

    template <typename T>
    std::array<Array<T>, 3> read(const std::string &filename) {

        std::ifstream file(filename, std::ios::binary);
        if (!file) { throw std::runtime_error("Cannot open OVF file: " + filename); }

        dims_t dims = {0, 0, 0};
        std::string data_format;
        std::string line;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] != '#') { break; }

            if (line.size() > 2 && line[1] == ' ') { line = line.substr(2); }

            size_t colon_pos = line.find(':');
            if (colon_pos == std::string::npos) continue;

            std::string key = line.substr(0, colon_pos);
            std::string value = line.substr(colon_pos + 1);

            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            if (key == "xnodes") {
                dims.n3 = std::stoul(value);
            } else if (key == "ynodes") {
                dims.n2 = std::stoul(value);
            } else if (key == "znodes") {
                dims.n1 = std::stoul(value);
            } else if (key == "Begin" && value.find("Data") != std::string::npos) {
                size_t space = value.find(' ');
                if (space != std::string::npos) {
                    data_format = value.substr(space + 1);
                    data_format.erase(data_format.find_last_not_of(" \t\n\r") + 1);
                }
                break;
            }
        }

        if (dims.n1 == 0 || dims.n2 == 0 || dims.n3 == 0) {
            throw std::runtime_error("Invalid dimensions in OVF file");
        }

        Array<T> mx(dims);
        Array<T> my(dims);
        Array<T> mz(dims);
        size_t ntotal = dims.n1 * dims.n2 * dims.n3;

        if (data_format == "Text") {
            for (size_t i = 0; i < ntotal; ++i) {
                while (std::getline(file, line)) {
                    if (!line.empty() && line[0] != '#') { break; }
                }
                std::istringstream iss(line);
                T x, y, z;
                if (!(iss >> x >> y >> z)) {
                    throw std::runtime_error(
                        "Failed to parse text data from OVF file");
                }
                // transpose the data to match (nz, ny, nx) order
                mx[i] = x;
                my[i] = y;
                mz[i] = z;
            }
        } else {
            throw std::runtime_error("Unsupported OVF data format: " + data_format);
        }

        return {std::move(mx), std::move(my), std::move(mz)};
    }

} // namespace tomocam::ovf
