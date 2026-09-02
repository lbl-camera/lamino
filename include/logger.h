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

#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

namespace tomocam {

    enum class LogMode { SILENT, STDOUT, LOGFILE, BOTH };

    class Logger {
      private:
        uint8_t mode_ = 0x0;
        std::unique_ptr<std::ofstream> file_stream_ = nullptr;

      public:
        Logger(LogMode mode = LogMode::STDOUT, const std::string &filename = "") {

            if (mode == LogMode::SILENT) return;
            if (mode == LogMode::STDOUT || mode == LogMode::BOTH) { mode_ |= 0x1; }
            if (mode == LogMode::LOGFILE || mode == LogMode::BOTH) { mode_ |= 0x2; }

            if (mode_ & 0x2) {
                if (filename.empty()) {
                    auto now = std::chrono::floor<std::chrono::seconds>(
                        std::chrono::system_clock::now());
                    std::string ts = std::format("{:%Y%m%d_%H%M%S}", now);
                    std::string fname = "log_" + ts + ".log";
                    file_stream_ = std::make_unique<std::ofstream>(fname);
                } else {
                    file_stream_ = std::make_unique<std::ofstream>(filename);
                }
                if (!file_stream_->is_open()) {
                    std::cerr
                        << "Warning: Unable to open log file. Falling back to STDOUT.\n";
                    mode_ = 0x1;
                    file_stream_.reset();
                }
            }
        }

        void log(const std::string &message) {
            if (!mode_ || message.empty()) return;

            bool has_newline = message.back() == '\n';
            if (mode_ & 0x1) {
                std::cout << message << (has_newline ? "" : "\n");
            }
            if (mode_ & 0x2)
                *file_stream_ << message << (has_newline ? "" : "\n");
        }
    };

} // namespace tomocam

#endif // LOGGER_H
