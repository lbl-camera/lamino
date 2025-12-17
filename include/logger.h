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

#ifndef LOGGER__H
#define LOGGER__H

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

namespace tomocam {

    enum class LogMode { SILENT, STDOUT, FILE };

    class Logger {
      private:
        LogMode mode_;
        std::unique_ptr<std::ofstream> file_stream_;

      public:
        Logger(LogMode mode = LogMode::STDOUT, const std::string &filename = "")
            : mode_(mode) {
            if (mode_ == LogMode::FILE) {
                file_stream_ = std::make_unique<std::ofstream>(filename);
                if (!file_stream_->is_open()) {
                    std::cerr << "Warning: Unable to open log file '" << filename
                              << "'. Falling back to STDOUT.\n";
                    mode_ = LogMode::STDOUT;
                }
            }
        }

        template <typename... Args>
        void log(const std::string &message) {
            if (mode_ == LogMode::SILENT) return;

            if (mode_ == LogMode::STDOUT) {
                std::cout << message;
            } else if (mode_ == LogMode::FILE && file_stream_) {
                *file_stream_ << message;
                file_stream_->flush();
            }
        }

        LogMode get_mode() const { return mode_; }

        ~Logger() {
            if (file_stream_ && file_stream_->is_open()) { file_stream_->close(); }
        }
    };

} // namespace tomocam

#endif // LOGGER__H
