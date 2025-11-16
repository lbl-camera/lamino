// clang-format off
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
 //clang-format on

#include <complex>
#include <string>

#include "array.h"

#ifndef FILTERS__H
    #define FILTERS__H

namespace tomocam {
    template <typename T>
    void apply_filter(Array<std::complex<T>> &arr, const std::string &filter_name) {

        // apply filter
        T rmax = static_cast<T>(arr.ncols()) / 2;
        for (size_t i = 0; i < arr.nslices(); ++i) {
            for (size_t j = 0; j < arr.nrows(); ++j) {
                for (size_t k = 0; k < arr.ncols(); ++k) {

                    T filt = 0;
                    T f = std::abs(static_cast<T>(k) - rmax) / rmax;
                    if (filter_name == "ram-lak") {
                        if (f < 0.75) {
                            filt = f;
                        } else {
                            filt = 0;
                        }
                    } else if (filter_name == "shepp-logan") {
                        if (f == 0) {
                            filt = 1;
                        } else {
                            filt = f * std::sin(M_PI * f) / (M_PI * f);
                        }
                    } else if (filter_name == "cosine") {
                        filt = f * std::cos(M_PI * f / 2);
                    } else if (filter_name == "hamming") {
                        filt = f * (0.54 + 0.46 * std::cos(M_PI * f));
                    } else if (filter_name == "ramp") {
                        filt = f;
                    } else {
                        throw std::invalid_argument("Unknown filter name: " +
                                                    filter_name);
                    }

                    arr[{i, j, k}] *= filt;
                }
            }
        }
    }
} // namespace tomocam
#endif
