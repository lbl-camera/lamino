
#include <complex>
#include <string>

#include "array.h"

#ifndef FILTERS__H
    #define FILTERS__H

namespace tomocam {
    template <typename T>
    void apply_filter(Array<std::complex<T>> &arr,
        const std::string &filter_name) {

        // allocate return values
        auto out = Array<std::complex<T>>(arr.dims());

        // apply filter
        T rmax = static_cast<T>(arr.nrows()) / 2;
        for (uint64_t i = 0; i < arr.nslices(); ++i) {
            for (uint64_t j = 0; j < arr.nrows(); ++j) {

                T filt = 0;
                T f = std::abs(static_cast<T>(j) - rmax) / rmax;
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
                } else {
                    throw std::invalid_argument(
                        "Unknown filter name: " + filter_name);
                }

                for (uint64_t k = 0; k < arr.ncols(); ++k) {
                    out[{i, j, k}] = filt * arr[{i, j, k}];
                }
            }
        }
        arr = std::move(out);
    }
} // namespace tomocam
#endif
