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
#ifndef TOMOCAM_PRECOND_H
#define TOMOCAM_PRECOND_H

#include <algorithm>
#include <execution>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "fft.h"

namespace tomocam::opt {
    template <typename T>
    class RampPreconditioner {
      private:
        Array<T> filter_;

      public:
        RampPreconditioner(dims_t dims) {

            dims_t filter_dims = {1, dims.n1, dims.n2 / 2 + 1};
            filter_ = Array<T>(filter_dims);

            // After transposing (n1,n2,n3) -> (n3,n1,n2), each 2D slice is (n1,n2).
            // r2c FFT of (n1,n2) gives (n1, n2/2+1), so x-freq spans n2 and y-freq spans n1.

            // Create 1-d ramp filter in x-direction (corresponds to n2 after transpose)
            int n2 = static_cast<int>(dims.n2);
            std::vector<T> xfreq(dims.n2 / 2 + 1, 0);
            for (int i = 0; i < n2 / 2 + 1; ++i) { xfreq[i] = (T)i / (T)n2; }

            // Create 1-d ramp filter in y-direction (corresponds to n1 after transpose)
            int n1 = static_cast<int>(dims.n1);
            std::vector<T> yfreq(dims.n1, 0);
            for (int i = 0; i < n1; ++i) {
                yfreq[i] = (i <= n1 / 2) ? (T)i / (T)n1 : (T)(i - n1) / (T)n1;
            }

            // create a 2D ramp filter by outer addition
            for (size_t j = 0; j < dims.n1; ++j) {
                for (size_t i = 0; i < dims.n2 / 2 + 1; ++i) {
                    auto f = std::sqrt(xfreq[i] * xfreq[i] + yfreq[j] * yfreq[j]);
                    filter_[{0, j, i}] = f;
                }
            }
#ifdef DEBUG
            for (auto f : filter_) {
                if (std::isnan(f) || std::isinf(f)) {
                    throw std::runtime_error(
                        "Error: Ramp filter contains NaN/Inf values.\n");
                }
            }
#endif
        }

        Array<T> apply(const Array<T> &input) const {

            // Apply the ramp filter in frequency domain
            auto dims = input.dims();
            auto scale = 1.0 / (dims.n1 * dims.n2);

            // transpose input to make y-z plane contiguous in memory
            auto inout_T = array::transpose(input, {2, 0, 1});
            auto tdims = inout_T.dims(); // (n3, n1, n2)

            auto fft_input = fft::fft2_r2c(inout_T);
            for (size_t i = 0; i < tdims.n1; ++i) {
                auto slice = fft_input.slice(i, i + 1);
                std::transform(std::execution::seq, slice.begin(), slice.end(),
                               filter_.begin(), slice.begin(),
                               std::multiplies<std::complex<T>>());
            }
            auto filtered = fft::fft2_c2r(fft_input, tdims);

            // transpose back to original layout and apply scaling
            return array::transpose(filtered * scale, {1, 2, 0});
        }
    };
} // namespace tomocam::opt
#endif // TOMOCAM_PRECOND_H
