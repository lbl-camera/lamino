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
#ifndef TOMOCAM_PRECOND__H
#define TOMOCAM_PRECOND__H

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

            // creae 2D ramp filter in frequency domain of size (1, n2, n3+1)

            dims_t filter_dims = {1, dims.n2, dims.n3 / 2 + 1};
            filter_ = Array<T>(filter_dims);

            // Create 1-d ramp filter along the x-axis (n3 direction)
            std::vector<T> freqx(dims.n3 / 2 + 1, 0);
            for (size_t i = 0; i < freqx.size(); ++i) {
                T u = (T)i / (T)dims.n3;
                freqx[i] = std::abs(u);
            }
            // Create 1-d ramp filter along the y-axis (n2 direction)
            std::vector<T> freqy(dims.n2, 0);
            for (size_t i = 0; i < freqy.size(); ++i) {
                T v = i < dims.n2 / 2 + 1 ? (T)i / (T)dims.n2
                                          : (T)(i - dims.n2) / (T)dims.n2;
                freqy[i] = std::abs(v);
            }

            T fx_max = *std::max_element(freqx.begin(), freqx.end());
            T fy_max = *std::max_element(freqy.begin(), freqy.end());
            T f_max = std::sqrt(fx_max * fx_max + fy_max * fy_max);
            T eps = 1.e-6 * f_max; // small value to avoid zero filter response

            // create a 2D ramp filter by outer addition
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t i = 0; i < dims.n3 / 2 + 1; ++i) {
                    auto f = std::sqrt(freqy[j] * freqy[j] + freqx[i] * freqx[i]);
                    filter_[{0, j, i}] = std::max(eps, f);
                }
            }
        }

        Array<T> apply(const Array<T> &input) const {

            // Apply the ramp filter in frequency domain
            auto dims = input.dims();
            auto fft_input = fft::fft2_r2c<T>(input);
            auto filtered = Array<std::complex<T>>(fft_input.dims());
            for (size_t i = 0; i < dims.n1; ++i) {
                for (size_t j = 0; j < dims.n2; ++j) {
                    for (size_t k = 0; k < dims.n3 / 2 + 1; ++k) {
                        filtered[{i, j, k}] =
                            fft_input[{i, j, k}] * filter_[{0, j, k}];
                    }
                }
            }
            auto output = fft::fft2_c2r<T>(filtered, input.dims());
            return output / (dims.n2 * dims.n3); // normalize by the number of pixels
        }
    };
} // namespace tomocam::opt
#endif // TOMOCAM_PRECOND__H
