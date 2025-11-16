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
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "fft.h"

#ifndef TOMOCAM_PRECOND__H
    #define TOMOCAM_PRECOND__H

template <typename T>
class RampPreconditioner {
  private:
    Array<T> filter_;

  public:
    RampPreconditioner(dims_t dims) {

        dims_t filter_dims = {1, dims.n2, dims.n3};
        filter_ = Array<T>(filter_dims);

        // Create 1-d ramp filter in frequency domain
        std::vector<T> freq(dims.n3, 0);
        for (size_t i = 0; i < dims.n3; ++i) {
            T u = (i < dims.n3 / 2) ? i : i - dims.n3;
            freq[i] = std::abs(u);
        }

        T f_max = *std::max_element(freq.begin(), freq.end());
        T eps = ((T)1 / (T)dims.n3) * std::sqrt(2) * f_max;

        // create a 2D ramp filter by outer addition
        for (size_t j = 0; j < dims.n2; ++j) {
            for (size_t i = 0; i < dims.n3; ++i) {
                auto f = std::sqrt(freq[i] * freq[i] + freq[j] * freq[j]);
                filter_[{0, j, i}] = std::max(eps, f);
            }
        }
    }

    Array<T> operator()(const Array<T> &input) const {

        // Apply the ramp filter in frequency domain
        auto fft_input = fft2D(input);
        auto filtered = filter_ * fft_input;
        auto output = ifft2D(filtered);
        return array::to_real(output);
    }
};

#endif // TOMOCAM_PRECOND__H
