/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 *National Laboratory (subject to receipt of any required approvals from the
 *U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 *IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 *the U.S. Government has been granted for itself and others acting on its
 *behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 *to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#include "array.h"

namespace tomocam::opt {

    constexpr float WEIGHT[27] = {0.0302f, 0.037f, 0.0302f, 0.037f, 0.0523f, 0.037f,
                                  0.0302f, 0.037f, 0.0302f, 0.037f, 0.0523f, 0.037f,
                                  0.0532f, 0.f,    0.0523f, 0.037f, 0.0523f, 0.037f,
                                  0.0302f, 0.037f, 0.0302f, 0.037f, 0.0523f, 0.037f,
                                  0.0302f, 0.037f, 0.0302f};

    template <typename T>
    T weight(size_t x, size_t y, size_t z) {
        return static_cast<T>(WEIGHT[x * 9 + y * 3 + z]);
    }
    const double MRF_Q = 2;
    const double MRF_C = 0.001;

    template <typename T>
    T sign(T x) {
        if (x < 0.f)
            return -1;
        else
            return 1;
    }

    template <typename T>
    T d_potential_fcn(T delta, T sigma, T p) {

        T Q = static_cast<T>(MRF_Q);
        T C = static_cast<T>(MRF_C);

        T sigma_q = std::pow(sigma, Q);
        T sigma_q_p = std::pow(sigma, Q - p);

        T temp1 = std::pow(std::abs(delta), Q - p) / sigma_q_p;
        T temp2 = std::pow(std::abs(delta), Q - 1);
        T temp3 = C + temp1;

        if (std::abs(delta) > 0.f) {
            return ((sign(delta) * temp2 / (temp3 * sigma_q)) *
                    (MRF_Q - ((MRF_Q - p) * temp1) / (temp3)));
        } else {
            return MRF_Q / (sigma_q * MRF_C);
        }
    }

    template <typename T>
    void qggmrf(const Array<T> &input, Array<T> &output, T sigma, T p) {

        // dims
        auto dims = input.dims();

        // write out for clarity
        size_t nslices = dims.x();
        size_t nrows = dims.y();
        size_t ncols = dims.z();

#pragma omp parallel for collapse(3)
        for (size_t i = 1; i < nslices - 1; i++) {
            for (size_t j = 1; j < nrows - 1; j++) {
                for (size_t k = 1; k < ncols - 1; k++) {
                    T sum_du = 0;
                    T u = input[{i, j, k}];
                    for (int x = -1; x < 2; x++) {
                        for (int y = -1; y < 2; y++) {
                            for (int z = -1; z < 2; z++) {
                                T v = input.at(i + x, j + y, k + z);
                                T w = weight<T>(x + 1, y + 1, z + 1);
                                sum_du += w * d_potential_fcn(v - u, sigma, p);
                            }
                        }
                    }
                    output[{i, j, k}] = sum_du;
                }
            }
        }
    }
    // template instantiations
    template void qggmrf<float>(const Array<float> &, Array<float> &, float, float);
    template void qggmrf<double>(const Array<double> &, Array<double> &, double,
                                 double);
} // namespace tomocam::opt
