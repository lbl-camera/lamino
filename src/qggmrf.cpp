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

    // sigma_q  = pow(sigma, Q)
    // sigma_qp = pow(sigma, Q - p)
    // Both are precomputed once per qggmrf call and passed in.
    template <typename T>
    T d_potential_fcn(T delta, T sigma_q, T sigma_qp, T p) {

        T C = static_cast<T>(MRF_C);

        T abs_delta = std::abs(delta);
        T temp1 = std::pow(abs_delta, MRF_Q - p) / sigma_qp;
        T temp3 = C + temp1;

        if (abs_delta > 0.f) {
            // With Q=2: sign(delta)*pow(abs(delta), Q-1) == delta
            return (delta / (temp3 * sigma_q)) *
                   (MRF_Q - ((MRF_Q - p) * temp1) / temp3);
        } else {
            return MRF_Q / (sigma_q * MRF_C);
        }
    }

    template <typename T>
    void qggmrf(Array<T> &gradient, const Array<T> &x, T sigma, T p) {

        // dims
        auto dims = x.dims();

        // write out for clarity
        size_t nslices = dims.x();
        size_t nrows = dims.y();
        size_t ncols = dims.z();

        // Precompute sigma terms once; they are constant for the entire call.
        const T sigma_q  = std::pow(sigma, static_cast<T>(MRF_Q));
        const T sigma_qp = std::pow(sigma, static_cast<T>(MRF_Q) - p);

#pragma omp parallel for collapse(3)
        for (size_t i = 1; i < nslices - 1; i++) {
            for (size_t j = 1; j < nrows - 1; j++) {
                for (size_t k = 1; k < ncols - 1; k++) {
                    T sum_du = 0;
                    T u = x[{i, j, k}];
                    for (int dx = -1; dx < 2; dx++) {
                        for (int dy = -1; dy < 2; dy++) {
                            for (int dz = -1; dz < 2; dz++) {
                                // center voxel weight is 0; skip it
                                if (dx == 0 && dy == 0 && dz == 0) continue;
                                T v = x.at(i + dx, j + dy, k + dz);
                                T w = weight<T>(dx + 1, dy + 1, dz + 1);
                                sum_du += w * d_potential_fcn(v - u, sigma_q, sigma_qp, p);
                            }
                        }
                    }
                    gradient[{i, j, k}] += sum_du;
                }
            }
        }
    }
    // template instantiations
    template void qggmrf<float>(Array<float> &, const Array<float> &, float, float);
    template void qggmrf<double>(Array<double> &, const Array<double> &, double,
                                 double);
} // namespace tomocam::opt
