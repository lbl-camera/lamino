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

#include <cuda/std/complex>
#include <thrust/device_vector.h>

#include "gpu/device_ptr.h"
#include "gpu/phase_shift.h"
#include "gpu/utils.h"

constexpr double TOMOCAM_PHASE_SHIFT_PI = 3.14159265358979323846;

namespace tomocam::gpu {

    template <typename T>
    __global__ void phase_shift2d_kernel(DevicePtr<Complex<T>> data, const T *dx,
                                         const T *dy) {
        using cuda::std::cos;
        using cuda::std::sin;

        auto dims = data.dims();
        auto idx = Index3D();
        if (idx < dims) {
            T dqx = T(2) * T(TOMOCAM_PHASE_SHIFT_PI) / static_cast<T>(dims.n3);
            T dqy = T(2) * T(TOMOCAM_PHASE_SHIFT_PI) / static_cast<T>(dims.n2);
            T qx = (idx.z + T(0.5)) * dqx - T(TOMOCAM_PHASE_SHIFT_PI);
            T qy = (idx.y + T(0.5)) * dqy - T(TOMOCAM_PHASE_SHIFT_PI);
            T phase = -(qx * dx[idx.x] + qy * dy[idx.x]);
            data[idx] *= Complex<T>(cos(phase), sin(phase));
        }
    }

    template <typename T>
    void phase_shift2d(DeviceArray<Complex<T>> &data,
                       const thrust::device_vector<T> &dx,
                       const thrust::device_vector<T> &dy) {
        dim3 threads(1, 16, 16);
        auto g = make_grid(data.dims(), threads);
        phase_shift2d_kernel<T><<<g, threads>>>(
            data, thrust::raw_pointer_cast(dx.data()),
            thrust::raw_pointer_cast(dy.data()));
        SAFE_CALL(cudaGetLastError());
    }

    template void phase_shift2d<float>(DeviceArray<Complex<float>> &,
                                       const thrust::device_vector<float> &,
                                       const thrust::device_vector<float> &);
    template void phase_shift2d<double>(DeviceArray<Complex<double>> &,
                                        const thrust::device_vector<double> &,
                                        const thrust::device_vector<double> &);

} // namespace tomocam::gpu
