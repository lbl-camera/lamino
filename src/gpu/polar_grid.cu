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
#include <cmath>
#include <cuda_runtime.h>

#include "gpu/utils.h"
#include "gpu/device_ptr.h"
#include "gpu/device_array.h"
#include "gpu/polar_grid_kernel.h"

namespace tomocam::gpu {

    template <typename T>
    __global__ void polar_grid_kenrel(DevicePtr<T> x, DevicePtr<T> y,
        DevicePtr<T> z, T *angles, T gamma) {

        auto dims = x.dims();
        uint3 idx = Index3D();
        if (idx < dims) {
            T cos_t = cos(angles[idx.x]);
            T sin_t = sin(angles[idx.x]);
            T cos_g = cos(gamma);
            T sin_g = sin(gamma);

            T dx = (2 * M_PI) / static_cast<T>(dims.n3);
            T dr = (2 * M_PI) / static_cast<T>(dims.n2);

            auto xcrd = idx.z * dx;
            auto ycrd = (idx.y * dr - M_PI) * cos_t;
            auto zcrd = (idx.y * dr - M_PI) * sin_t;
            // gama rotation
            x[idx] = xcrd * cos_g - ycrd * sin_g;
            y[idx] = xcrd * sin_g + ycrd * cos_g;
            z[idx] = zcrd;
        }
    }

    template <typename T>
    void calc_polar_grid(DeviceArray<T> &x, DeviceArray<T> &y,
        DeviceArray<T> &z, std::vector<T> &angles, T gamma) {

        dim3 block(1, 32, 32);
        dim3 grid = make_grid(x.dims(), block);
        polar_grid_kenrel<T><<<grid, block>>>(x, y, z, angles.data(), gamma);
        SAFE_CALL(cudaGetLastError());
    }
    // explicit template instantiation
    template void calc_polar_grid<float>(DeviceArray<float> &x, DeviceArray<float> &y,
        DeviceArray<float> &z, std::vector<float> &angles, float gamma);
    template void calc_polar_grid<double>(DeviceArray<double> &x, DeviceArray<double> &y,
        DeviceArray<double> &z, std::vector<double> &angles, double gamma);

} // namespace tomocam::gpu
