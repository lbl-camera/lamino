#include <cmath>
#include <cuda_runtime.h>

#include "gpu/device_ptr.cuh"
#include "gpu/utils.cuh"

namespace tomocam::gpu {

    template <typename T>
    __global__ void make_polar_grid(DevicePtr<T> x, DevicePtr<T> y,
        DevicePtr<T> z, float *angles) {

        auto dims = x.dims();
        uint3 idx = Index3D();
        if (idx < dims) {
            T cos_t = cos(angles[idx.x]);
            T sin_t = sin(angles[idx.x]);

            T dx = (2 * M_PI) / static_cast<T>(dims.x());
            T dr = (2 * M_PI) / static_cast<T>(dims.y());

            x[idx] = dims.z() * dx - M_PI;
            y[idx] = (dims.y() * dr - M_PI) * sin_t;
            z[idx] = (dims.y() * dr - M_PI) * cos_t;
        }
    }

    template <typename T>
    __global__ void rotate_polar_grid(DevicePtr<T> x, DevicePtr<T> y, T gamma) {

        auto dims = x.dims();
        uint3 idx = Index3D();
        if (idx < dims) {
            auto cos_g = cos(gamma);
            auto sin_g = sin(gamma);
            auto rx = x[idx];
            auto ry = y[idx];

            x[idx] = rx * cos_g - ry * sin_g;
            y[idx] = rx * sin_g + ry * cos_g;
        }
    }

} // namespace tomocam::gpu
