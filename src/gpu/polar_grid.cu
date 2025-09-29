#include <__clang_cuda_runtime_wrapper.h>
#include <cmath>
#include <cuda_runtime.h>

#include "device_ptr.cuh"
#include "utils.cuh"

namespace tomocam::gpu {
    __global__ void make_polar_grid(DevicePtr x, DevicePtr y, DevicePtr z,
        float *angles) {

        auto dims = x.dims();
        uint3 idx = Index3D();
        if (idx < dims) {
            float cos_t = cos(angles[idx.x()]);
            float sin_t = sin(angles[idx.x()]);

            float dx = (2 * M_PI) / static_cast<float>(dims.x());
            float dr = (2 * M_PI) / static_cast<float>(dims.y());

            x[idx] = dims.z() * dx - M_PI;
            y[idx] = (dims.y() * dr - M_PI) * sin_t;
            z[idx] = (dims.y() * dr - M_PI) * cos_t;
        }
    }

    __global__ void rotate_polar_grid(DevicePtr x, DevicePtr y, float gamma) {

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
