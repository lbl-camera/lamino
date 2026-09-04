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

#include <cuda_runtime.h>
#include <thrust/device_vector.h>

#include "gpu/device_array.h"
#include "gpu/device_ptr.h"
#include "gpu/polar_grid.h"
#include "gpu/utils.h"

constexpr double PI = 3.14159265358979323846;

namespace tomocam::gpu {

    // Literal port of tomocam::RotationTranspose(theta,gamma,beta) applied to
    // (qX, qY, 0) -- see include/rotation.h. Returns transpose(Rx(theta) *
    // Ry(beta) * Rz(gamma)) * (qX, qY, 0).
    template <typename T>
    __device__ __forceinline__ void rotate_transpose(T theta, T gamma, T beta,
                                                      T qX, T qY, T &qx, T &qy,
                                                      T &qz) {
        T ct = cos(theta), st = sin(theta);
        T cb = cos(beta), sb = sin(beta);
        T cg = cos(gamma), sg = sin(gamma);

        qx = qX * (cb * cg) + qY * (ct * sg + st * sb * cg);
        qy = qX * (-cb * sg) + qY * (ct * cg - st * sb * sg);
        qz = qX * sb + qY * (-st * cb);
    }

    template <typename T>
    __global__ void make_polar_grid_kernel(T *theta, T *gamma, T *beta,
                                           DevicePtr<T> x, DevicePtr<T> y,
                                           DevicePtr<T> z, DevicePtr<T> w) {

        auto dims = x.dims();
        auto idx = Index3D();

        T dX = 2 * PI / (T)dims.n3;
        T dY = 2 * PI / (T)dims.n2;
        if (idx < dims) {
            // Offset 0.5 gives qX[N-1-k] = -qX[k] for ODD N, ensuring
            // Hermitian-symmetric NUFFT output for real inputs.
            // Even N is not supported: FINUFFT's asymmetric mode range
            // [-N/2, N/2-1] breaks Hermitian symmetry at the Nyquist.
            T qX = (idx.z + 0.5) * dX - PI;
            T qY = (idx.y + 0.5) * dY - PI;

            T qx, qy, qz;
            rotate_transpose(theta[idx.x], gamma[idx.x], beta[idx.x], qX, qY, qx,
                             qy, qz);
            x[idx] = qx;
            y[idx] = qy;
            z[idx] = qz;
            w[idx] = (fabs(qx) <= PI && fabs(qy) <= PI && fabs(qz) <= PI) ? T(1)
                                                                          : T(0);
        }
    }

    template <typename T>
    void make_polar_grid(thrust::device_vector<T> &d_angles,
                         thrust::device_vector<T> &d_gammas,
                         thrust::device_vector<T> &d_betas, DeviceArray<T> &x,
                         DeviceArray<T> &y, DeviceArray<T> &z, DeviceArray<T> &w) {

        auto dims = x.dims();
        dim3 blockSize(1, 16, 16);
        dim3 gridSize;
        gridSize.x = (dims.n1 + blockSize.x - 1) / blockSize.x;
        gridSize.y = (dims.n2 + blockSize.y - 1) / blockSize.y;
        gridSize.z = (dims.n3 + blockSize.z - 1) / blockSize.z;
        make_polar_grid_kernel<T><<<gridSize, blockSize>>>(
            thrust::raw_pointer_cast(d_angles.data()),
            thrust::raw_pointer_cast(d_gammas.data()),
            thrust::raw_pointer_cast(d_betas.data()), x, y, z, w);
        SAFE_CALL(cudaGetLastError());
    }

    template <typename T>
    PolarGrid<T>::PolarGrid(const std::vector<T> &theta_, size_t nrows,
                            size_t ncols, T gamma, T beta) {
        theta = theta_;
        gammas = std::vector<T>(theta.size(), gamma);
        betas = std::vector<T>(theta.size(), beta);

        auto dims = dims_t{theta.size(), nrows, ncols};
        npts = dims.n1 * dims.n2 * dims.n3;
        x = DeviceArray<T>(dims);
        y = DeviceArray<T>(dims);
        z = DeviceArray<T>(dims);
        w = DeviceArray<T>(dims);
        angles = thrust::device_vector<T>(theta);
        d_gammas = thrust::device_vector<T>(gammas);
        d_betas = thrust::device_vector<T>(betas);
        make_polar_grid(angles, d_gammas, d_betas, x, y, z, w);
    }

    template struct PolarGrid<float>;
    template struct PolarGrid<double>;

} // namespace tomocam::gpu
