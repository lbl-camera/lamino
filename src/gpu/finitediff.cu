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

#include "gpu/device_array.h"
#include "gpu/device_ptr.h"
#include "gpu/utils.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu::opt {

    // Forward differences with a zero (Neumann) difference at the last index of
    // each axis:
    //   (D u)[i] = u[i + 1] - u[i]  for i < n - 1,  (D u)[n - 1] = 0.
    // Its exact adjoint is D^T p = -(backward difference):
    //   (D^T p)[i] = p[i - 1] - p[i], with p[-1] = p[n - 1] = 0.
    // grad() applies D, and divergence() returns -D^T (backward differences), so
    // neg_divergence() = D^T and neg_laplacian() = D^T D (positive semidefinite).
    //
    // axis 0, 1, 2 differentiate along the fastest-varying n3, n2 and n1 index
    // respectively. `adjoint` selects D (false) or -D^T (true); `update` adds into
    // the output instead of overwriting it.
    template <typename T, int axis, bool adjoint, bool update>
    __global__ void diff_kernel(DevicePtr<const T> u, DevicePtr<T> out) {

        auto idx = Index3D();
        if (idx < u.dims()) {
            const unsigned n = (axis == 0)   ? u.dims().n3
                               : (axis == 1) ? u.dims().n2
                                             : u.dims().n1;
            const unsigned c = (axis == 0) ? idx.z : (axis == 1) ? idx.y : idx.x;

            // value of u shifted by `off` along the differentiated axis
            auto at = [&](int off) -> T {
                unsigned k = c + off;
                if (axis == 0) return u(idx.x, idx.y, k);
                if (axis == 1) return u(idx.x, k, idx.z);
                return u(k, idx.y, idx.z);
            };

            T val;
            if (!adjoint) {
                val = (c + 1 < n) ? at(1) - at(0) : T(0);
            } else {
                T cur = (c + 1 < n) ? at(0) : T(0);
                T prev = (c > 0) ? at(-1) : T(0);
                val = cur - prev;
            }
            if (update) {
                out[idx] += val;
            } else {
                out[idx] = val;
            }
        }
    }

    template <typename T>
    VecArray<T> grad(const DeviceArray<T> &u) {
        DeviceArray<T> du_x(u.dims());
        DeviceArray<T> du_y(u.dims());
        DeviceArray<T> du_z(u.dims());

        dim3 blockSize(1, 8, 32);
        dim3 gridSize = make_grid(u.dims(), blockSize);

        diff_kernel<T, 0, false, false><<<gridSize, blockSize>>>(u, du_x);
        diff_kernel<T, 1, false, false><<<gridSize, blockSize>>>(u, du_y);
        diff_kernel<T, 2, false, false><<<gridSize, blockSize>>>(u, du_z);
        SAFE_CALL(cudaGetLastError());
        SAFE_CALL(cudaDeviceSynchronize());

        VecArray<T> du;
        du[0] = std::move(du_x);
        du[1] = std::move(du_y);
        du[2] = std::move(du_z);
        return du;
    }
    // Explicit template instantiations
    template VecArray<float> grad(const DeviceArray<float> &u);
    template VecArray<double> grad(const DeviceArray<double> &u);

    // divergence := -D^T (backward differences)
    template <typename T>
    DeviceArray<T> divergence(const VecArray<T> &u) {

        auto dims = u[0].dims();
        DeviceArray<T> div_u(dims);

        dim3 blockSize(1, 8, 32);
        dim3 gridSize = make_grid(dims, blockSize);

        diff_kernel<T, 0, true, false><<<gridSize, blockSize>>>(u[0], div_u);
        diff_kernel<T, 1, true, true><<<gridSize, blockSize>>>(u[1], div_u);
        diff_kernel<T, 2, true, true><<<gridSize, blockSize>>>(u[2], div_u);
        SAFE_CALL(cudaGetLastError());
        return div_u;
    }
    // Explicit template instantiations
    template DeviceArray<float> divergence(const VecArray<float> &u);
    template DeviceArray<double> divergence(const VecArray<double> &u);

    // laplacian := divergence(grad(u)) = -D^T D (negative semidefinite)
    template <typename T>
    DeviceArray<T> laplacian(const DeviceArray<T> &u) {
        auto temp = grad(u);
        return divergence(temp);
    }
    // Explicit template instantiations
    template DeviceArray<float> laplacian(const DeviceArray<float> &u);
    template DeviceArray<double> laplacian(const DeviceArray<double> &u);

} // namespace tomocam::gpu::opt
