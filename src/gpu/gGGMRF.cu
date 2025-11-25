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

#include <cuda_runtime.h>

#include "dev_array.h"
#include "gpu/utils.cuh"

#include "potential_function.cuh"

namespace tomocam::gpu {

    constexpr int N1 = 2;
    constexpr int N2 = 16;
    constexpr int N3 = 16; 
    constexpr SHAMEM = (N1 + 2) * (N2 + 2) * (N3 + 2);

    template <typename T>
    __global__ void qGGMRF_kernel(DevicePtr<T> f_ext, DevicePtr<T> g, T sigma, T p) {

        // shared memory for the block
        __shared__ T s_val[SHAMEM];

        // get 3D index
        auto idx = Index3D();
        if (idx < g.dims()) {

            // last block with shift less thet block size
            int i_shift = min(blockDim.z, g.dims().n1 - blockIdx.z * blockDim.z);
            int j_shift = min(blockDim.y, g.dims().n2 - blockIdx.y * blockDim.y);
            int k_shift = min(blockDim.x, g.dims().n3 - blockIdx.x * blockDim.x);
            auto flatidx = [&] (int i, int j, int k) {
               return (i * (N2 + 2) + j) * (N3 + 2) + k; 
            };

            /* copy values into shared memory. */
            for (int i = threadIdx.z; i < N1 + 2; i += i_shift) {
                for (int j = threadIdx.y; j < N2 + 2; j += j_shift) {
                    for (int k = threadIdx.x; k < N3 + 2; k += k_shift) {
                        int x = (int)(blockIdx.z * blockDim.z) + i - 1;
                        int y = (int)(blockIdx.y * blockDim.y) + j - 1;
                        int z = (int)(blockIdx.x * blockDim.x) + k - 1;
                        int idx = flatidx(i, j, k);
                        s_val[idx] = f_ext.at(x, y, z);
                    }
                }
            }
            __syncthreads();

            // compute the qGGMRF contribution
            T v = s_val[threadIdx.z + 1][threadIdx.y + 1][threadIdx.x + 1];
            T temp = 0.f;
            for (int ix = 0; ix < 3; ix++) {
                for (int iy = 0; iy < 3; iy++) {
                    for (int iz = 0; iz < 3; iz++) {
                        if (ix == 1 && iy == 1 && iz == 1)
                            continue;
                        auto idx = flatidx(threadIdx.z + ix, threadIdx.y + iy,
                                              threadIdx.x + iz);
                        auto delta = v - s_val[idx];
                        auto tv = d_pot_func(delta, p, sigma);
                        temp += weight(ix, iy, iz) * tv;
                    }
                }
            }
            g[idx] += temp;
        }
        __syncthreads();
    }

    template <typename T>
    void add_qGRRMRF(const DeviceArray<T> &sol, DeviceArray<T> &grad, T sigma,
                        T p) {

        // data size
        auto dims = grad.dims();
        // CUDA kernel parameters
        dim3 block(N3, N2, N1);
        dim3 grid;
        grid.x = divup(dims.n3, block.x);
        grid.y = DIVUP(dims.n2, block.y);
        grid.z = DIVUP(dims.n1, block.z);

        // update gradients inplace
        qGGMRF_kernel<T><<<grid, block>>>(sol, grad, sigma, p);
        SAFE_CALL(cudaGetLastError());
    }

    // instantiate the template
    template void add_qGRRMRF<float>(const DeviceArray<float> &sol,
                                    DeviceArray<float> &grad, float sigma,
                                    float p);
    template void add_qGRRMRF<double>(const DeviceArray<double> &sol,
                                     DeviceArray<double> &grad, double sigma,
                                     double p);

} // namespace tomocam::gpu
