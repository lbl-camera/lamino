#include <cuda/std/complex>
#include <thrust/device_vector.h>

#include "gpu/device_ptr.h"
#include "gpu/polar_grid.h"
#include "gpu/utils.h"
#include "gpu/vec_array.h"
#include "gpu/xmcdprojs.h"

namespace tomocam::gpu {

    template <typename T>
    using Complex = cuda::std::complex<T>;

    // applied during in polar domain for forward and adjoint XMCD projections.;
    template <typename T>
    __global__ void projection_kernel(DevicePtr<Complex<T>> m_comp,
                                      const T *proj_angles, const T *proj_gammas,
                                      const T *proj_betas, int comp) {
        auto idx = Index3D();
        if (idx < m_comp.dims()) {
            T theta = proj_angles[idx.x];
            T gamma = proj_gammas[idx.x];
            T beta = proj_betas[idx.x];
            T e = beam_dir_vector(theta, gamma, beta, comp);
            m_comp[idx] *= e;
        }
    }

    template <typename T>
    void project_component(DeviceArray<Complex<T>> &m_comp,
                           const gpu::PolarGrid<T> &grid, int comp) {
        dim3 threads(1, 16, 16);
        auto g = make_grid(m_comp.dims(), threads);
        projection_kernel<T><<<g, threads>>>(
            m_comp, thrust::raw_pointer_cast(grid.angles.data()),
            thrust::raw_pointer_cast(grid.d_gammas.data()),
            thrust::raw_pointer_cast(grid.d_betas.data()), comp);
        SAFE_CALL(cudaGetLastError());
    }
    // explicit instantiations for projection of individual components.
    template void project_component<float>(DeviceArray<Complex<float>> &,
                                         const gpu::PolarGrid<float> &, int);
    template void project_component<double>(DeviceArray<Complex<double>> &,
                                          const gpu::PolarGrid<double> &, int);

    // Applies the (e ⊗ e) projection matrix in-place for every element.
    // idx.x is the projection-angle index; e is a real vector, data is complex.
    template <typename T>
    __global__ void
    xmcd_projection_kernel(DevicePtr<Complex<T>> mx, DevicePtr<Complex<T>> my,
                           DevicePtr<Complex<T>> mz, const T *proj_angles,
                           const T *proj_gammas, const T *proj_betas) {
        auto idx = Index3D();
        if (idx < mx.dims()) {
            T theta = proj_angles[idx.x];
            T gamma = proj_gammas[idx.x];
            T beta = proj_betas[idx.x];
            Complex<T> tx = mx[idx];
            Complex<T> ty = my[idx];
            Complex<T> tz = mz[idx];
            T e[3];
            for (int i = 0; i < 3; ++i) {
                e[i] = beam_dir_vector(theta, gamma, beta, i);
            }
            mx[idx] = tx * e[0] * e[0] + ty * e[0] * e[1] + tz * e[0] * e[2];
            my[idx] = tx * e[1] * e[0] + ty * e[1] * e[1] + tz * e[1] * e[2];
            mz[idx] = tx * e[2] * e[0] + ty * e[2] * e[1] + tz * e[2] * e[2];
        }
    }

    template <typename T>
    void xmcd_projection(VecArray<Complex<T>> &m, const gpu::PolarGrid<T> &grid) {
        dim3 threads(1, 16, 16);
        auto g = make_grid(m[0].dims(), threads);
        xmcd_projection_kernel<T><<<g, threads>>>(
            m[0], m[1], m[2], thrust::raw_pointer_cast(grid.angles.data()),
            thrust::raw_pointer_cast(grid.d_gammas.data()),
            thrust::raw_pointer_cast(grid.d_betas.data()));
        SAFE_CALL(cudaGetLastError());
    }

    template void xmcd_projection<float>(VecArray<Complex<float>> &,
                                         const gpu::PolarGrid<float> &);
    template void xmcd_projection<double>(VecArray<Complex<double>> &,
                                          const gpu::PolarGrid<double> &);

} // namespace tomocam::gpu
