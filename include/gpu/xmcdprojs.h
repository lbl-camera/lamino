

#include <cuda/std/complex>

#include "gpu/utils.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu {

    template <typename T>
    using Complex = cuda::std::complex<T>;

    // Literal port of tomocam::beam_dir_vector(theta,gamma,beta) (see
    // include/projection.h): third column of RotationTranspose(theta,gamma,beta).
    template <typename T>
    __device__ __forceinline__ T beam_dir_vector(T theta, T gamma, T beta,
                                                 int comp) {
        using cuda::std::cos;
        using cuda::std::sin;
        T st = sin(theta), ct = cos(theta);
        T sb = sin(beta), cb = cos(beta);
        T sg = sin(gamma), cg = cos(gamma);
        switch (comp) {
            case 0: return st * sg - ct * sb * cg;
            case 1: return st * cg + ct * sb * sg;
            case 2: return ct * cb;
            default: return T(0);
        }
    }

    template <typename T>
    void project_component(DeviceArray<Complex<T>> &m_comp,
                           const gpu::PolarGrid<T> &grid, int comp);

    template <typename T>
    void xmcd_projection(VecArray<Complex<T>> &m, const gpu::PolarGrid<T> &grid);
} // namespace tomocam::gpu
