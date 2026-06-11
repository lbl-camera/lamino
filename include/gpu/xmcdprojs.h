

#include <cuda/std/complex>

#include "gpu/utils.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu {

    template <typename T>
    using Complex = cuda::std::complex<T>;

    template <typename T>
    void project_component(DeviceArray<Complex<T>> &m_comp,
                           const gpu::PolarGrid<T> &grid, T gamma, int comp);

    template <typename T>
    void xmcd_projection(VecArray<Complex<T>> &m, const gpu::PolarGrid<T> &grid,
                         T gamma);
} // namespace tomocam::gpu
