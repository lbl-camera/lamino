/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 *National Laboratory (subject to receipt of any required approvals from the
 *U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 *IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 *the U.S. Government has been granted for itself and others acting on its
 *behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 *to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#include <cuda/std/complex>

#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/nufft.h"
#include "gpu/polar_grid.h"
#include "gpu/vec_array.h"
#include "gpu/xmcdprojs.h"

namespace tomocam::gpu {

    template <typename T>
    using Complex = cuda::std::complex<T>;

    template <typename T>
    VecArray<T> sysmat(const VecArray<T> &x, const gpu::PolarGrid<T> &grid,
                       T gamma) {

        VecArray<T> Ax;

        // Uniform -> Polar (nufft3d2 output is complex)
        VecArray<Complex<T>> tmp;
        for (size_t i = 0; i < 3; ++i) {
            auto xcmplx = gpu::array::to_complex(x[i]);
            auto ccmplx = DeviceArray<Complex<T>>(grid.dims());
            gpu::nufft::nufft3d2(ccmplx, xcmplx, grid);
            tmp[i] = std::move(ccmplx);
        }
        // Apply (e ⊗ e) projection matrix in the polar domain
        xmcd_projection(tmp, grid, gamma);

        // Polar -> Uniform (nufft3d1), take real part
        for (size_t i = 0; i < 3; ++i) {
            auto xcmplx = DeviceArray<Complex<T>>(x[i].dims());
            gpu::nufft::nufft3d1(tmp[i], xcmplx, grid);
            Ax[i] = gpu::array::to_real(xcmplx);
        }
        return Ax;
    }
    // explicit instantiations
    template VecArray<float> sysmat(const VecArray<float> &x,
                                    const gpu::PolarGrid<float> &grid, float gamma);
    template VecArray<double> sysmat(const VecArray<double> &x,
                                     const gpu::PolarGrid<double> &grid,
                                     double gamma);

} // namespace tomocam::gpu
