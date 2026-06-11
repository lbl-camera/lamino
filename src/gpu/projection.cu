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

#include <array>
#include <cuda/std/complex>

#include "dtypes.h"

#include "gpu/cufinufft_plan_cache.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/fft.h"
#include "gpu/fftshift.h"
#include "gpu/nufft.h"
#include "gpu/polar_grid.h"
#include "gpu/utils.h"

namespace tomocam::gpu {

    template <typename T>
    using complex = cuda::std::complex<T>;

    template <typename T>
    DeviceArray<T> forward(const std::array<DeviceArray<T>, 3> &m,
                           const gpu::PolarGrid<T> &pg) {

        // cast to complex
        auto dims = pg.dims();
        T scale = static_cast<T>(dims.n2 * dims.n3);
        auto proj = DeviceArray<complex<T>>(dims);

        for (size_t i = 0; i < 3; ++i) {
            auto mcmplx = array::to_complex(m[i]);
        
            DeviceArray<complex<T>> C(dims);

            // NUFFT type-2: uniform -> nonuniform
            nufft::nufft3d2(C, mcmplx, pg);

            // ifft with fftshift
            C = gpu::fftshift2(C);
            C = gpu::fft::ifft2d(C);
            C = gpu::ifftshift2(C);

            // accumulate projections

        return gpu::array::to_real(C);
    }

    // Explicit instantiations for forward
    template DeviceArray<float> forward(const DeviceArray<float> &,
                                        const gpu::PolarGrid<float> &);
    template DeviceArray<double> forward(const DeviceArray<double> &,
                                         const gpu::PolarGrid<double> &);

    // -------------------------------------------------------------------------
    // backward: projections -> volume
    // -------------------------------------------------------------------------

    template <typename T>
    DeviceArray<T> adjoint(const DeviceArray<T> &proj, const gpu::PolarGrid<T> &pg,
                           const dims_t &recon_dims, T gamma) {
        // cast to complex
        auto C = array::to_complex(proj);

        // 2D Fourier transform with fftshift
        C = gpu::fftshift2(C);
        C = gpu::fft::fft2d(C);
        C = gpu::ifftshift2(C);

        DeviceArray<complex<T>> F(recon_dims);
        nufft::nufft3d1(C, F, pg);

        // scale and return real part
        auto result = array::to_real(F);
        T scale = static_cast<T>(proj.nrows() * proj.ncols());
        return result / scale;
    }
    template DeviceArray<float> adjoint(const DeviceArray<float> &,
                                        const gpu::PolarGrid<float> &,
                                        const dims_t &, float);
    template DeviceArray<double> adjoint(const DeviceArray<double> &,
                                         const gpu::PolarGrid<double> &,
                                         const dims_t &, double);
} // namespace tomocam::gpu
