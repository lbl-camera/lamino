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

#ifndef GPU_NUFFT_H
#define GPU_NUFFT_H

#include <cmath>
#include <cuda/std/complex>
#include <vector>

#include "array.h"
#include "cufinufft_plan.h"
#include "cufinufft_plan_cache.h"
#include "dtypes.h"
#include "polar_grid.h"

namespace tomocam::nufft {

    // 3D Type-1 NUFFT (GPU): nonuniform points to uniform grid
    template <typename T>
    void gpu_nufft3d1(const Array<std::complex<T>> &cz, Array<std::complex<T>> &fz,
                      const PolarGrid<T> &pg) {

        using complex_t = cuda::std::complex<T>;

        std::array<int64_t, 3> n_modes = {(int64_t)fz.ncols(), (int64_t)fz.nrows(),
                                          (int64_t)fz.nslices()};
        auto &plan = plans::cu_cache<T>.get_plan(1, 3, n_modes, 1);
        plan.set_points(pg);

        // input array on device
        DeviceArray<complex_t> d_cz(cz.dims());
        gpu::copy_to_device(d_cz.begin(), cz.begin(), cz.bytes());

        // output array on device
        DeviceArray<complex_t> d_fz(fz.dims());

        int ierr = plan.execute((std::complex<T> *)cz.begin(),
                                (std::complex<T> *)fz.begin());
        // copy
        gpu::copy_to_host(fz.begin(), d_fz.begin(), fz.bytes());

        if (ierr != 0) { throw std::runtime_error("Error in cufinufft_execute"); }
    }

    // 3D Type-2 NUFFT (GPU): uniform grid to nonuniform points
    template <typename T>
    void gpu_nufft3d2(Array<std::complex<T>> &cz, const Array<std::complex<T>> &fz,
                      const PolarGrid<T> &pg) {

        std::array<int64_t, 3> n_modes = {(int64_t)fz.ncols(), (int64_t)fz.nrows(),
                                          (int64_t)fz.nslices()};
        auto &plan = plans::cu_cache<T>.get_plan(2, 3, n_modes, -1);
        plan.set_points(pg);

        int ierr = plan.execute((std::complex<T> *)cz.begin(),
                                (std::complex<T> *)fz.begin());
        if (ierr != 0) { throw std::runtime_error("Error in cufinufft_execute"); }
    }

} // namespace tomocam::nufft

#endif // GPU_NUFFT__H
