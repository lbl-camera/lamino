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

#ifndef NUFFT__H
#define NUFFT__H

#include <cmath>
#include <complex>
#include <vector>

#include "array.h"
#include "dtypes.h"
#include "finufft_plan.h"
#include "finufft_plan_cache.h"
#include "polar_grid.h"

namespace tomocam::nufft {
    // 3D Type-1 NUFFT: nonuniform points to uniform grid
    template <typename T>
    void nufft3d1(const Array<std::complex<T>> &cz, Array<std::complex<T>> &fz,
                  const PolarGrid<T> &pg) {

        std::array<int64_t, 3> n_modes = {(int64_t)fz.ncols(), (int64_t)fz.nrows(),
                                          (int64_t)fz.nslices()};
        const int NUFFT_TYPE = 1;
        const int ndims = 3;
        const int sign = 1;
        auto &plan = plans::cache<T>.get_plan(NUFFT_TYPE, ndims, n_modes, sign);
        plan.set_points(pg);

        int ierr = plan.execute((std::complex<T> *)cz.begin(),
                                (std::complex<T> *)fz.begin());
        if (ierr != 0) { throw std::runtime_error("Error in finufft_execute"); }
    }

    // 3D Type-2 NUFFT: uniform grid to nonuniform points
    template <typename T>
    void nufft3d2(Array<std::complex<T>> &cz, const Array<std::complex<T>> &fz,
                  const PolarGrid<T> &pg) {

        std::array<int64_t, 3> n_modes = {(int64_t)fz.ncols(), (int64_t)fz.nrows(),
                                          (int64_t)fz.nslices()};
        const int NUFFT_TYPE = 2;
        const int ndims = 3;
        const int sign = -1;
        auto &plan = plans::cache<T>.get_plan(NUFFT_TYPE, ndims, n_modes, sign);
        plan.set_points(pg);

        int ierr = plan.execute((std::complex<T> *)cz.begin(),
                                (std::complex<T> *)fz.begin());
        if (ierr != 0) { throw std::runtime_error("Error in finufft_execute"); }
    }
} // namespace tomocam::nufft

#endif // NUFFT__H
