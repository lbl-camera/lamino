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
#include <finufft.h>
#include <vector>

#include "array.h"
#include "dtypes.h"
#include "polar_grid.h"

namespace tomocam::nufft {
    // 3D Type-1 NUFFT: nonuniform points to uniform grid
    template <typename T>
    void nufft3d1(const Array<std::complex<T>> &cz, Array<std::complex<T>> &fz,
                  const PolarGrid<T> &pg) {
        int M = pg.npts;
        int N1 = fz.nslices();
        int N2 = fz.nrows();
        int N3 = fz.ncols();
        finufft_opts *opts = new finufft_opts;

        if (std::is_same_v<T, double>) {
            finufft_default_opts(opts);
            opts->upsampfac = 2.0;
            T tol = 1e-14;
            double *x = (double *)pg.x.begin();
            double *y = (double *)pg.y.begin();
            double *z = (double *)pg.z.begin();
            std::complex<double> *cptr = (std::complex<double> *)cz.begin();
            std::complex<double> *fptr = (std::complex<double> *)fz.begin();
            finufft3d1(M, x, y, z, cptr, 1, tol, N3, N2, N1, fptr, opts);
        } else {
            finufftf_default_opts(opts);
            opts->upsampfac = 2.0;
            T tol = 1.2e-06;
            float *x = (float *)pg.x.begin();
            float *y = (float *)pg.y.begin();
            float *z = (float *)pg.z.begin();
            std::complex<float> *cptr = (std::complex<float> *)cz.begin();
            std::complex<float> *fptr = (std::complex<float> *)fz.begin();
            finufftf3d1(M, x, y, z, cptr, 1, tol, N3, N2, N1, fptr, opts);
        }
        delete opts;
    }

    // 3D Type-2 NUFFT: uniform grid to nonuniform points
    template <typename T>
    void nufft3d2(Array<std::complex<T>> &cz, const Array<std::complex<T>> &fz,
                  const PolarGrid<T> &pg) {
        auto M = pg.npts;
        auto N1 = fz.nslices();
        auto N2 = fz.nrows();
        auto N3 = fz.ncols();
        finufft_opts *opts = new finufft_opts;

        if (std::is_same_v<T, double>) {
            finufft_default_opts(opts);
            opts->upsampfac = 2.0;
            T tol = 1e-14;
            double *x = (double *)pg.x.begin();
            double *y = (double *)pg.y.begin();
            double *z = (double *)pg.z.begin();
            std::complex<double> *cptr = (std::complex<double> *)cz.begin();
            std::complex<double> *fptr = (std::complex<double> *)fz.begin();
            finufft3d2(M, x, y, z, cptr, -1, tol, N3, N2, N1, fptr, opts);
        } else {
            finufftf_default_opts(opts);
            opts->upsampfac = 2.0;
            T tol = 1.2e-06;
            float *x = (float *)pg.x.begin();
            float *y = (float *)pg.y.begin();
            float *z = (float *)pg.z.begin();
            std::complex<float> *cptr = (std::complex<float> *)cz.begin();
            std::complex<float> *fptr = (std::complex<float> *)fz.begin();
            finufftf3d2(M, x, y, z, cptr, -1, tol, N3, N2, N1, fptr, opts);
        }
        delete opts;
    }
} // namespace tomocam::nufft

#endif // NUFFT__H
