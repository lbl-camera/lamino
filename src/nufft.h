
#include <cmath>
#include <complex>
#include <finufft.h>
#include <vector>

#include "array.h"
#include "dtypes.h"
#include "polar_grid.h"

#ifndef NUFFT__H
#define NUFFT__H

namespace tomocam {
    namespace nufft {

        template <typename T>
        void nufft3d1(const Array<std::complex<T>> &cz,
            Array<std::complex<T>> &fz, const PolarGrid<T> &pg) {

            int M = pg.npts();
            int N1 = fz.nslices();
            int N2 = fz.nrows();
            int N3 = fz.ncols();
            finufft_opts *opts = new finufft_opts;

            if (std::is_same_v<T, double>) {
                finufft_default_opts(opts);
                opts->upsampfac = 2.0;

                T tol = 1e-15;
                double *x = (double *)pg.x.begin();
                double *y = (double *)pg.y.begin();
                double *z = (double *)pg.z.begin();
                std::complex<double> *cptr = (std::complex<double> *)cz.begin();
                std::complex<double> *fptr = (std::complex<double> *)fz.begin();
                finufft3d1(M, x, y, z, cptr, 1, tol, N1, N2, N3, fptr, opts);
            } else {
                finufftf_default_opts(opts);
                opts->upsampfac = 2.0;

                T tol = 1e-07;
                float *x = (float *)pg.x.get();
                float *y = (float *)pg.y.get();
                float *z = (float *)pg.z.get();
                std::complex<float> *cptr = (std::complex<float> *)cz.begin();
                std::complex<float> *fptr = (std::complex<float> *)fz.begin();
                finufftf3d1(M, x, y, z, cptr, 1, tol, N1, N2, N3, fptr, opts);
            }
            delete[] opts;
        }

        template <typename T>
        void nufft3d2(const Array<std::complex<T>> &cz,
            const Array<std::complex<T>> &fz, const PolarGrid<T> &pg) {

            int M = pg.npts();
            int N1 = fz.nslices();
            int N2 = fz.nrows();
            int N3 = fz.ncols();
            finufft_opts *opts = new finufft_opts;

            if (std::is_same_v<T, double>) {
                finufft_default_opts(opts);
                opts->upsampfac = 2.0;

                T tol = 1e-15;
                double *x = (double *)pg.x.get();
                double *y = (double *)pg.y.get();
                double *z = (double *)pg.z.get();
                std::complex<double> *cptr = (std::complex<double> *)cz.begin();
                std::complex<double> *fptr = (std::complex<double> *)fz.begin();
                finufft3d2(M, x, y, z, cptr, tol, N1, N2, N3, fptr);
            } else {
                finufftf_default_opts(opts);
                opts->upsampfac = 2.0;

                T tol = 1e-07;
                float *x = (float *)pg.x.get();
                float *y = (float *)pg.y.get();
                float *z = (float *)pg.z.get();
                std::complex<float> *cptr = (std::complex<float> *)cz.begin();
                std::complex<float> *fptr = (std::complex<float> *)fz.begin();
                finufftf3d2(M, x, y, z, cptr, -1, tol, N1, N2, N3, fptr, opts);
            }

            delete[] opts;
        }

    } // namespace nufft
} // namespace tomocam

#endif // NUFFT__H
