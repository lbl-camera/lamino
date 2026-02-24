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

#include <complex>
#include <fftw3.h>
#include <new>
#include <type_traits>

#include "array.h"
#include "dtypes.h"

#ifndef FFTDEFS_H
#define FFTDEFS_H

namespace tomocam::fft {

    template <typename T>
    Array<std::complex<T>> fft2(const Array<std::complex<T>> &input) {

        // double or float
        bool is_double = std::is_same_v<T, double>;

        // create return array
        dims_t dims = input.dims();
        Array<std::complex<T>> output(input.dims());

        int rank = 2;
        int n[] = {static_cast<int>(dims.n2), static_cast<int>(dims.n3)};
        int idist = static_cast<int>(dims.n2 * dims.n3);
        int odist = static_cast<int>(dims.n2 * dims.n3);
        int istride = 1;
        int ostride = 1;
        int *inembed = NULL;
        int *onembed = NULL;
        int batches = static_cast<int>(dims.n1);

        if (is_double) {
            auto *idata = (fftw_complex *)input.begin();
            auto *odata = (fftw_complex *)output.begin();
            fftw_plan plan = fftw_plan_many_dft(
                rank, n, batches, idata, inembed, istride, idist, odata, onembed,
                ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
            if (!plan) { throw std::bad_alloc(); }
            fftw_execute(plan);
            fftw_destroy_plan(plan);
        } else {
            auto *idata = (fftwf_complex *)input.begin();
            auto *odata = (fftwf_complex *)output.begin();
            fftwf_plan plan = fftwf_plan_many_dft(
                rank, n, batches, idata, inembed, istride, idist, odata, onembed,
                ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        }
        return output;
    }

    template <typename T>
    Array<std::complex<T>> ifft2(const Array<std::complex<T>> &input) {
        // double or float
        bool is_double = std::is_same_v<T, double>;

        // create return array
        dims_t dims = input.dims();
        Array<std::complex<T>> output(input.dims());

        int rank = 2;
        int n[] = {static_cast<int>(dims.n2), static_cast<int>(dims.n3)};
        int idist = static_cast<int>(dims.n2 * dims.n3);
        int odist = static_cast<int>(dims.n2 * dims.n3);
        int istride = 1;
        int ostride = 1;
        int *inembed = NULL;
        int *onembed = NULL;
        int batches = static_cast<int>(dims.n1);

        if (is_double) {
            auto *idata = (fftw_complex *)input.begin();
            auto *odata = (fftw_complex *)output.begin();
            fftw_plan plan = fftw_plan_many_dft(
                rank, n, batches, idata, inembed, istride, idist, odata, onembed,
                ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);
        } else {
            auto *idata = (fftwf_complex *)input.begin();
            auto *odata = (fftwf_complex *)output.begin();
            fftwf_plan plan = fftwf_plan_many_dft(
                rank, n, batches, idata, inembed, istride, idist, odata, onembed,
                ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        }
        return output;
    }

    template <typename T>
    Array<std::complex<T>> fft2_r2c(const Array<T> &input) {

        // double or float
        bool is_double = std::is_same_v<T, double>;

        // create return array
        dims_t dims = input.dims();
        Array<std::complex<T>> output(dims.n1, dims.n2, dims.n3 / 2 + 1);

        int rank = 2;
        int n[] = {static_cast<int>(dims.n2), static_cast<int>(dims.n3)};
        int idist = static_cast<int>(dims.n2 * dims.n3);
        int odist = static_cast<int>(dims.n2 * (dims.n3 / 2 + 1));
        int istride = 1;
        int ostride = 1;
        int *inembed = NULL;
        int *onembed = NULL;
        int batches = static_cast<int>(dims.n1);

        if (is_double) {
            auto *idata = (double *)input.begin();
            auto *odata = (fftw_complex *)output.begin();
            fftw_plan plan = fftw_plan_many_dft_r2c(rank, n, batches, idata, inembed,
                                                    istride, idist, odata, onembed,
                                                    ostride, odist, FFTW_ESTIMATE);
            if (!plan) { throw std::bad_alloc(); }
            fftw_execute(plan);
            fftw_destroy_plan(plan);
        } else {
            auto *idata = (float *)input.begin();
            auto *odata = (fftwf_complex *)output.begin();
            fftwf_plan plan = fftwf_plan_many_dft_r2c(
                rank, n, batches, idata, inembed, istride, idist, odata, onembed,
                ostride, odist, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        }
        return output;
    }

    template <typename T>
    Array<T> fft2_c2r(const Array<std::complex<T>> &input, dims_t output_dims) {
        // double or float
        bool is_double = std::is_same_v<T, double>;

        // create return array
        Array<T> output(output_dims);

        int rank = 2;
        int n[] = {static_cast<int>(output_dims.n2),
                   static_cast<int>(output_dims.n3)};
        int idist = static_cast<int>(output_dims.n2 * (output_dims.n3 / 2 + 1));
        int odist = static_cast<int>(output_dims.n2 * output_dims.n3);
        int istride = 1;
        int ostride = 1;
        int *inembed = NULL;
        int *onembed = NULL;
        int batches = static_cast<int>(output_dims.n1);

        if (is_double) {
            auto *idata = (fftw_complex *)input.begin();
            auto *odata = (double *)output.begin();
            fftw_plan plan = fftw_plan_many_dft_c2r(rank, n, batches, idata, inembed,
                                                    istride, idist, odata, onembed,
                                                    ostride, odist, FFTW_ESTIMATE);
            if (!plan) { throw std::bad_alloc(); }
            fftw_execute(plan);
            fftw_destroy_plan(plan);
        } else {
            auto *idata = (fftwf_complex *)input.begin();
            auto *odata = (float *)output.begin();
            fftwf_plan plan = fftwf_plan_many_dft_c2r(
                rank, n, batches, idata, inembed, istride, idist, odata, onembed,
                ostride, odist, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
        }
        return output;
    }
} // namespace tomocam::fft

#endif // FFTDEFS_H
