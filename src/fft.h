#include <complex>
#include <fftw3.h>
#include <new>
#include <type_traits>

#include "array.h"
#include "dtypes.h"

#ifndef FFTDEFS__H
#define FFTDEFS__H

namespace tomocam {
    namespace fft {

        template <typename T>
        Array<std::complex<T>> fft1(const Array<std::complex<T>> &input) {

            // double or float
            bool is_double = std::is_same_v<T, double>;

            // create return array
            dims_t dims = input.dims();
            Array<std::complex<T>> output(dims);

            int rank = 1;
            int n[] = {dims.z};
            int idist = dims.z;
            int odist = dims.z;
            int istride = 1;
            int ostride = 1;
            int *inembed = NULL;
            int *onembed = NULL;
            int batches = dims.x * dims.y;

            if (is_double) {
                fftw_complex *idata = (fftw_complex *)input.begin();
                fftw_complex *odata = (fftw_complex *)output.begin();
                fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_FORWARD, FFTW_ESTIMATE);
                if (!plan) { throw std::bad_alloc(); }
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            } else {
                fftwf_complex *idata = (fftwf_complex *)input.begin();
                fftwf_complex *odata = (fftwf_complex *)output.begin();
                fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_FORWARD, FFTW_ESTIMATE);
                fftwf_execute(plan);
                fftwf_destroy_plan(plan);
            }
            return output;
        }

        template <typename T>
        Array<std::complex<T>> ifft1(const Array<std::complex<T>> &arr) {

            // double or float
            bool is_double = std::is_same_v<T, double>;

            // create return array
            dims_t dims = arr.dims();
            Array<std::complex<T>> output(arr.dims());

            int rank = 1;
            int n[] = {dims.z};
            int idist = dims.z;
            int odist = dims.z;
            int istride = 1;
            int ostride = 1;
            int *inembed = NULL;
            int *onembed = NULL;
            int batches = dims.x * dims.y;

            if (is_double) {
                fftw_complex *idata = (fftw_complex *)arr.begin();
                fftw_complex *odata = (fftw_complex *)output.begin();
                fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_BACKWARD, FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            } else {
                fftwf_complex *idata = (fftwf_complex *)arr.begin();
                fftwf_complex *odata = (fftwf_complex *)output.begin();
                fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_BACKWARD, FFTW_ESTIMATE);
                fftwf_execute(plan);
                fftwf_destroy_plan(plan);
            }
            return output;
        }

        template <typename T>
        Array<std::complex<T>> fft2(const Array<std::complex<T>> &input) {

            // double or float
            bool is_double = std::is_same_v<T, double>;

            // create return array
            dims_t dims = input.dims();
            Array<std::complex<T>> output(input.dims());

            int rank = 2;
            int n[] = {dims.y, dims.z};
            int idist = dims.y * dims.z;
            int odist = dims.y * dims.z;
            int istride = 1;
            int ostride = 1;
            int *inembed = NULL;
            int *onembed = NULL;
            int batches = dims.x;

            if (is_double) {
                fftw_complex *idata = (fftw_complex *)input.begin();
                fftw_complex *odata = (fftw_complex *)output.begin();
                fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_FORWARD, FFTW_ESTIMATE);
                if (!plan) { throw std::bad_alloc(); }
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            } else {
                fftwf_complex *idata = (fftwf_complex *)input.begin();
                fftwf_complex *odata = (fftwf_complex *)output.begin();
                fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_FORWARD, FFTW_ESTIMATE);
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
            int n[] = {dims.y, dims.z};
            int idist = dims.y * dims.z;
            int odist = dims.y * dims.z;
            int istride = 1;
            int ostride = 1;
            int *inembed = NULL;
            int *onembed = NULL;
            int batches = dims.x;

            if (is_double) {
                fftw_complex *idata = (fftw_complex *)input.begin();
                fftw_complex *odata = (fftw_complex *)output.begin();
                fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_BACKWARD, FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            } else {
                fftwf_complex *idata = (fftwf_complex *)input.begin();
                fftwf_complex *odata = (fftwf_complex *)output.begin();
                fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata,
                    inembed, istride, idist, odata, onembed, ostride, odist,
                    FFTW_BACKWARD, FFTW_ESTIMATE);
                fftwf_execute(plan);
                fftwf_destroy_plan(plan);
            }
            return output;
        }

        template <typename T>
        Array<std::complex<T>> fft3(const Array<std::complex<T>> &input) {

            // double or float
            bool is_double = std::is_same_v<T, double>;

            // create return array
            dims_t dims = input.dims();
            Array<std::complex<T>> output(input.dims());

            if (is_double) {
                fftw_complex *idata = (fftw_complex *)input.begin();
                fftw_complex *odata = (fftw_complex *)output.begin();
                fftw_plan plan = fftw_plan_dft_3d(dims.x, dims.y, dims.z, idata,
                    odata, -1, FFTW_ESTIMATE);
                if (!plan) { throw std::bad_alloc(); }
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            } else {
                fftwf_complex *idata = (fftwf_complex *)input.begin();
                fftwf_complex *odata = (fftwf_complex *)output.begin();
                fftwf_plan plan = fftwf_plan_dft_3d(dims.x, dims.y, dims.z,
                    idata, odata, -1, FFTW_ESTIMATE);
                if (!plan) { throw std::bad_alloc(); }
                fftwf_execute(plan);
                fftwf_destroy_plan(plan);
            }
            return output;
        }

        template <typename T>
        Array<std::complex<T>> ifft3(const Array<std::complex<T>> &input) {

            // double or float
            bool is_double = std::is_same_v<T, double>;

            dims_t dims = input.dims();
            Array<std::complex<T>> output(dims);

            if (is_double) {
                fftw_complex *idata = (fftw_complex *)input.begin();
                fftw_complex *odata = (fftw_complex *)output.begin();
                fftw_plan plan = fftw_plan_dft_3d(dims.x, dims.y, dims.z, idata,
                    odata, 1, FFTW_ESTIMATE);
                if (!plan) { throw std::bad_alloc(); }
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            } else {
                fftwf_complex *idata = (fftwf_complex *)input.begin();
                fftwf_complex *odata = (fftwf_complex *)output.begin();
                fftwf_plan plan = fftwf_plan_dft_3d(dims.x, dims.y, dims.z,
                    idata, odata, 1, FFTW_ESTIMATE);
                if (!plan) { throw std::bad_alloc(); }
                fftwf_execute(plan);
                fftwf_destroy_plan(plan);
            }

            return output;
        }

    } // namespace fft
} // namespace tomocam
#endif // FFTDEFS__H
