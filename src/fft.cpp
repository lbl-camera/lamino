#include <fftw3.h>
#include <omp.h>

#include <complex>
#include <iostream>

#include "dtypes.h"
#include "fft.h"

const complex_t J1(0, 1);

Array<complex_t> fft1(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftw_complex *idata = (fftw_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftw_complex *odata = (fftw_complex *) output.ptr();
    int rank = 1;
    int n[] = { dims.y };
    int idist = dims.y;
    int odist = dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return output;
}

Array<complex_t> ifft1(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftw_complex *idata = (fftw_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftw_complex *odata = (fftw_complex *) output.ptr();
    int rank = 1;
    int n[] = { dims.y };
    int idist = dims.y;
    int odist = dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return output;
}

Array<complex_t> fft2(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftw_complex *idata = (fftw_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftw_complex *odata = (fftw_complex *) output.ptr();
    int rank = 2;
    int n[] = { dims.x, dims.y };
    int idist = dims.x * dims.y;
    int odist = dims.x * dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = 1;
    fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return output;
}

Array<complex_t> ifft2(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftw_complex *idata = (fftw_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftw_complex *odata = (fftw_complex *) output.ptr();
    int rank = 2;
    int n[] = { dims.x, dims.y };
    int idist = dims.x * dims.y;
    int odist = dims.x * dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = 1;
    fftw_plan plan = fftw_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return output;
}

Array<complex_t> fft2r2c(Array<double> input) {
    dim2_t dims = input.dims();
    double *idata = (double *) input.ptr();

    // allocate return array
    Array<complex_t> output(dims.x, dims.y/2+1);
    fftw_complex * odata = (fftw_complex *) output.ptr();

    fftw_plan plan = fftw_plan_dft_r2c_2d(dims.x, dims.y, idata, odata, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return output;
}

Array<double> fft2c2r(Array<complex_t> input) {
    dim2_t dims = input.dims();
    fftw_complex * idata = (fftw_complex *) input.ptr();

    // allocate return array
    Array<double> output(dims.x, 2*(dims.y-1));
    double * odata = (double *) output.ptr();
    fftw_plan plan = fftw_plan_dft_c2r_2d(dims.x, 2*(dims.y-1), idata, odata, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return output;
}
