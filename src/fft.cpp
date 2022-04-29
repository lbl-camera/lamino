#include <fftw3.h>
#include <omp.h>

#include <complex>
#include <iostream>

#include "dtypes.h"
#include "fft.h"

const complex_t J1(0, 1);

Array<complex_t> fft1(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftwf_complex *idata = (fftwf_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftwf_complex *odata = (fftwf_complex *) output.ptr();
    int rank = 1;
    int n[] = { dims.y };
    int idist = dims.y;
    int odist = dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    return output;
}

Array<complex_t> ifft1(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftwf_complex *idata = (fftwf_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftwf_complex *odata = (fftwf_complex *) output.ptr();
    int rank = 1;
    int n[] = { dims.y };
    int idist = dims.y;
    int odist = dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    return output;
}

Array<complex_t> fft2(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftwf_complex *idata = (fftwf_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftwf_complex *odata = (fftwf_complex *) output.ptr();
    int rank = 2;
    int n[] = { dims.x, dims.y };
    int idist = dims.x * dims.y;
    int odist = dims.x * dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = 1;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    return output;
}

Array<complex_t> ifft2(const Array <complex_t> &input) {

    dim2_t dims = input.dims();
    fftwf_complex *idata = (fftwf_complex *) input.ptr();
    Array<complex_t> output(dims.x, dims.y);
    fftwf_complex *odata = (fftwf_complex *) output.ptr();
    int rank = 2;
    int n[] = { dims.x, dims.y };
    int idist = dims.x * dims.y;
    int odist = dims.x * dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = 1;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, idata, inembed, istride,
                                          idist, odata, onembed, ostride, odist,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    return output;
}
