#include <fftw3.h>
#include <omp.h>

#include <complex>
#include <iostream>

#include "dtypes.h"
#include "fft.h"

const complex_t J1(0, 1);

void fft1(Array < complex_t > &input) {

    dim2_t dims = input.dims();
    fftwf_complex *data = (fftwf_complex *) input.ptr();
    int rank = 1;
    int n[] = { dims.y };
    int idist = dims.y;
    int odist = dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed, istride,
                                          idist, data, onembed, ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}

void ifft1(Array < complex_t > &input) {

    dim2_t dims = input.dims();
    fftwf_complex *data = (fftwf_complex *) input.ptr();
    int rank = 1;
    int n[] = { dims.y };
    int idist = dims.y;
    int odist = dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed, istride,
                                          idist, data, onembed, ostride, odist,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}

void fft2(Array < complex_t > &input) {

    dim2_t dims = input.dims();
    fftwf_complex *data = (fftwf_complex *) input.ptr();
    int rank = 2;
    int n[] = { dims.x, dims.y };
    int idist = dims.x * dims.y;
    int odist = dims.x * dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = 1;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed, istride,
                                          idist, data, onembed, ostride, odist,
                                          FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}

void ifft2(Array < complex_t > &input)
{

    dim2_t dims = input.dims();
    fftwf_complex *data = (fftwf_complex *) input.ptr();
    int rank = 2;
    int n[] = { dims.x, dims.y };
    int idist = dims.x * dims.y;
    int odist = dims.x * dims.y;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = 1;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed, istride,
                                          idist, data, onembed, ostride, odist,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}

void fftshift1(Array < complex_t > &input) {
    dim2_t dims = input.dims();
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
			input(i, j) *= std::pow(-1, j & 1);
}

void fftshift2(Array < complex_t > &input) {
    dim2_t dims = input.dims();
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
			input(i, j) *= std::pow(-1, (i + j) & 1);
}

void xfftshift(Array < complex_t > &input, float x) {
    dim2_t dims = input.dims();
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++) {
			float w = (2 * M_PI * x * j) / (float) dims.y;
			input(i, j) *= std::exp(-J1 * w);
		}
}
