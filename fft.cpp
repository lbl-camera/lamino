#include <fftw3.h>
#include <omp.h>

#include <iostream>
#include <complex>

#include "dtypes.h"
#include "fft.h"


void fft1(Array<complex_t> &input) {

    dim3_t dims = input.dims();
    fftwf_complex * data = (fftwf_complex *) input.ptr();
    int rank = 1;
    int n[] = { dims.z };
    int idist = dims.z;
    int odist = dims.z;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x * dims.y;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed,
        istride, idist, data, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}
void ifft1(Array<complex_t> &input) {

    dim3_t dims = input.dims();
    fftwf_complex * data = (fftwf_complex *) input.ptr();
    int rank = 1;
    int n[] = { dims.z };
    int idist = dims.z;
    int odist = dims.z;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x * dims.y;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed,
        istride, idist, data, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}
void fft2(Array<complex_t> &input) {

    dim3_t dims = input.dims();
    fftwf_complex * data = (fftwf_complex *) input.ptr();
    int rank = 2;
    int n[] = { dims.y, dims.z };
    int idist = dims.y * dims.z;
    int odist = dims.y * dims.z;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed,
        istride, idist, data, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}
void ifft2(Array<complex_t> &input) {

    dim3_t dims = input.dims();
    fftwf_complex * data = (fftwf_complex *) input.ptr();
    int rank = 2;
    int n[] = { dims.y, dims.z };
    int idist = dims.y * dims.z;
    int odist = dims.y * dims.z;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftwf_plan plan = fftwf_plan_many_dft(rank, n, batches, data, inembed,
        istride, idist, data, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}

void fftshift1(Array<complex_t> &input) {
    dim3_t dims = input.dims();
    #pragma omp parallel for collapse(2)
    for (int i = 0; i< dims.x; i++)
        for (int j = 0; j < dims.y; j++) 
                for (int k = 0; k < dims.z; k++)
                    input(i,j,k) *= std::pow(-1, k & 1);
}

void fftshift2(Array<complex_t> &input) {
    dim3_t dims = input.dims();
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
            for (int k = 0; k < dims.z; k++)
                input(i,j,k) *= std::pow(-1, (j + k) & 1);
}
