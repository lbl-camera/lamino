#include <fftw3.h>
#include <omp.h>

#include <iostream>
#include <complex>

#include "dtypes.h"
#include "fft.h"


void fft1d(fftw_complex * input, dim3_t dims) {
    int rank = 1;
    int n[] = { dims.z };
    int idist = dims.z;
    int odist = dims.z;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x * dims.y;
    fftw_plan plan = fftw_plan_many_dft(rank, n, batches, input, inembed,
        istride, idist, input, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

void ifft2d(fftw_complex * input, dim3_t dims) {
    int rank = 2;
    int n[] = { dims.y, dims.z };
    int idist = dims.y * dims.z;
    int odist = dims.y * dims.z;
    int istride = 1;
    int ostride = 1;
    int *inembed = NULL;
    int *onembed = NULL;
    int batches = dims.x;
    fftw_plan plan = fftw_plan_many_dft(rank, n, batches, input, inembed,
        istride, idist, input, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

void fftshift1d(complex_t * input, dim3_t dims) {
    int batches = dims.x * dims.y;

    #pragma omp parallel for
    for (int i = 0; i< batches; i++)
        for (int j = 0; j < dims.z; j++) {
            double a = std::pow(-1, j & 1);
            int id = i * dims.z + j;
            input[id] *= a;
        }
}

void fftshift2d(complex_t * input, dim3_t dims) {
    int ny = dims.y;
    int nz = dims.z;
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) 
        for (int k = 0; k < nz; k++) {
            double a = std::pow(-1, (j+k)&1);
            int id = j * nz + k;
            input[id] *= a;
        }
}

void xfftshift(complex_t * input, dim3_t dims, double shift) {
    int batches = dims.x * dims.y;
    #pragma omp parallel for
    for (int i = 0; i< batches; i++)
        for (int j = 0; j < dims.z; j++) {
            double d = 2 * M_PI * shift * j /((double) dims.z);
            int id = i * dims.z + j;
            complex_t I1 = complex_t(0, 1);
            input[id] = input[id] * exp(-I1 * d);
        }
}
