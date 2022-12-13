
#include <complex>
#include <fftw3.h>

#include "array.h"
#include "dtypes.h"

#ifndef FFTDEFS__H
#define FFTDEFS__H

Array<complex_t> fft1(const Array <complex_t> &);
Array<complex_t> ifft1(const Array <complex_t> &);
Array<complex_t> fft2(const Array <complex_t> &);
Array<complex_t> ifft2(const Array <complex_t> &);

Array<complex_t> fft2r2c(Array <float>);
Array<float> fft2c2r(Array <complex_t>);

template <typename T> 
void fftshift1(Array<T> & arr) {
    auto dims = arr.dims();
    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
            arr(i,j) *= std::pow(-1, j & 1);
}

template <typename T> 
void fftshift2(Array<T> & arr) {
    auto dims = arr.dims();
    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
            arr(i,j) *= std::pow(-1, (i + j) & 1);
}

inline void xfftshift(Array <complex_t> &arr, float d) {
    auto dims = arr.dims();
    #pragma omp parallel for 
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++) {
			float w = (2 * M_PI * d * j) / (float) dims.y;
			arr(i, j) *= std::exp(-I * w);
		}
}
#endif // FFTDEFS__H
