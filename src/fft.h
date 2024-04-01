
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

Array<complex_t> fft2r2c(Array <double>);
Array<double> fft2c2r(Array <complex_t>);


inline int fftwsize(int n) {
    int primes[] = {3, 5, 7, 11};
    int pows[] =   {0, 0, 0, 0};
    int np = 4;

    while (true) {
        auto N = n;
        for (int i = 0; i < np; i++) {
            auto a = primes[i];
            while (N % a == 0) {
                N /= a;
                pows[i]++;
            }
            if (N == 1) break;
        }
        if (N > 1) {
             n++;
             memset(pows, 0, sizeof(int)*np);
        } else if (pows[np-2] + pows[np-1] > 1) {
            n++;
             memset(pows, 0, sizeof(int)*np);
        } else
            break;
    }
    return n;
}

template <typename T> 
void fftshift1(Array<T> & arr) {
    auto dims = arr.dims();
    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
            arr(i,j) *= std::pow(-1, (j & 1));
}

template <typename T> 
void fftshift2(Array<T> & arr) {
    auto dims = arr.dims();
    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
            arr(i,j) *= std::pow(-1, (i + j) & 1);
}


template <typename T>
Array<T> fftshift(const Array<T> & arr) {

    auto dims = arr.dims();
    int dy = dims.y / 2;

    Array<T> out(dims.x, dims.y);
    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++) {
        for (int j = 0; j < dims.y; j++) {
            int jj = (j + dy) % dims.y;
            out(i,jj) = arr(i,j);
        }
    }
    return out;
}

template <typename T>
Array<T> ifftshift(const Array<T> & arr) {

    auto dims = arr.dims();
    int dy = dims.y / 2;

    Array<T> out(dims.x, dims.y);
    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++) {
        for (int j = 0; j < dims.y; j++) {
            int jj = (j + dy) % dims.y;
            out(i,j) = arr(i,jj);
        }
    }
    return out;
}   

inline void xfftshift(Array <complex_t> &arr, double d) {
    auto dims = arr.dims();
    #pragma omp parallel for 
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++) {
			double w = (2 * M_PI * d * j) / (double) dims.y;
			arr(i, j) *= std::exp(I * w);
		}
}
#endif // FFTDEFS__H
