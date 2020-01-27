#ifndef FFTDEFS__H
#define FFTDEFS__H

#include <complex>
#include <fftw3.h>

#include "dtypes.h"


void fft1d(fftw_complex *, dim3_t);

void ifft2d(fftw_complex *, dim3_t);

void fftshift1d(complex_t *, dim3_t);

void fftshift2d(complex_t *, dim3_t);

void xfftshift(complex_t *, dim3_t, double);


#endif // FFTDEFS__H
