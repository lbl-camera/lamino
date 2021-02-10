#ifndef FFTDEFS__H
#define FFTDEFS__H

#include <complex>
#include <fftw3.h>

#include "array.h"
#include "dtypes.h"

void fft1(Array < complex_t > &);
void ifft1(Array < complex_t > &);
void fft2(Array < complex_t > &);
void ifft2(Array < complex_t > &);

void fftshift1(Array < complex_t > &);
void fftshift2(Array < complex_t > &);
void xfftshift(Array < complex_t > &, float);

#endif                          // FFTDEFS__H
