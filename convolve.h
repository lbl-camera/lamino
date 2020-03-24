

#include "dtypes.h"
#include "array.h"
#include "kernel.h"

#ifndef CONVOLVE__H
#define CONVOLVE__H


void convolve_p2c(Array_ref<complex_t>, Array_ref<complex_t>, Array_ref<float>, Kernel);
void convolve_c2p(Array_ref<complex_t>, Array_ref<complex_t>, Array_ref<float>, Kernel);
void deconvolve(Array_ref<complex_t>, Kernel);

void forward(Array_ref<complex_t>, Array_ref<complex_t>, Array_ref<float>, Kernel);
void backward(Array_ref<complex_t>, Array_ref<complex_t>, Array_ref<float>, Kernel);

double calc_error(float *, Array_ref<complex_t>, int);
void write_output(Array_ref<complex_t>, int);

#endif // CONVOLVE__H
