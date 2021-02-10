
#include "array.h"
#include "dtypes.h"
#include "kernel.h"

#ifndef CONVOLVE__H
#define CONVOLVE__H

void convolve_p2c(Array < complex_t > &, Array < complex_t > &,
                  Array < float > &, Kernel, float);
void convolve_c2p(Array < complex_t > &, Array < complex_t > &,
                  Array < float > &, Kernel, float);
void deconvolve(Array < complex_t > &, Kernel);

void forward(Array < complex_t > &, Array < complex_t > &,
             Array < float > &, Kernel, float center = 0.f);
void backward(Array < complex_t > &, Array < complex_t > &,
              Array < float > &, Kernel, float center = 0.f);


void ramp_filter(Array <complex_t> &, float);
void remove_padding(float *, Array < complex_t > &, int, int);
void normalize(float *, int);
void write_output(float *, int);
float calc_error(float *, float *, int);

#endif // CONVOLVE__H
