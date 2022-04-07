
#include "array.h"
#include "dtypes.h"

#ifndef TOMOCAM__H
#define TOMOCAM__H

Array<float> forward(Array <float> &, const Array <float> &, float);
Array<float> backward(Array <float> &, const Array <float> &, float);

void ramp_filter(Array <complex_t> &, float);
Array<complex_t> add_padding(const Array <float> &, float);
Array<float> remove_padding(const Array <complex_t> &, float);
Array<float> fftconvolve(const Array <float> &, const Array<float> &);
void normalize(Array<float> &, int);
void write_output(const char *, char *, size_t);
float calc_error(const Array<float> &, const Array<float> &);


#endif // TOMOCAM__H
