

#include "dtypes.h"


#ifndef CONVOLVE__H
#define CONVOLVE__H

void deapod(complex_t *, dim3_t);
void convolve(complex_t *, complex_t *, dim3_t , dim3_t, float *);


#endif // CONVOLVE__H
