#include <cmath>
#include <iostream>
#include <tr1/cmath>

#include "array.h"
#include "dtypes.h"
#include "kernel.h"

void deconvolve(Array < complex_t > &arr, Kernel kernel) {

    dim2_t dims = arr.dims();

	float cen = dims.y / 2;

    // bit hacky, but wtf
    bool twod = false;
    if (dims.x == dims.y)
        twod = true;

	if (twod) {
		#pragma omp parallel for
		for (int i = 0; i < dims.x; i++) {
			float ky = (i - cen) / static_cast<float>(dims.x);
			float wy = kernel.kaiserFT(ky);
			for (int j = 0; j < dims.y; j++) {
				float kx = (j - cen) / static_cast<float>(dims.y);
				float wx = kernel.kaiserFT(kx);
				arr(i, j) /= (wx * wy);
			}
		}
	} else {
		#pragma omp parallel for
		for (int i = 0; i < dims.x; i++) 
			for (int j = 0; j < dims.y; j++)  {
				float kx = (j - cen) / static_cast<float>(dims.y);
				arr(i, j) /= kernel.kaiserFT(kx);
			}
	}
}
