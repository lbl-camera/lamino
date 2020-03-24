#include <cmath>
#include <tr1/cmath>
#include <iostream>

#include "dtypes.h"
#include "array.h"
#include "kernel.h"

void deconvolve(Array<complex_t> & arr, Kernel kernel) {

    dim3_t dims = arr.dims();

    // bit hacky, but wtf
    bool twod = false;
    if (dims.y == dims.z) 
        twod = true;

    float cen = dims.z / 2;
    float wy, wz;
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++) {
            if (twod)
                wy = kernel.kaiserFT((j - cen) / dims.y);
            else 
                wy = 1.f; 
            for (int k = 0; k < dims.z; k++) {
                wz = kernel.kaiserFT((k - cen) / dims.z);
                arr(i,j,k) /= wy * wz;
            }
        }
}
