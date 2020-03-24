

#include <iostream>

#include "array.h"
#include "dtypes.h"
#include "kernel.h"

void convolve_p2c(Array_ref<complex_t> input, Array_ref<complex_t> output, Array_ref<float> angles, Kernel kernel) {

    dim3_t idims = input.dims();
    dim3_t odims = output.dims();
    float center = (float) (idims.z / 2);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < idims.x; i++) {
        for (int j = 0; j < idims.y; j++) {
            float angle = angles[j];
            float cosa = std::cos(angle);
            float sina = std::sin(angle);

            for (int k = 0; k < idims.z; k++){
                float y = (k - center) * sina + center;
                float z = (k - center) * cosa + center;
                complex_t val = input(i, j, k);

                float r = kernel.radius();
                int iy = std::max((int) (y - r), 0);
                int iymax = std::min((int) (y + r), odims.y-1);
                int izmin = std::max((int) (z - r), 0);
                int izmax = std::min((int) (z + r), odims.z-1);
 
                for (; iy < iymax; iy++) {
                    complex_t temp = val * kernel.weight(y - iy);
                    for (int iz = izmin; iz < izmax; iz++) {
                        #pragma omp critical
                        output(i, iy, iz) += temp * kernel.weight(z - iz);
                    } // iz
                } // iy
            } // k
        } // j
    } // i
}

void convolve_c2p(Array_ref<complex_t> input, Array_ref<complex_t> output, Array_ref<float> angles, Kernel kernel) {

    dim3_t idims = input.dims();
    dim3_t odims = output.dims();
    float center = (float) (idims.z / 2);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < idims.x; i++) {
        for (int j = 0; j < idims.y; j++) {
            float angle = angles[j];
            float cosa = std::cos(angle);
            float sina = std::sin(angle);

            for (int k = 0; k < idims.z; k++){
                float y = (k - center) * sina + center;
                float z = (k - center) * cosa + center;

                float r = kernel.radius();
                int iy = std::max((int) (y - r), 0);
                int iymax = std::min((int) (y + r), odims.y-1);
                int izmin = std::max((int) (z - r), 0);
                int izmax = std::min((int) (z + r), odims.z-1);
 
                for (; iy < iymax; iy++) {
                    float wy = kernel.weight(y - iy);
                    for (int iz = izmin; iz < izmax; iz++) {
                        float wz = kernel.weight(z - iz);
                        #pragma omp critical
                        output(i, j, k) += input(i, iy, iz) * wy * wz;
                    } // ix
                } // iy
            } // k
        } // j
    } // i
}
