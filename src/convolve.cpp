
#include <iostream>

#include "array.h"
#include "dtypes.h"
#include "kernel.h"

void convolve_p2c(Array<complex_t> & input, Array<complex_t> & output,
    Array<float> & angles, Kernel kernel, float center) {

    dim2_t idims = input.dims();
    dim2_t odims = output.dims();
    float r = kernel.radius();

    #pragma omp parallel for 
    for (int i = 0; i < idims.x; i++) {
        float angle = angles[i];
        float cosa = std::cos(angle);
        float sina = std::sin(angle);

        for (int j = 0; j < idims.y; j++) {
            float x = (j - center) * cosa + center;
            float y = (j - center) * sina + center;
            complex_t val = input(i, j);

            int iy = (int) std::max(std::floor(y - r), 0.f);
            int iymax = (int) std::min(std::ceil(y + r), (float) (odims.x - 1));

            int ixmin = (int) std::max(std::floor(x - r), 0.f);
            int ixmax = (int) std::min(std::ceil(x + r), (float) (odims.y - 1));

            for (; iy < iymax; iy++) {
                complex_t temp = val * kernel.weight(y - iy);

                for (int ix = ixmin; ix < ixmax; ix++) {
                    #pragma omp critical
                        output(iy, ix) += temp * kernel.weight(x - ix);
                }
            }
        }
    }
}

void convolve_c2p(Array<complex_t> & input, Array<complex_t> & output,
    Array<float> & angles, Kernel kernel, float center) {

    dim2_t idims = input.dims();
    dim2_t odims = output.dims();
    float r = kernel.radius();

    #pragma omp parallel for
    for (int i = 0; i < odims.x; i++) {
        float angle = angles[i];
        float cosa = std::cos(angle);
        float sina = std::sin(angle);

        for (int j = 0; j < odims.y; j++) {
            float x = (j - center) * cosa + center;
            float y = (j - center) * sina + center;

            int iy = (int) std::max(std::floor(y - r), 0.f);
            int iymax = (int) std::min(std::ceil(y + r), (float) (idims.x - 1));
            int ixmin = (int) std::max(std::floor(x - r), 0.f);
            int ixmax = (int) std::min(std::ceil(x + r), (float) (idims.y - 1));

            for (; iy < iymax; iy++) {
                float wy = kernel.weight(y - iy);

                for (int ix = ixmin; ix < ixmax; ix++)
                    #pragma omp critical
                    output(i, j) += input(iy, ix) * wy * kernel.weight(x - ix);
            }
        }
    }
}
