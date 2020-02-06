

#include <cmath>
#include <tr1/cmath>
#include <iostream>

#include "dtypes.h"

const double W = 5;
const double beta = 3 * M_PI;

double kaiser(double r) {
    double RMAX = 0.5 * (W - 1);
    double weight = 0.;
    if (std::abs(r) <  RMAX) {
        double t0 = 2 * r / W;
        double t1 = std::tr1::cyl_bessel_i(0, beta * std::sqrt(1 - t0 * t0));
        double t2 = std::tr1::cyl_bessel_i(0, beta);
        weight = t1/t2;
    }
    return weight;
}

void convolve(complex_t *input, complex_t *output, dim3_t idims, dim3_t odims, float *angles) {

    double center = (double) (idims.z / 2);

    #pragma omp parallel for
    for (int j = 0; j < idims.y; j++) {
        double angle = angles[j];
        double cosa = std::cos(angle);
        double sina = std::sin(angle);
        for (int k = 0; k < idims.z; k++){
            // global id (input)
            int ii = j * idims.z + k;

            double x = (k - center) * cosa + center;
            double y = (k - center) * sina + center;
            complex_t val = input[ii];

            int iy = std::max((int) y - 2, 0);
            int ixmin = std::max((int) x - 2, 0);
            int iymax = std::min((int) y + 3, odims.y-1);
            int ixmax = std::min((int) x + 3, odims.z-1);
 
            for (; iy < iymax; iy++) {
                double wy = kaiser(y - iy);
                if (wy < 1.0E-10) continue;
                complex_t temp = val * wy;
                for (int ix = ixmin; ix < ixmax; ix++) {
                    double wx = kaiser(x - ix);
                    if (wx < 1.0E-10) continue;
                    // global id (output)
                    int jj = iy * odims.z + ix;
                    #pragma omp critical
                    output[jj] += temp * wx;
                }
            }
        }
    }
}
