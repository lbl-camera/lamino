

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
        double t1 = std::tr1::cyl_bessel_i(0, beta * std::sqrt(1 - std::pow(2*r/W, 2)));
        double t2 = std::tr1::cyl_bessel_i(0, beta);
        weight = t1/t2;
    }
    return weight;
}

coord_t<double> cartesian(int ix, double angle, double center) {
    coord_t<double> r;
    r.x = (double) ((ix - center) * std::cos(angle) + center);
    r.y = (double) ((ix - center) * std::sin(angle) + center);
    return r;
}

void convolve(complex_t * input, complex_t * output, dim3_t idims, dim3_t odims) {

    double center = (double) idims.z / 2.;
    double dtheta = M_PI / (double) idims.y ;

    #pragma omp parallel for
    for (int j = 0; j < idims.y; j++) {
        double angle = j * dtheta;
        for (int k = 0; k < idims.z; k++){
            coord_t<double> r = cartesian(k, angle, center);
            int ii = j * idims.z + k;
            complex_t val = input[ii];

            int iy = std::max((int) r.y - 3, 0);
            int ixmin = std::max((int) r.x - 3, 0);
            int iymax = std::min((int) r.y + 3, odims.y-1);
            int ixmax = std::min((int) r.x + 3, odims.z-1);
 
            for (; iy < iymax; iy++) {
                double wy = kaiser(r.y - iy);
                if (wy < 1.0E-10) continue;
                complex_t temp = val * wy;
                for (int ix = ixmin; ix < ixmax; ix++) {
                    double wx = kaiser(r.x - ix);
                    if (wx < 1.0E-10) continue;
                    int jj = iy * odims.z + ix;
                    #pragma omp critical
                    output[jj] += temp * wx;
                }
            }
        }
    }
}
