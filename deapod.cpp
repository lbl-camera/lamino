#include <cmath>
#include <tr1/cmath>
#include <iostream>

#include "dtypes.h"

double kaiserFT(double W, double beta, double x) {
    double w = 1;
    double t1 = (beta * beta) - (x * x * W * W * M_PI * M_PI);
    if (t1 > 0) {
        t1 = std::sqrt(t1);
        double t2 = W / std::tr1::cyl_bessel_i(0, beta);
        w = t2 * std::sinh(t1)/t1;
    }
    return w;
} 


void deapod(complex_t *arr, dim3_t dims) {

    double W = 5.;
    double beta = 12.5663706144;

    double cen = dims.z / 2;
    for (int i = 0; i < dims.x; i++) {
        for (int j = 0; j < dims.y; j++) {
            double y = (j - cen)/dims.y;
            double wy = kaiserFT(W, beta, y);
            for (int k = 0; k < dims.z; k++) {
                double x = (k - cen)/dims.z;
                int idx = i * dims.y * dims.z + j * dims.z + k;
                arr[idx] /= (wy * kaiserFT(W, beta, x));
            }
        }
    }
}
