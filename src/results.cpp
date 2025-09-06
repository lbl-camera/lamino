#include <fstream>
#include <iostream>

#include "array.h"

double calc_error(const Array<double> &array, const Array<double> &data) {
    double err = 0.f;

    for (int i = 0; i < array.size(); i++)
        err += pow(array[i] - data[i], 2);
    return std::sqrt(err);
}

void write_output(char * filename, char *buffer, size_t bytes) {
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    out.write(buffer, bytes);
    out.close();
}

void ramp_filter(Array<complex<double>> & array, double fmax) {

    dim2_t dims = array.dims();
    double cen = static_cast<double>(dims.y) / 2.f;
    double a = static_cast<double>(dims.y) / M_PI;

    int jmin = static_cast<int>((1. - fmax) * dims.y) / 2;
    int jmax = dims.y - jmin;
    for (int i = 0; i < dims.x; i++) {
        for (int j = 0; j < dims.y; j++) {
            if (j < jmin) array(i, j) = 0.f;
            else if (j > jmax) array(i, j) = 0.f;
            else {
                double t = static_cast<double>(j - cen) / static_cast<double>(dims.y);
                double w = a * std::abs(std::sin(t * M_PI));
                array(i, j) *= w;
            }
        }
    }
}
