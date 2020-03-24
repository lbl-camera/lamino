#include <iostream>
#include <fstream>
#include "array.h"

double calc_error(float *data, Array_ref<complex_t> array, int ipad) {
    dim3_t dims = array.dims();
    double error = 0.f;
    int DIM = dims.z - 2 * ipad;
    for (int j = 0; j < DIM; j++) 
        for (int k = 0; k < DIM; k++) {
            float v = array(0, j+ipad, k+ipad).real();
            error += (double) std::pow(v - data[j * DIM + k], 2);
        }
    return std::sqrt(error);
}

void write_output(Array_ref<complex_t> array, int ipad) {
    dim3_t dims = array.dims();
    int DIM = dims.z - 2 * ipad;
    float * image = new float[DIM * DIM];

    for (int j = 0; j < DIM; j++) 
        for (int k = 0; k < DIM; k++) 
            image[j * DIM + k] = array(0, j + ipad, k + ipad).real();

    std::ofstream real("real.bin", std::ios::out | std::ios::binary);
    size_t bytes = DIM * DIM * sizeof(float);
    real.write((char *) image, bytes);
    real.close();
}
