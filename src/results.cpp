#include "array.h"
#include <fstream>
#include <iostream>


void remove_padding(float *data, Array < complex_t > &array, int ipad, int jpad) {
    dim2_t dims = array.dims();
    int nrow = dims.x - 2 * ipad;
    int ncol = dims.y - 2 * jpad;

    for (int i = 0; i < nrow; i++)  
        for (int j = 0; j < ncol; j++) {
			int id = i *  ncol + j;
			data[id] = array(i + ipad, j + jpad).real();
		}
}

void normalize(float * array, int size) {
    float MAX = 1.0E-10;
    for (int i = 0; i < size; i++)
        if (array[i] > MAX) MAX = array[i];

    std::cout << "dividing by " << MAX << std::endl;
    #pragma omp parallel for   
    for (int i = 0; i < size; i++)
        array[i] /= MAX;
}

float calc_error(float * array, float * data, int size) {
    float err = 0.f;
    for (int i = 0; i < size; i++)
        err += pow(array[i] - data[i], 2);
    return std::sqrt(err);
}

void write_output(float *array, int size) {
    size_t bytes = size * sizeof(float);
    std::ofstream real("output.bin", std::ios::out | std::ios::binary);
    real.write((char *) array, bytes);
    real.close();
}


void ramp_filter(Array<complex_t> & array, float fmax) {

	dim2_t dims = array.dims();
	float cen = static_cast<float>(dims.y) / 2.f;
	float a = static_cast<float>(dims.y) / M_PI;	

	int jmin = static_cast<int>((1.-fmax) * dims.y) / 2;
	int jmax = dims.y - jmin;
	for (int i = 0; i < dims.x; i++) {
		for (int j = 0; j < dims.y; j++) {
			if (j < jmin) array(i, j) = 0.f;
			else if (j > jmax) array(i,j) = 0.f;
			else {
				float t = static_cast<float>(j-cen) / static_cast<float>(dims.y);
				float w = a * std::abs(std::sin(t * M_PI));
				array(i, j) *= w;
			}
		}
	}
}
