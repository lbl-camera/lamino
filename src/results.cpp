#include <fstream>
#include <iostream>

#include "array.h"

float calc_error(const Array<float> &array, const Array<float> &data) {
    float err = 0.f;

    for (int i = 0; i < array.size(); i++)
        err += pow(array[i] - data[i], 2);
    return std::sqrt(err);
}

void write_output(char * filename, char *buffer, size_t bytes) {
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    out.write(buffer, bytes);
    out.close();
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
