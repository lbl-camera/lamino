#include <fstream>
#include <iostream>

#include "array.h"
#include "convolve.h"
#include "kernel.h"

const char * FILENAME = "/home/dkumar/data/shepp_logan/sino400.bin";
const int dims[] = { 400, 400 };
const int num_angles = 400;

int main(int argc, char **argv)
{

    // read data
    size_t size = num_angles * dims[1];
    float *data = new float[size];

    // read sinogram
    std::fstream fp(FILENAME, std::ios::in | std::ios::binary );
    if (!fp.is_open()) {
        std::cerr << "error! unable to open data file." << std::endl;
        exit(1);
    }
    fp.read((char *)data, sizeof(float) * size);
    fp.close();

    // angles
    Array < float >angles(1, num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / (num_angles - 1);

    // oversampling
    float oversample = 1.5;
    if (argc >= 2) {
        oversample = atof(argv[1]);
    }
    // convolution kernel radius
    float radius = 2.f;
    if (argc >= 3)
        radius = atof(argv[2]);
    // kaiser window scaling  (default 4 * PI)
    float beta = 12.566371;
    if (argc >= 4)
        beta = atof(argv[3]);

    std::cout << "Over-sampling : " << oversample << ", " << std::endl;
    std::cout << "Kernel radius : " << radius << ", " << std::endl;
    std::cout << "Kernel scaling: " << beta << "." << std::endl;
    Kernel kernel(radius, beta);

    // input and output with padding
    int ipad = (int)((oversample - 1.f) * dims[1] / 2);
    int padded = dims[1] + 2 * ipad;

    // working arrays
    Array < complex_t > sino(num_angles, padded);
    Array < complex_t > image(padded, padded);

    dim2_t idims = sino.dims();

	#pragma omp parallel for
    for (int i = 0; i < num_angles; i++)
        for (int j = 0; j < dims[1]; j++)
            sino(i, j + ipad) = complex_t(data[i * dims[1] + j], 0.f);

    float center = 200 + ipad;
    backward(sino, image, angles, kernel, center);

	// remove padding
	size = dims[0] * dims[1];
	float * recn = new float[size]; 
	remove_padding(recn, image, ipad, ipad);	

	// write output
	write_output(recn, size);

    return 0;
}
