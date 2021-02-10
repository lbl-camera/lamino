#include <fstream>
#include <iostream>

#include "array.h"
#include "convolve.h"
#include "kernel.h"

const int dims[] = { 400, 400 };

const int num_angles = 360;
int main(int argc, char **argv)
{

    // read data
    size_t size = dims[0] * dims[1];
    float *data = new float[size];

    // read sinogram
    std::fstream fp("shepp_logan.bin", std::ifstream::in);
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
    Array < complex_t > image(padded, padded);
    Array < complex_t > sino(num_angles, padded);

    dim2_t idims = image.dims();
#pragma omp parallel for
    for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
            image(i + ipad, j + ipad) = data[i * dims[1] + j];

    float center = 200.f;
    forward(image, sino, angles, kernel, center);

    // reset zero-padding
    idims = sino.dims();
    for (int i = 0; i < idims.x; i++)
        for (int j = 0; j < idims.y; j++)
			if ((j < ipad) || (j >= idims.y - ipad))
				sino(i, j) = 0;

    image.clear();
    backward(sino, image, angles, kernel, center);

    float * output = new float [size];
    remove_padding(output, image, ipad, ipad);
    write_output(output, size);

    delete [] output;
    delete [] data;
    return 0;
}
