#include <iostream>
#include <fstream>

#include "dtypes.h"
#include "fft.h"
#include "convolve.h"

const int dims[] = { 1, 1536, 1280 };

void write(complex_t * arr, dim3_t d) {
    std::ofstream out("slice.out", std::ofstream::out);
    for (int i = 0; i < d.y; ++i)
        for (int j = 0; j < d.z; ++j) {
            int id =  i * d.z + j;
            out << arr[id].real() << " ";
        }
    out.close();
}

int main(int argc, char **argv) {

    // read data
    size_t size = dims[0] * dims[1] * dims[2];
    float * data = new float[size];
    float * angs = new float[dims[1]];

    // read sinogram
    std::fstream fp("sino_00016.bin", std::ifstream::in);
    if (! fp.is_open()) {
        std::cerr << "error! unable to open data file." << std::endl;
        exit(1);
    }
    fp.read((char *) data, sizeof(float) * size);
    fp.close();

    // read angles
    std::fstream fp2("angles.bin", std::ifstream::in);
    if (! fp2.is_open() ) {
        std::cerr << "error! unable to open angles file." << std::endl;
        exit(1);
    }
    fp2.read((char *) angs, sizeof(float) * dims[1]);
    fp2.close();

    double center = 640;
    double oversample = 1.5;
    if (argc == 3) {
        center = atof(argv[1]);
        oversample = atof(argv[2]);
    }

    // input and output with padding
    int padded = (int) (oversample * dims[2]);
    int ipad = (padded - dims[2])/2;
    dim3_t idims(dims[0], dims[1], padded);
    dim3_t odims(dims[0], padded, padded);
    complex_t *input = new complex_t [idims.x * idims.y * idims.z];
    complex_t *output = new complex_t [odims.x * odims.y * odims.z];
    
    int cols = dims[2];
    #pragma omp parallel for
    for (int i = 0; i < idims.x * idims.y; i++)
        for (int j = 0; j < dims[2]; j++)
            input[i*idims.z + j] = complex_t((double) data[i*cols + j], 0);

    // fftshift1D
    fftshift1d(input, idims);

    // 1-D fft
    fft1d((fftw_complex *) input, idims);

    // center-shift
    xfftshift(input, idims, -center);

    // convolution
    convolve(input, output, idims, odims, angs);

    // fftshift2d
    fftshift2d(output, odims);

    // ifft2
    ifft2d((fftw_complex *) output, odims);

    // fftshift2d
    fftshift2d(output, odims);

    // de-apodize

    // write out
    write(output, odims);

    return 0;
}
