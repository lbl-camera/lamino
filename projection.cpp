
#include "dtypes.h"
#include "array.h"
#include "kernel.h"
#include "fft.h"
#include "convolve.h"

void backward(Array_ref<complex_t> input, Array_ref<complex_t> output, Array_ref<float> angles, Kernel kernel) {

    // fftshift
    fftshift1(input);

    // 1-D fft
    fft1(input);

// rescale
    float nn = std::pow(input.dims().z, 2);
#pragma omp parallel for
    for (int i = 0; i < input.size(); i++) 
        input[i] /= nn;

    // convolution
    convolve_p2c(input, output, angles, kernel);

    // fftshift
    fftshift2(output);

    // ifft2
    ifft2(output);

    // fftshift2d
    fftshift2(output);

    // deconvolve
    deconvolve(output, kernel);
}

void forward(Array_ref<complex_t> input, Array_ref<complex_t> output, Array_ref<float> angles, Kernel kernel) {

    // 2D fft-shift (chess-board)
    fftshift2(input);

    // fft2 on oversampled grid
    fft2(input);

    // shift, again
    fftshift2(input);

    // rescale
    float nn = std::pow(input.dims().z, 2);
#pragma omp parallel for
    for (int i = 0; i < input.size(); i++) input[i] /= nn;

    // NUFFT convolution
    convolve_c2p(input, output, angles, kernel);

    // 1D shift
    fftshift1(output);

    // ifft
    ifft1(output);

    // deconvolve
    deconvolve(output, kernel);
}
