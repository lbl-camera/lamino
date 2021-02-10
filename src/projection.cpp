#include <iostream>

#include "array.h"
#include "convolve.h"
#include "dtypes.h"
#include "fft.h"
#include "kernel.h"

void backward(Array < complex_t > & input,
              Array < complex_t > & output,
              Array < float > & angles, Kernel kernel, float center)
{

    // fftshift
    fftshift1(input);

    // 1-D fft
    fft1(input);

	// put zero-frequency in the center
    //fftshift1(input);

    // correct center offset
    xfftshift(input, center);

	// apply ramp-filter
    ramp_filter(input, 0.65);

    // rescale
    float nn = std::pow(input.dims().y, 2);
#pragma omp parallel for
    for (int i = 0; i < input.size(); i++)
        input[i] /= nn;

    // convolution
    convolve_p2c(input, output, angles, kernel, center);

    // fftshift
    fftshift2(output);

    // ifft2
    ifft2(output);

    // fftshift2d
    fftshift2(output);

    // deconvolve
    deconvolve(output, kernel);
}

void forward(Array < complex_t > & input,
             Array < complex_t > & output,
             Array < float > &angles, Kernel kernel, float center)
{

    // 2D fft-shift (chess-board)
    fftshift2(input);

    // fft2 on oversampled grid
    fft2(input);

    // shift, again
    fftshift2(input);

    // rescale
    float nn = std::pow(input.dims().y, 2);
#pragma omp parallel for
    for (int i = 0; i < input.size(); i++)
        input[i] /= nn;

    // NUFFT convolution
    convolve_c2p(input, output, angles, kernel, center);

    // 1D shift
    //xfftshift(output, center);
    fftshift1(output);

    // ifft
    ifft1(output);

    // 1D shift
    fftshift1(output);

    // deconvolve
    deconvolve(output, kernel);
}
