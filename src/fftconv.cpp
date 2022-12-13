#include <iostream>

#include <fftw3.h>
#include <finufft.h>

#include "dtypes.h"
#include "array.h"
#include "fft.h"
#include "tomocam.h"


Array<float> calc_psf(const Array<float> &angles, int num_pixel) {
    int num_projs = angles.dims().x;

    complex_t v(1, 0);
    Array<complex_t> ones(num_projs, num_pixel, v);
    Array<complex_t> psf(2*num_pixel, 2*num_pixel);

    int64_t M = num_pixel * num_projs;
    float *x = new float[M];
    float *y = new float[M];
    float cen = num_pixel/2.0;
    float s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            float r = s * static_cast<float>(j - cen);
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        }

    // nufft params
    int N1 = 2*num_pixel;
    int N2 = 2*num_pixel;
    int iflag = 1;
    float tol = 1.0E-06;

    // numff opts
    nufft_opts *opts = new nufft_opts;
    finufftf_default_opts(opts);
    opts->upsampfac = 2.0;
    int ier = finufftf2d1(M, x, y, ones.ptr(), iflag, tol, N1, N2, psf.ptr(), opts);
    return psf.real();
}



Array<float> fftconvolve(Array<float> signal, Array<float> filter) {
    // pad array
    int oversample = 2;
    auto x = add_padding2d(signal, oversample);

    // Forward FFT
    auto Xt = fft2r2c(x) / static_cast<complex_t>(x.size());
    auto Ft = fft2r2c(filter) / static_cast<complex_t>(filter.size());

    // multiply with the FT of signal
    Xt *= Ft;

    // Backward FFT
    auto z = fft2c2r(Xt);
 
    // remove padding and return 
    return remove_padding2d(z, oversample);
}

