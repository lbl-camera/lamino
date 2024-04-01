#include <iostream>

#include <fftw3.h>
#include <finufft.h>

#include "dtypes.h"
#include "array.h"
#include "fft.h"
#include "tomocam.h"


Array<double> calc_psf(const Array<double> &angles, int num_pixel) {
    int num_projs = angles.size();

    complex_t v(1, 0);
    int N = 2 * num_pixel + 1;
    std::cout << N << std::endl;
    //int N = 2 * num_pixel;
    Array<complex_t> ones(num_projs, num_pixel, v);
    Array<complex_t> psf(N, N);

    int64_t M = num_pixel * num_projs;
    double *x = new double[M];
    double *y = new double[M];
    double cen = (num_pixel-1)/2.0;
    double s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            double r = s * static_cast<double>(j - cen);
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        }

    // nufft params
    int iflag = 1;
    double tol = 1.0E-12;

    // numff opts
    finufft_opts *opts = new finufft_opts;
    finufft_default_opts(opts);
    opts->upsampfac = 2.0;
    int ier = finufft2d1(M, x, y, ones.ptr(), iflag, tol, N, N, psf.ptr(), opts);

    delete []x;
    delete []y;
    return psf.real();
}



Array<double> fftconvolve(Array<double> signal, Array<double> filter) {
    // pad array
    auto x = addPadding(signal, filter.dims());

    // Forward FFT
    auto Xt = fft2r2c(x);
    auto Ft = fft2r2c(filter); 
    
    // multiply with the FT of signal
    Xt *= Ft;

    // Backward FFT
    auto z = fft2c2r(Xt)/static_cast<double>(filter.size());;
 
    // remove padding and return 
    return removePadding(z, signal.dims());
}

