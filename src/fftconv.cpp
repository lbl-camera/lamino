#include <iostream>

#include <fftw3.h>
#include <finufft.h>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "fft.h"
#include "nufft.h"
#include "padding.h"
#include "tomocam.h"

template <typename Real_t>
Array<Real_t> calc_psf(const PolarGrid<Real_t> &pg) {

    // get the dimensions
    int nprojs = pg.nprojs;
    int npixel = pg.npixel;
    int npts = nprojs * npixel;

    complex<Real_t> v(1, 0);
    Array<complex<Real_t>> ones(nprojs, npixel);
    for (int i = 0; i < npts; i++) ones[i] = v;

    // create psf array
    const int N = 2 * npixel - 1;
    Array<complex<Real_t>> psf(N, N);

    // type-1 NUFFT
    nufft::nufft2d1(ones, psf, pg);
    auto psf2 = real<Real_t>(psf);
    return ifftshift2(psf2);
}
template Array<double> calc_psf<double>(const PolarGrid<double> &);
template Array<float> calc_psf<float>(const PolarGrid<float> &);

template <typename Real_t>
Array<Real_t> fftconvolve(const Array<Real_t> &signal, const Array<Real_t> &filter) {

    // pad array
    int filter_size = filter.ncols();
    int signal_size = signal.ncols();

    //int s = filter_size + signal_size - 1;
    int s = filter_size;

    int padx = s - signal_size;
    auto x = pad2d<Real_t>(signal, padx, PadType::SYMMETRIC);
    x = ifftshift2<Real_t>(x);

    // Forward FFT
    auto Xt = fft2<Real_t>(to_complex<Real_t>(x));
    auto Ft = fft2<Real_t>(to_complex<Real_t>(filter));

    // multiply with the FT of signal
    Xt *= Ft;

    // Backward FFT
    auto z = ifft2<Real_t>(Xt);
    auto z2 = real<Real_t>(z) / static_cast<Real_t>(s * s);
    z2 = fftshift2<Real_t>(z2);
    z2 = crop2d<Real_t>(z2, padx, PadType::SYMMETRIC);

    // remove padding and return
    Real_t scale = (Real_t) std::pow(signal_size, 3);
    return z2 / scale;
}

template Array<double> fftconvolve(const Array<double> &, const Array<double> &);
template Array<float> fftconvolve(const Array<float> &, const Array<float> &);
