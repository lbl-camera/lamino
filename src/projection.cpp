#include <iostream>

#include <finufft.h>

#include "array.h"
#include "tomocam.h"
#include "dtypes.h"
#include "fft.h"
#include "kernel.h"



Array<float> backward(Array <float> & sinogram,
                      const Array <float> & angles, 
                      float center) {

    // add zero-padding
    Array<float> c = add_padding1d<float>(sinogram, 2.0);
    int cen_shift = (c.dims().y - sinogram.dims().y)/2;
    center += cen_shift;

    // dimensions
    dim2_t dims = c.dims();
    int num_projs = dims.x;
    int num_pixel = dims.y;

    // 1-D Fourier transforms
    fftshift1<float>(c);
    auto cz = c.cmplx();
    Array<complex_t> C = fft1(cz) / static_cast<complex_t>(num_pixel);
    xfftshift(C, -center); 

   // create points from [-pi, pi]
    int64_t M = num_projs * num_pixel;
    float * x = new float[M];
    float * y = new float[M];
    float cen = num_pixel / 2.f;
    float s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            float r = s * static_cast<float>(j - cen);
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        } 

    // nufft inputs
    int N1 = num_pixel;
    int N2 = num_pixel;
    int iflag = 1;
    float tol = 1.0e-06;
 
    nufft_opts *opts = new nufft_opts;
    finufftf_default_opts(opts);
    opts->upsampfac = 2.0;

    // working arrays
    Array<complex_t> F(num_pixel, num_pixel);
    
    // backward FFT
    int ier = finufftf2d1(M, x, y, C.ptr(), iflag, tol, N1, N2, F.ptr(), opts);

    delete [] x;
    delete [] y;
    delete [] opts;

    auto rv = F.real();
    
    return  remove_padding2d(rv, 2.0);
}

Array<float> forward(Array <float> & image, 
                    const Array <float> &angles, float center) {

    Array<float> F = add_padding2d(image, 2.);
    int cen_shift = (F.dims().y - image.dims().y)/2;
    center += cen_shift;
    int num_projs = angles.dims().x;
    int num_pixel = F.dims().y;
 
    // create points from [-pi, pi]
    int64_t M = num_pixel * num_projs;
    float * x = new float[M];
    float * y = new float[M];
    float cen = num_pixel / 2.0;
    float s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            float r = s * static_cast<float>(j - cen);
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        } 

    // nufft paramters 
    int N1 = num_pixel;
    int N2 = num_pixel;
    int iflag = -1;
    float tol = 1.0e-06;
    
    // cast float array to complex
    Array<complex_t> Ft = F.cmplx();
    Array<complex_t> C(num_projs, num_pixel);    

    // execute nufft
    nufft_opts *opts = new nufft_opts;
    finufftf_default_opts(opts);
    opts->upsampfac = 2.0;
    int ier = finufftf2d2(M, y, x, C.ptr(), iflag, tol, N1, N2, Ft.ptr(), opts);

    C = C / static_cast<complex_t>(N1 * N2);

    // axis shift
    xfftshift(C, center);

    // ifft
    Array<complex_t> Ct = ifft1(C);

    // shift
    fftshift1(Ct);

    delete [] x;
    delete [] y;
    delete [] opts;

    auto rv = Ct.real();
    return remove_padding1d(rv, 2.0);
}



