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

    dim2_t dims = sinogram.dims();
    int num_projs = dims.x;
    int num_pixel = dims.y;

    // 1-D Fourier transforms
    fftshift1<float>(sinogram);
    Array<complex_t> C = fft1(sinogram.cmplx());
    xfftshift<complex_t>(C, center); 

   // create points from [-pi, pi]
    int64_t M = num_projs * num_pixel;
    float * x = new float[M];
    float * y = new float[M];
    float s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            float r = s * static_cast<float>(j - center);
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        } 

    // nufft inputs
    int N1 = num_pixel;
    int N2 = num_pixel;
    int iflag = 1;
    float tol = 1e-06;
 
    nufft_opts *opts = new nufft_opts;
    finufftf_default_opts(opts);
    opts->upsampfac = 2.0;

    // working arrays
    Array<complex_t> F(num_pixel, num_pixel);
    
    // backward FFT
    int ier = finufftf2d1(M, y, x, C.ptr(), iflag, tol, N1, N2, F.ptr(), opts);

    delete [] x;
    delete [] y;
    delete [] opts;

    auto rv = F.real() / (float) (num_pixel * num_pixel);
    return  rv;
}

Array<float> forward(Array <float> & image, 
                    const Array <float> &angles, float center) {

    
    int num_projs = angles.dims().y;
    int num_pixel = image.dims().y;
    Array<complex_t> C(num_projs, num_pixel);    
 
    // create points from [-pi, pi]
    int64_t M = num_pixel * num_projs;
    float * x = new float[M];
    float * y = new float[M];
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            float r = 2 * M_PI * static_cast<float>(j - center) / num_pixel;
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        } 

    // nufft paramters 
    int N1 = num_pixel;
    int N2 = num_pixel;
    int iflag = -1;
    float tol = 1e-06;
    
    // cast float array to complex
    Array<complex_t> F = image.cmplx();

    // execute nufft
    nufft_opts *opts = new nufft_opts;
    finufftf_default_opts(opts);
    opts->upsampfac = 2.0;
    int ier = finufftf2d2(M, y, x, C.ptr(), iflag, tol, N1, N2, F.ptr(), opts);

    // axis shift
    xfftshift(F, -center);

    // ifft
    ifft1(C);

    // shift
    fftshift1(C);

    delete [] x;
    delete [] y;
    delete [] opts;

    return C.real();
}
