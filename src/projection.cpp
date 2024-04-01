#include <iostream>

#include <finufft.h>

#include "array.h"
#include "tomocam.h"
#include "dtypes.h"
#include "fft.h"


constexpr double tol = 1.0E-12;

Array<double> backward(Array <double> & sinogram,
                      const Array <double> & angles, 
                      double center) {

    // dimensions
    dim2_t dims = sinogram.dims();
    int num_projs = dims.x;
    int num_pixel = dims.y;
    double cen = (num_pixel-1) / 2.0;

    // 1-D Fourier transforms
    auto cz = sinogram.cmplx();
    cz  = ifftshift(cz);
    Array<complex_t> C = fft1(cz);
    //C = C / static_cast<complex_t>(num_pixel);
    C = fftshift(C);

   // create points from [-pi, pi]
    int64_t M = num_projs * num_pixel;
    double * x = new double[M];
    double * y = new double[M];
    double s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++) {
        for (int j = 0; j < num_pixel; j++) {
            double r = s * (j-cen);
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        } 
    }
    // nufft inputs
    int N1 = num_pixel;
    int N2 = num_pixel;
    int iflag = 1;
 
    finufft_opts *opts = new finufft_opts;
    finufft_default_opts(opts);
    opts->upsampfac = 2.0;

    // working arrays
    Array<complex_t> F(num_pixel, num_pixel);
    
    // backward FFT
    int ier = finufft2d1(M, x, y, C.ptr(), iflag, tol, N1, N2, F.ptr(), opts);
    
    delete [] x;
    delete [] y;
    delete [] opts;

    return F.real();
}

Array<double> forward(Array <double> & image, 
                    const Array <double> &angles, double center) {

    auto F = image;
    int num_projs = angles.dims().x;
    int num_pixel = F.dims().y;
    double cen = (num_pixel-1)/2.0;
 
    // create points from [-pi, pi]
    int64_t M = num_pixel * num_projs;
    double * x = new double[M];
    double * y = new double[M];
    double s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            double r = s * (j - cen);
            x[i * num_pixel + j] = r * std::cos(angles[i]);
            y[i * num_pixel + j] = r * std::sin(angles[i]);
        } 

    // nufft paramters 
    int N1 = num_pixel;
    int N2 = num_pixel;
    int iflag = -1;
    
    // cast double array to complex
    Array<complex_t> Ft = F.cmplx();
    Array<complex_t> C(num_projs, num_pixel);    

    // execute nufft
    finufft_opts *opts = new finufft_opts;
    finufft_default_opts(opts);
    opts->upsampfac = 2.0;
    opts->modeord = 0;
    int ier = finufft2d2(M, x, y, C.ptr(), iflag, tol, N1, N2, Ft.ptr(), opts);
    
    // axis shift
    C = ifftshift(C); 

    // ifft
    Array<complex_t> Ct = ifft1(C) / static_cast<complex_t>(num_pixel);

    // shift
    Ct = fftshift(Ct);

    delete [] x;
    delete [] y;
    delete [] opts;

    return Ct.real();
}
