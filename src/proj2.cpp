#include <iostream>
#include <finufft.h>

#include "array.h"
#include "tomocam.h"
#include "dtypes.h"
#include "fft.h"


constexpr double tol = 1.0E-12;

Array<double> gradient(const Array<double> & f, 
                       const Array<double> &theta, double center) {

    int num_projs = theta.size();
    int num_pixel = f.ncols();

    // create points from [-pi, pi]
    double cen = num_pixel / 2;
    int64_t M = num_pixel * num_projs;
    double * x = new double[M];
    double * y = new double[M];
    double s = 2 * M_PI / num_pixel;
    for (int i = 0; i < num_projs; i++)
        for (int j = 0; j < num_pixel; j++) {
            double r = s * (j - cen);
            x[i * num_pixel + j] = r * std::cos(theta[i]);
            y[i * num_pixel + j] = r * std::sin(theta[i]);
        } 

    // nufft paramters 
    int N1 = num_pixel;
    int N2 = num_pixel;
 
    // cast to complex
    Array<complex_t> fz = f.cmplx();
    Array<complex_t> cz(num_projs, num_pixel);

    // execute nufft
    finufft_opts *opts = new finufft_opts;
    finufft_default_opts(opts);
    opts->upsampfac = 2.0;
    opts->modeord = 0;
    int ier = finufft2d2(M, x, y, cz.ptr(), -1, tol, N1, N2, fz.ptr(), opts);
    ier = finufft2d1(M, x, y, cz.ptr(), 1, tol, N1, N2, fz.ptr(), opts);
 
    delete [] x;
    delete [] y;
    delete [] opts;

    return fz.real();
}



