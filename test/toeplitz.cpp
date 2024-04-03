#include <fstream>
#include <iostream>
#include <ctime>
#include <array>

#include "array.h"
#include "tomocam.h"
#include "nesterov.h"
#include "finufft.h"
#include "fft.h"
#include "reader.h"
#include "writer.h"

const int MAX_ITERS = 4;
const char * FILENAME = "/home/dkumar/data/shepp_logan/shepp_logan.h5";


template <typename T>
T random() {
    auto a = static_cast<T>(rand());
    auto b = static_cast<T>(RAND_MAX);
    return (a/b);
}

int main() {

    // read data
    double center = 256;

    // read sinogram
    tomocam::H5Reader reader(FILENAME);
    tomocam::H5Writer writer("toeplitz.h5");

    reader.setDataset("projs"); 
    auto sinogram = reader.read_sinogram(0);

    int num_pixels = sinogram.ncols();
    int num_angles = sinogram.nrows();

    std::cout << "Num of projs: " << num_angles << std::endl;

    // angles
    auto anglesf = reader.read_angles("angs");
    if (num_angles != anglesf.dims().x) {
        std::cerr << "booo!" << std::endl;
        std::exit(1);
    }
    Array<double> angles(anglesf.nrows(), anglesf.ncols());
    for (int i = 0; i < anglesf.size(); i++)
        angles[i] = anglesf[i];

    Array<double> recon(num_pixels, num_pixels);
    for (int i = 0; i < recon.size(); i++)
        recon[i] = random<double>();

    auto psf = calc_psf(angles, num_pixels);
    writer.write("psf", psf);
   
    const int M = num_angles * num_pixels;
    int N = num_pixels;
    double * x = new double[M];
    double * y = new double[M];
    double cen = num_pixels / 2.f;
    double s = 2 * M_PI / num_pixels;
    for (int i = 0; i < num_angles; i++)
        for (int j = 0; j < num_pixels; j++) {
            double r = s * static_cast<double>(j - cen);
            x[i * num_pixels + j] = r * std::cos(angles[i]);
            y[i * num_pixels + j] = r * std::sin(angles[i]);
        } 

    finufft_opts *opts = new finufft_opts;
    finufftf_default_opts(opts);
    opts->upsampfac = 2.0;
    double tol = 6.0E-12; 
    complex_t *C = new complex_t[M];
    auto Ft = recon.cmplx();
    int ier = finufft2d2(M, x, y,  C, -1, tol, N, N, Ft.ptr(), opts);
    ier = finufft2d1(M, x, y, C, 1, tol, N, N, Ft.ptr(), opts);
    auto g1 = Ft.real();
    writer.write("nufft", g1);
    delete []x;
    delete []y;
    delete []C;

    // gradient 2
    auto x2 = recon;
    auto g2 = fftconvolve(x2, psf);
    writer.write("toeplitz", g2);
    
    // compare
    auto err = (g1 - g2).norm();
    std::cout << "Error: " << err/g1.norm() << std::endl;
    return 0;
}
