#include <fstream>
#include <iostream>
#include <ctime>
#include <array>

#include "array.h"
#include "tomocam.h"
#include "nesterov.h"
#include "fft.h"
#include "nufft.h"
#include "reader.h"
#include "writer.h"
#include "timer.h"

typedef std::complex<double> complex_t;
template <typename T>
T random() {
    auto a = static_cast<T>(rand());
    auto b = static_cast<T>(RAND_MAX);
    return (a / b);
}

int main() {

    // Open an hdf5 file
    tomocam::H5Writer writer("toeplitz.h5");

    int num_pixels = 2047;
    int num_angles = 360;

    // set seed
    srand(101);

    // angles
    std::vector<double> angles(num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / static_cast<double>(num_angles);

    Array<double> recon(num_pixels, num_pixels);
    for (int i = 0; i < recon.size(); i++)
        recon[i] = random<double>();

    auto pg = nufft::polar_grid<double>(angles, num_pixels);
    auto psf = calc_psf(pg);
    writer.write("psf", psf);

    tomocam::Timer t1;
    t1.start();
    auto g1 = gradient(recon, pg);
    t1.stop();
    writer.write("nufft", g1);

    // gradient 2
    auto x2 = recon;
    tomocam::Timer t2;
    t2.start();
    auto g2 = fftconvolve(x2, psf);
    t2.stop();
    writer.write("toeplitz", g2);

    std::cout << "NUFFT time (ms): " << t1.seconds() << " seconds" << std::endl;
    std::cout << "Toeplitz time (ms): " << t2.seconds() << " seconds" << std::endl;
    // compare
    auto err = norm(g1 - g2)/norm(g1);
    std::cout << "Error: " << err << std::endl;
    return 0;
}
