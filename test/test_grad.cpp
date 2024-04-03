#include <fstream>
#include <iostream>

#include "array.h"
#include "tomocam.h"
#include "nesterov.h"

const int MAX_ITERS = 11;
const int num_angles = 400;
const int num_pixels = 400;
const char * FILENAME = "/home/dkumar/data/shepp_logan/sino400.bin";

int main() {

    // read data
    size_t size = num_pixels * num_pixels;
    float *data = new float[size];
    float center = 200;
    // read sinogram
    Array<float> sinogram(num_angles, num_pixels);
    sinogram.fromfile(FILENAME);
    sinogram /= sinogram.max();

    // angles
    Array <float> angles(1, num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / (num_angles - 1);

    // oversampling
    float oversample = 2.f;

    // transpose data
    auto sinoT = backward(sinogram, angles, center);

    // calculate point spread function
    auto psf = calc_psf(num_pixels, num_angles, angles, center);

    // initalize recon array
    Array<float> recon(num_pixels, num_pixels, 1);

    // calulate gradient normally
    auto t1 = forward(recon, angles, center);
    auto g1 = backward(t1-sinogram, angles, center);   
 
    // calculate grdient with Toeplitz method
    auto g2 = fftconvolve(recon, psf) - sinoT;
    
    return 0;
}
