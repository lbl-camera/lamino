#include <fstream>
#include <iostream>

#include "array.h"
#include "tomocam.h"
#include "adam.h"

const int MAX_ITERS = 10;
const int num_angles = 360;
const int num_pixels = 400;
const char * infilename = "/home/dkumar/";

int main() {

    // read data
    size_t size = num_pixels * num_pixels;
    float *data = new float[size];
    float center = 200;
    // read sinogram
    Array<float> sinogram(num_angles, num_pixels);
    sinogram.fromfile(infilename);

    // angles
    Array <float> angles(1, num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / (num_angles - 1);

    // oversampling
    float oversample = 2.f;

    // transpose data
    auto sinoT = backward(sinogram, angles, center);
    sinoT.tofile("output.bin");
    std::exit(1);

    // calculate point spread function
    Array<float> ones(num_angles, num_pixels, 1);
    auto psf = backward(ones, angles, center);

    // turn the crank  
    msd::AdamMinimizer adam(num_pixels, num_pixels);

    std::cout << "starting .... " << std::endl;
    Array<float> recon(num_pixels, num_pixels, 1);
    for (int iter = 0; iter < MAX_ITERS; iter++) {
        auto xf = forward(recon, angles, center);
        auto g = backward(xf - sinogram, angles, center); 
        auto err = g.norm();
        std::cout << "step " << iter << ", error:" << err << std::endl;
        recon -= adam.update(g);
    }
    return 0;
}
