#include <fstream>
#include <iostream>

#include "array.h"
#include "tomocam.h"
#include "nesterov.h"

const int MAX_ITERS = 1000;
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
    //auto sinoT = backward(sinogram, angles, center);

    // calculate point spread function
    //Array<float> ones(num_angles, num_pixels, 1);
    //auto psf = backward(ones, angles, center);

    // turn the crank  
    opt::NAGoptimizer opt(num_angles, num_pixels, angles, center);

    std::cout << "starting .... " << std::endl;
    Array<float> recon(num_pixels, num_pixels, 1);
    for (int iter = 0; iter < MAX_ITERS; iter++) {
        auto xf = forward(recon, angles, center);
        auto e = xf - sinogram;
        auto g = backward(e, angles, center); 
        auto err = g.norm();
        std::cout << "step " << iter << ", error:" << err << std::endl;
        opt.update(recon, g);
    }
    recon.tofile("output.bin");
    return 0;
}
