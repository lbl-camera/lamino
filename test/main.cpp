#include <fstream>
#include <iostream>

#include "array.h"
#include "tomocam.h"
#include "nesterov.h"
#include "fft.h"
#include "reader.h"

const int MAX_ITERS = 4;
//const char * FILENAME = "/home/dkumar/data/phantom_00016/phantom_00016/phantom_00016.h5";
const char * FILENAME = "/home/dkumar/data/shepp_logan/shepp_logan.h5";

int main() {

    // read data
    float center = 256;
    //float center = 640;

    // read sinogram
    tomocam::H5Reader reader(FILENAME);
    reader.setDataset("projs"); 
    Array<float> sinogram = reader.read_sinogram(0);
    sinogram /= sinogram.max();
    int num_pixels = sinogram.dims().y;
    int num_angles = sinogram.dims().x;

    // angles
    Array <float> angles = reader.read_angles("angs");
    if (num_angles != angles.dims().x) {
        std::cerr << "booo!" << std::endl;
        std::exit(1);
    }

    // transpose data
    auto sinoT = backward(sinogram, angles, center);

    // point spread function
    auto psf = calc_psf(angles, num_pixels);
   
    Array<float> recon(num_pixels, num_pixels, 1);

    // gradient 1
    auto t1 = forward(recon, angles, center) - sinogram;
    //t1.tofile("error.bin");
    
    auto g1 = backward(t1, angles, center);
    g1.tofile("gradient1.bin");

    // gradient 2
    auto g2 = fftconvolve(recon, psf) - sinoT;
    g2.tofile("gradient2.bin");
    return 0;
}
