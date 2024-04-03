#include <fstream>
#include <iostream>

#include "array.h"
#include "tomocam.h"
#include "nesterov.h"

const int MAX_ITERS = 10;
const int num_angles = 400;
const int num_pixels = 400;
const char * FILENAME = "/home/dkumar/data/shepp_logan/shepp400.bin";

const float guass[] = {0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625};
int main() {

    // read sinogram
    Array<float> signal(num_pixels, num_pixels);
    signal.fromfile(FILENAME);

    // filter
    Array<float> filter(num_pixels, num_pixels);
    int i0 = num_pixels/2 - 1;
    int j0 = num_pixels/2 - 1;
    
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            filter(i0+i, j0+j) = guass[i * 3 + j]; 

    // calculate grdient with Toeplitz method
    auto g2 = fftconvolve(signal, filter);
    
    g2.tofile("output.bin");
    return 0;
}
