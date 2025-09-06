#include <fstream>
#include <iostream>

#include "array.h"
#include "tomocam.h"
#include "nesterov.h"
#include "utils.h"

const int num_angles = 400;
const int num_pixels = 400;

const float guass[] = {0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625};
int main() {

    // read sinogram
    auto rng = cam::NumPyRandom();
    Array<float> signal(num_pixels, num_pixels);

    #pragma omp parallel for
    for (int i = 0; i < signal.size(); i++)
        signal(i) = rng.rand<float>();

    // filter
    Array<float> filter(num_pixels, num_pixels);
    int j0 = num_pixels / 2 - 1;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            filter(i0 + i, j0 + j) = guass[i * 3 + j];

    // calculate grdient with Toeplitz method
    auto g2 = fftconvolve(signal, filter);

    return 0;
}
