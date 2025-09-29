#include <fstream>
#include <iostream>
#include <ctime>

#include "array.h"
#include "array_ops.h"
#include "tomocam.h"
#include "fft.h"
#include "nufft.h"
#include "reader.h"
#include "writer.h"

const int MAX_ITERS = 4;
const char * FILENAME = "/home/dkumar/data/shepp_logan/shepp_logan.h5";

template <typename T>
T random() {
    auto a = static_cast<T>(rand());
    auto b = static_cast<T>(RAND_MAX);
    return (a / b);
}

int main() {

    // read data
    double center = 256;

    // read sinogram
    tomocam::H5Reader reader(FILENAME);
    tomocam::H5Writer writer("twocode.h5");

    reader.setDataset("projs");
    auto data = reader.read_sinogram(0);
    Array<double> sinogram(data.nrows(), data.ncols());
    for (int i = 0; i < data.size(); i++)
        sinogram[i] = static_cast<double>(data[i]);

    std::cout << "sino - min:" << min(sinogram) << " max: " << max(sinogram) << std::endl;
    int num_angles = sinogram.nrows();
    int num_pixels = sinogram.ncols();

    std::cout << "Num of projs: " << num_angles << std::endl;

    // angles
    auto anglesf = reader.read_angles("angs");
    std::vector<double> angles(anglesf.size());
    std::transform(std::execution::par_unseq, anglesf.begin(), anglesf.end(), angles.begin(),
        [](float f) { return static_cast<double>(f); });

    if (num_pixels % 2 == 0) {
        num_pixels -= 1;
    }

    Array<double> b(num_angles, num_pixels);
    for (int i = 0; i < num_angles; i++)
        for (int j = 0; j < num_pixels; j++)
            b(i, j) = random<double>();

    // transpose data
    auto sinoT = backward(b, angles);

    auto pg = nufft::polar_grid(angles, num_pixels);
    auto psf = calc_psf(pg);
    //writer.write("psf", psf);

    Array<double> recon(num_pixels, num_pixels);
    for (int i = 0; i < recon.size(); i++)
        recon[i] = random<double>();
    writer.write("input", recon);

    return 0;
}
