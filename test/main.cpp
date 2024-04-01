#include <fstream>
#include <iostream>
#include <ctime>

#include "array.h"
#include "tomocam.h"
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
    tomocam::H5Writer writer("twocode.h5");

    reader.setDataset("projs"); 
    auto data = reader.read_sinogram(0);
    Array<double> sinogram(data.nrows(), data.ncols());
    for (int i = 0; i < data.size(); i++)
        sinogram[i] = static_cast<double>(data[i]);

    std::cout << "sino - min:" << sinogram.min() << " max: " << sinogram.max() << std::endl;
    int num_angles = sinogram.nrows();
    int num_pixels = sinogram.ncols();

    std::cout << "Num of projs: " << num_angles << std::endl;

    // angles
    auto anglesf = reader.read_angles("angs");
    Array<double> angles(anglesf.nrows(), anglesf.ncols());
    for (int i = 0; i < anglesf.size(); i++)
        angles[i] = static_cast<double>(anglesf[i]);


    if (num_pixels % 2 == 0) {
        num_pixels -= 1;
    }
    Array<double> b(num_angles, num_pixels);
    for (int i = 0; i < num_angles; i++)
        for (int j = 0; j < num_pixels; j++)
            b(i, j) = random<double>();

    // transpose data
    auto sinoT = backward(b, angles, center);
    //writer.write("sinoT", sinoT);

    // point spread function
    auto psf = calc_psf(angles, num_pixels);
    //writer.write("psf", psf);
   
    
    Array<double> recon(num_pixels, num_pixels);
    for (int i = 0; i < recon.size(); i++)
        recon[i] = random<double>();
    writer.write("input", recon);
    auto x2 = recon;

    // gradient 1
    auto t1 = forward(recon, angles, center);
    auto err1 = (t1 - b).dot(t1 - b);
    auto g1 = backward(t1,  angles, center) - sinoT;
    writer.write("nufft", g1);

    // gradient 2
    auto t2 = fftconvolve(x2, psf);
    auto g2 = t2-sinoT;
    writer.write("toeplitz", g2);
    
    auto err = (g1 - g2).norm();
    std::cout << "Error: " << err/g1.norm() << std::endl;

    // function evaluation
    auto p1 = x2.dot(t2);
    auto p2 = 2 * x2.dot(sinoT);
    auto p3 = b.norm();
    std::cout << "F1 = " << err1 << std::endl;
    std::cout << "F2 = " << p1 - p2 + p3 << std::endl;
 
    return 0;
}
