#include <fstream>
#include <iostream>
#include <cstdlib>

#include <finufft.h>
#include "array.h"
#include "tomocam.h"
#include "dtypes.h"

const int MAX_ITERS = 1;
const int num_pixels = 400;
const int num_angles = 400;
const char * FILENAME = "/home/dkumar/data/shepp_logan/sino400.bin";

int main(int argc, char **argv) {

    // read data
    srand(100);
    Array<float> sinogram(num_angles, num_pixels);
    sinogram.fromfile(FILENAME);

    // normalize
    sinogram /= sinogram.max();

    //for (int i = 0; i < sinogram.size(); i++)
    //    sinogram[i] = (float) rand() / RAND_MAX;
   
    // angles
    Array < float >angles(1, num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / (num_angles - 1);

    // x-y coordinates for sinogram values
    int M1 = num_angles;
    int M2 = num_pixels;
    int M = num_angles * num_pixels;
    float scale = 2 * M_PI / M2;
    float center = M2 / 2.f;
    float * x = new float[M];
    float * y = new float[M];
    for(int i = 0; i < M1; i++)
        for(int j = 0; j < M2; j++){
            float r = scale * (j - center);
            x[i * M2 + j] =  r * std::cos(angles[i]);
            y[i * M2 + j] =  r * std::sin(angles[i]);
        }

    // working arrays
    Array <complex_t> C(num_pixels, num_pixels);
    Array <complex_t> F(num_angles, num_pixels);
    for (int i =0; i < M; i++)
        C[i] = complex_t(sinogram[i],0);

    double tol = 1e-06; // desired accuracy
    finufft_opts opts; finufftf_default_opts(&opts);
    opts.upsampfac = 2;
    
    int ier = finufftf2d1(M, y, x, C.ptr(), -1, tol, num_pixels, num_pixels, F.ptr(), &opts);

    // frequencies
    float * kx = new float[num_pixels];
    float * ky = new float[num_pixels];
    for (int j = 0; j < num_pixels; j++) {
        kx[j] = j - num_pixels / 2;
        ky[j] = kx[j];
    }
            

    // calculate exactly
    Array<complex_t> Ft(num_pixels, num_pixels);
    #pragma omp parallel for
    for (int i = 0; i < num_pixels; i++)
        for (int j = 0; j < num_pixels; j++) {
            Ft(i,j) = 0;
            for (int m = 0; m < M; m++) {
                auto t = I * (kx[j] * x[m] + ky[i] * y[m]);
                Ft(i,j) += sinogram[m] * std::exp(-t);
            }
        }

    
    float err = 0;
    #pragma omp parallel for
    for (int i = 0; i < F.size(); i++)
        err += std::pow(std::abs(Ft[i] - F[i]),2);
  
    std::cout << "error = " << err << std::endl;

    delete [] x;
    delete [] y;
    delete [] kx;
    delete [] ky;
    return 0;
}
