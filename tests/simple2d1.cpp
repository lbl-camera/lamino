// this is all you must include for the finufft lib...
#include <finufft.h>
#include <complex>

// also needed for this example...
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "tomocam.h"
#include "fft.h"
#include "array.h"

int main(int argc, char *argv[]){

/* Simple 2D type-1 example of calling the FINUFFT library from C++, using plain
   arrays of C++ complex numbers, with a math test. Double precision version. 

   Compile multithreaded with
   g++ -fopenmp simple2d1.cpp -I ../src ../lib-static/libfinufft.a -o simple2d1 -lfftw3 -lfftw3_omp -lm
   single core with:
   g++ simple2d1.cpp -I ../src ../lib-static/libfinufft.a -o simple2d1 -lfftw3 -lm
   
   Usage:  ./simple2d1
*/

  int M1 = 400;
  int M2 = 400;   
  int N1 = 400;
  int N2 = 400;  
  double tol = 1e-6;           // desired accuracy
  nufft_opts opts; finufftf_default_opts(&opts);

    

  // generate non-uniform points on (x,y) and complex strengths (c):
  float * angle = new float [M1];
  for (int i = 0; i < M1; i++)
    angle[i] = i * M_PI / (M1-1);
    
  int M = M1 * M2;
  float * x = new float [M];
  float * y = new float [M];

  float scale = 2 * M_PI / M2;
  float center = M2 / 2.f;

  for(int i = 0; i < M1; i++)
    for(int j = 0; j < M2; j++){
        float r = scale * (j - center);
        x[i * M2 + j] =  scale * std::cos(angle[i]);
        y[i * M2 + j] =  scale * std::sin(angle[i]);
    }


  // each component uniform random in [-1,1]
  Array<complex_t> C(M1, M2);
  for (int i = 0; i < M1; i++) 
    for (int j = 0; j < M2; j++) 
        C(i, j) = std::cos(x[i*M2+j]) * std::sin(y[i*M2+j]);

  C.real().tofile("output.bin");
  std::exit(1);
  // output array for the Fourier modes
  Array<complex_t> F(N1, N2);

  // call the NUFFT (with iflag += 1): note passing in pointers...
  opts.upsampfac = 2;
  int ier = finufftf2d1(M, x, y, C.ptr(), -1, tol, N1, N2, F.ptr(), &opts);

  fftshift2(F);
  ifft2(F);
  fftshift2(F);
  
  F.real().tofile("output.bin");
  return ier;
}
