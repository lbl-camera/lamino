#include <iostream>
#include <fstream>

#include "array.h"
#include "kernel.h"
#include "convolve.h"

const int dims[] = {400, 400};
const int num_angles = 360;
int main(int argc, char **argv) {

    // read data
    size_t size = dims[0] * dims[1];
    float * data = new float[size];

    // read sinogram
    std::fstream fp("shepp_logan.bin", std::ifstream::in);
    if (! fp.is_open()) {
        std::cerr << "error! unable to open data file." << std::endl;
        exit(1);
    }
    fp.read((char *) data, sizeof(float) * size);
    fp.close();

    // angles
    Array<float> angles(1, 1, num_angles);
    for (int i = 0; i < num_angles; i++) 
        angles[i] = i * (M_PI / num_angles-1);

    // oversampling
    float oversample = 1.5;
    if (argc >= 2) {
        oversample = atof(argv[1]);
    }

    // convolution kernel radius
    float radius = 2.f;
    if (argc >= 3)
        radius = atof(argv[2]);
    // kaiser window scaling  (default 4 * PI)
    float beta = 12.566371; 
    if (argc >= 4)
        beta = atof(argv[3]);

    Kernel kernel(radius, beta);

    // input and output with padding
    int ipad = (int) ((oversample - 1.f) * dims[1] / 2);
    int padded = dims[2] + 2 * ipad;
  
    // working arrays  
    Array<complex_t> image(1, padded, padded);
    Array<complex_t> sino(1, num_angles, padded);

    dim3_t idims = image.dims();
    #pragma omp parallel for 
    for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
            image(0, i+ipad, j+ipad) = data[i * dims[1] + j];

    forward(image, sino, angles, kernel);

    // reset zero-padding
    idims = sino.dims();
    for (int i = 0; i < idims.x; i++) 
        for (int j = 0; j < idims.y; j++) 
            for (int k = 0; k < idims.z; k++)
                if ((k < ipad) || (k >= idims.z + ipad))
                    sino(i,j,k) = 0;
    image.clear();

    backward(sino, image, angles, kernel);
    auto e = calc_error(data, image, ipad);
    write_output(image, ipad);
    std::cout << "Error = " << e << std::endl;

    return 0;
}
