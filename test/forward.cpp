#include <fstream>
#include <iostream>

#include "tomocam.h"
#include "array.h"
#include "adam.h"

const int MAX_ITERS = 1;
const int num_angles = 400;
const int num_pixels = 400;
const char * FILENAME = "/home/dkumar/data/shepp_logan/shepp400.bin";

int main(int argc, char **argv)
{

    // read data
    Array<float> image(num_pixels, num_pixels);
    image.fromfile(FILENAME);

    Array < float >angles(1, num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / (num_angles - 1);

    // radon transform
    float center = 200;
	auto sino = forward(image, angles, center);
		
    sino.tofile("output.bin");
    return 0;
}
