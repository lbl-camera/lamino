#include <fstream>
#include <iostream>

#include "tomocam.h"
#include "array.h"

const char * FILENAME = "/home/dkumar/data/shepp_logan/sino400.bin";
const int num_pixel = 400;
const int num_angles = 400;

int main(int argc, char **argv) {
    // read sinogram
    Array<float> sino(num_angles, num_pixel);
    sino.fromfile(FILENAME);
    // normalize
    sino /= sino.max();

    // angles
    Array < float >angles(1, num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / (num_angles - 1);

    float center = 250;
    auto image = backward(sino, angles, center);

	// write output
	image.tofile("output.bin");

    return 0;
}
