#include <fstream>
#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "tomocam.h"
#include "nufft.h"

const int num_angles = 360;
const int num_pixels = 511;

template <typename T>
T radrand() {
    T rmax = static_cast<T>(RAND_MAX);
    return static_cast<T>(rand()) / rmax;
}

int main() {

    // set seed
    std::srand(100);

    // angles
    std::vector<double> angles(num_angles);
    for (int i = 0; i < num_angles; i++)
        angles[i] = i * M_PI / num_angles;

    Array<double> recon(num_pixels, num_pixels);
    for (int i = 0; i < recon.size(); i++) recon[i] = 1.0;

    // sinogram
    Array<double> sino(num_angles, num_pixels);
    for (int i = 0; i < sino.size(); i++) sino[i] = radrand<double>();

    auto pg = nufft::polar_grid(angles, num_pixels);

    // calulate gradient normally
    auto t1 = forward(recon, angles);
    auto t2 = t1 - sino;
    auto g1 = backward(t2, angles);

    auto sinoT = backward(sino, angles);
    auto t3 = gradient(recon, pg);
    auto g2 = t3 - sinoT;
    std::cout << "g2.norm = " << norm(g2) << std::endl;
    g2.print();
    std::exit(1);
    // check if they match
    double a = dot(g1, g2) / dot(g1, g1);
    double e = norm(g2 - a * g1) / g1.size();
    std::cout << "a = " << a << ", e = " << e << std::endl;

    // error
    auto e1 = norm(t2);
    e1 = e1 * e1;

    // using gradient
    auto p1 = dot(recon, t3);
    auto p2 = 2 * dot(recon, sinoT);
    auto p3 = dot(sino, sino);
    auto e3 = (p1 - p2) / num_pixels + p3;
    std::cout << "e1 = " << e1 << ", e2 = " << e3
        << ", (|e2-e3|)/e3 = " << std::abs(e1 - e3) / e3 << std::endl;
    return 0;
}
