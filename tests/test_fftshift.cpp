#include <iostream>

#include "array.h"
#include "fftutils.h"

using namespace tomocam;

void print_array(const Array<float> &a) {
    for (size_t i = 0; i < a.dims().n2; ++i) {
        for (size_t j = 0; j < a.dims().n3; ++j) {
            std::cout << a[{0, i, j}] << " ";
        }
        std::cout << "\n";
    }
}

int main() {

    dims_t dims = {1, 1, 13};
    Array<float> a(dims);
    for (size_t j = 0; j < a.size(); ++j) { a[j] = static_cast<float>(j); }
    std::cout << "O-freq at the beginning: \n";
    print_array(a);

    std::cout << "\n";
    std::cout << "fftshift2: \n";
    auto b = fft::fftshift2(a);
    print_array(b);

    std::cout << "\n";
    std::cout << "ifftshift2: \n";
    auto c = fft::ifftshift2(a);
    print_array(c);

    // create a centered array
    Array<float> d(dims);
    for (size_t j = 0; j < a.size(); ++j) {
        d[j] = static_cast<float>(j) - dims.n3 / 2;
    }

    std::cout << "\n\n";
    std::cout << "O-freq at the center: \n";
    print_array(d);

    std::cout << "\n";
    std::cout << "fftshift2: \n";
    auto e = fft::fftshift2(d);
    print_array(e);

    std::cout << "\n";
    std::cout << "ifftshift2: \n";
    auto f = fft::ifftshift2(d);
    print_array(f);

    return 0;
}
