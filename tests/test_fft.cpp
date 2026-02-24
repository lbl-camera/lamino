#include <cstdint>
#include <iostream>

#include "array.h"
#include "array_ops.h"
#include "fft.h"
#include "padding.h"
#include "tomocam.h"

using namespace tomocam;

template <typename T>
T rand() {
    return static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
}

int main() {

    // test 1-d array double
    Array<double> a(1, 1, 8);
    for (int i = 0; i < a.size(); i++) { a[i] = rand<double>(); }

    double factor = std::sqrt(2);

    // test 2-d array float
    Array<float> d(1, 8, 8);
    for (int i = 0; i < d.size(); i++) { d[i] = rand<float>(); }

    // pad -> fft -> ifft -> unpad
    auto dp = pad2d<float>(d, factor, PadType::RIGHT);
    auto m = static_cast<std::complex<float>>(dp.size());
    auto dc = array::to_complex(dp);
    auto e = fft::fft2(dc);
    auto ep = fft::ifft2(e) / m;
    auto f = array::to_real(ep);

    auto crop_size = dp.dims() - d.dims();
    f = crop2d(f, crop_size, PadType::RIGHT);

    // check error
    bool err_2d_float = false;
    for (int i = 0; i < d.size(); i++) {
        if (std::abs(d[i] - f[i]) > 1e-6) {
            err_2d_float = true;
            break;
        }
    }
    std::cout << "test: 2-d array<float>: pad,fft,ifft,crop: .... ";
    if (err_2d_float) {
        std::cout << "failed" << std::endl;
    } else {
        std::cout << "passed" << std::endl;
    }

    return 0;
}
