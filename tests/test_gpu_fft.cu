// cuFFT plans must match the precision of the data: a double-precision array
// has to be transformed completely (an FFT plan of the wrong precision only
// touches half of the bytes and leaves garbage/NaN in the rest).
// Checks, for float and double: batched fft2d -> ifft2d round trip, and the
// zero-frequency term of fft2d against the sum of each slice.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "array.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/fft.h"

using namespace tomocam;

template <typename T>
static bool run(const char *name, double tol) {
    const size_t nb = 4, n = 16;
    Array<T> h(nb, n, n);
    for (size_t i = 0; i < h.size(); ++i) {
        h[i] = static_cast<T>(std::rand()) / static_cast<T>(RAND_MAX);
    }

    gpu::DeviceArray<T> d(h);
    auto c = gpu::array::to_complex(d);
    auto f = gpu::fft::fft2d(c);
    auto r = gpu::fft::ifft2d(f);
    auto back = (gpu::array::to_real(r) / static_cast<T>(n * n)).to_host();
    auto fh = gpu::array::to_real(f).to_host(); // real part of the spectrum

    double err = 0, ref = 0;
    for (size_t i = 0; i < h.size(); ++i) {
        err += std::pow(static_cast<double>(back[i]) - h[i], 2);
        ref += std::pow(static_cast<double>(h[i]), 2);
    }
    err = std::sqrt(err / ref);

    // DC term of every slice = sum of the slice
    double dc_err = 0;
    for (size_t s = 0; s < nb; ++s) {
        double sum = 0;
        for (size_t i = 0; i < n * n; ++i) { sum += h[s * n * n + i]; }
        dc_err = std::max(dc_err,
                          std::abs(static_cast<double>(fh[s * n * n]) - sum) / sum);
    }

    bool ok = std::isfinite(err) && err < tol && dc_err < tol;
    std::cout << "test_gpu_fft<" << name << ">: round-trip err " << err
              << ", DC err " << dc_err << (ok ? "  passed\n" : "  FAILED\n");
    return ok;
}

int main() {
    bool ok = run<float>("float", 1e-5);
    ok &= run<double>("double", 1e-12);
    return ok ? 0 : 1;
}
