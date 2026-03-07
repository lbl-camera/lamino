/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals from the
 * U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 * IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 * the U.S. Government has been granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 * to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#include <cstdint>
#include <iostream>

#include "array.h"
#include "fft.h"
#include "tomocam.h"

using namespace tomocam;

template <typename T>
T rand_val() {
    return static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
}

int main() {

    // test 2-d real-to-complex with double
    Array<double> a(3, 5, 7);
    for (int i = 0; i < a.size(); i++) { a[i] = rand_val<double>(); }

    auto n = static_cast<double>(a.dims().n2 * a.dims().n3);
    auto b = fft::fft2_r2c(a);
    auto c = fft::fft2_c2r(b, a.dims());
    c /= n;

    bool err_2d_double = false;
    double max_err = 0.0;
    for (int i = 0; i < a.size(); i++) {
        double err = std::abs(a[i] - c[i]);
        max_err = std::max(max_err, err);
        if (err > 1e-10) { err_2d_double = true; }
    }
    std::cout << "test: 2-d array<double>: fft2_r2c,fft2_c2r: .... ";
    if (err_2d_double) {
        std::cout << "failed (max error: " << max_err << ")" << std::endl;
    } else {
        std::cout << "passed" << std::endl;
    }

    // test 2-d real-to-complex with float
    Array<float> d(5, 16, 16);
    for (int i = 0; i < d.size(); i++) { d[i] = rand_val<float>(); }

    auto m = static_cast<float>(d.dims().n2 * d.dims().n3);
    auto e = fft::fft2_r2c(d);
    auto f = fft::fft2_c2r(e, d.dims());
    f /= m;

    bool err_2d_float = false;
    for (int i = 0; i < d.size(); i++) {
        if (std::abs(d[i] - f[i]) > 1e-5) {
            err_2d_float = true;
            break;
        }
    }
    std::cout << "test: 2-d array<float>: fft2_r2c,fft2_c2r: .... ";
    if (err_2d_float) {
        std::cout << "failed" << std::endl;
    } else {
        std::cout << "passed" << std::endl;
    }

    // test with larger array (multiple batches) with double
    // Note: using n1=1 to work with current implementation
    Array<double> g(5, 16, 32);
    for (int i = 0; i < g.size(); i++) { g[i] = rand_val<double>(); }

    auto p = static_cast<double>(g.dims().n3 * g.dims().n2);
    auto h = fft::fft2_r2c(g);
    auto k = fft::fft2_c2r(h, g.dims());
    k /= p;

    bool err_3d_double = false;
    double max_err_3d = 0.0;
    for (int i = 0; i < g.size(); i++) {
        double err = std::abs(g[i] - k[i]);
        max_err_3d = std::max(max_err_3d, err);
        if (err > 1e-10) { err_3d_double = true; }
    }
    std::cout << "test: larger array<double>: fft2_r2c,fft2_c2r: .... ";
    if (err_3d_double) {
        std::cout << "failed (max error: " << max_err_3d << ")" << std::endl;
    } else {
        std::cout << "passed" << std::endl;
    }

    return 0;
}
