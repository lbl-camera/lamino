#include <array>
#include <cstddef>
#include <format>
#include <iostream>
#include <vector>

#include "array.h"
#include "optimize.h"
#include "tiff.h"
#include "tomocam.h"

typedef tomocam::Array<float> Array;

int main(int argc, char **argv) {

    size_t N = 100;
    if (argc == 2) { N = (size_t)atoi(argv[1]); }
    tomocam::dims_t dims = {1, 1, N};

    // create vector representing a diagonal matrix of size N
    std::vector<float> diagMat(N);
    for (size_t i = 0; i < N; i++) { diagMat[i] = static_cast<float>(i + 1); }

    auto x_true = tomocam::Array<float>::random(dims);

    // create b = A*x_true
    tomocam::Array<float> b(dims);
    for (size_t i = 0; i < N; i++) { b[{0, 0, i}] = diagMat[i] * x_true[{0, 0, i}]; }

    // Conjugate Gradient solver for diagonal matrix (vector version)
    using VecArray = std::array<tomocam::Array<float>, 3>;
    std::function<VecArray(const VecArray &)> A = [&](const VecArray &v) {
        VecArray Av;
        for (size_t j = 0; j < 3; j++) {
            Av[j] = tomocam::Array<float>(dims);
            for (size_t i = 0; i < N; i++) {
                Av[j][{0, 0, i}] = diagMat[i] * v[j][{0, 0, i}];
            }
        }
        return Av;
    };

    // call CG solver
    // initialize guess
    VecArray x;
    VecArray b_vec;
    for (size_t j = 0; j < 3; j++) {
        x[j] = tomocam::Array<float>::random(dims);
        b_vec[j] = tomocam::Array<float>(dims);
        for (size_t i = 0; i < N; i++) {
            b_vec[j][{0, 0, i}] = diagMat[i] * x_true[{0, 0, i}];
        }
    }
    auto x_sol = tomocam::opt::cgsolver(A, b_vec, x, 1000, 1e-6f, 0.0f);

    return 0;
}
