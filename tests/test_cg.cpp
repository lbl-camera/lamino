#include <array>
#include <cstddef>
#include <format>
#include <iostream>
#include <vector>

#include "array.h"
#include "tiff.h"
#include "tomocam.h"

typedef tomocam::Array<float> Array;

struct Kernel {
    std::array<size_t, 2> dims;
    std::vector<float> weights;
    Kernel() {
        dims = {5, 5};
        float w[5] = {0.1f, 0.2f, 0.5f, 0.2f, 0.1f};
        float sum = 0.f;
        weights.resize(25, 0.0f);
        for (size_t i = 0; i < dims[0]; i++) {
            for (size_t j = 0; j < dims[1]; j++) {
                weights[i * dims[1] + j] = w[i] * w[j];
                sum += w[i] * w[j];
            }
        }
        for (auto &w : weights) { w /= sum; }
    }

    float &operator()(size_t i, size_t j) { return weights[i * dims[1] + j]; }
    const float &operator()(size_t i, size_t j) const {
        return weights[i * dims[1] + j];
    }
};

void create_ground_truth(Array &x_true) {
    /**
     * Create an image with paritally overlapping rectangles as ground
     * truth. rect: (x0, y0), dx, dy, value rect0: 50, 50), 100, 150, 0.5
     *  rect1: (120, 80), 200, 100,
     *  rect2: (255, 255), 100, 200,
     *
     * Args:
     *   x_true: 3D array to be filled with the ground truth pattern
     */
    size_t nslices = x_true.nslices();
    size_t nrows = x_true.nrows();
    size_t ncols = x_true.ncols();

    for (size_t i = 0; i < nslices; i++) {
        for (size_t j = 0; j < nrows; j++) {
            for (size_t k = 0; k < ncols; k++) {
                x_true[{i, j, k}] = 0.0f;
                // rect0
                if (j >= 50 && j < 150 && k >= 50 && k < 200) {
                    x_true[{i, j, k}] += 0.5f;
                }
                // rect1
                if (j >= 80 && j < 180 && k >= 120 && k < 220) {
                    x_true[{i, j, k}] += 1.25f;
                }
                // rect2
                if (j >= 255 && j < 355 && k >= 255 && k < 455) {
                    x_true[{i, j, k}] += 1.f;
                }
            }
        }
    }
}

Array forward(const Array &input, const Kernel &kernel) {

    /**
     * Perform a forward pass of a convolutional layer.
     *
     * Args:
     *   input: 3D input array (nslices, nrows, ncols)
     *   kernel: 2D convolution kernel
     *
     * Returns:
     *   3D output array after applying the convolution
     */
    Array output(input.dims());
    size_t nslices = input.nslices();
    size_t nrows = input.nrows();
    size_t ncols = input.ncols();
    size_t wx = kernel.dims[0];
    size_t wy = kernel.dims[1];

    for (size_t i = 0; i < nslices; i++) {
        for (size_t j = 0; j < nrows; j++) {
            for (size_t k = 0; k < ncols; k++) {
                float sum = 0.0f;
                for (size_t ki = 0; ki < wx; ki++) {
                    for (size_t kj = 0; kj < wy; kj++) {
                        size_t jj = j + ki - wx / 2;
                        size_t kk = k + kj - wy / 2;
                        if (jj < nrows && kk < ncols) {
                            sum += input[{i, jj, kk}] * kernel(ki, kj);
                        }
                    }
                }
                output[{i, j, k}] = sum;
            }
        }
    }

    return output;
}

Array adjoint(const Array &input, const Kernel &kernel) {

    /**
     * Perform the adjoint (backward) pass of a convolutional layer.
     *
     * Args:
     *   input: 3D input array (nslices, nrows, ncols)
     *   kernel: 2D convolution kernel
     *
     * Returns:
     *   3D output array after applying the adjoint convolution
     */
    Array output(input.dims());

    size_t nslices = input.nslices();
    size_t nrows = input.nrows();
    size_t ncols = input.ncols();
    size_t wx = kernel.dims[0];
    size_t wy = kernel.dims[1];

    for (size_t i = 0; i < nslices; i++) {
        for (size_t j = 0; j < nrows; j++) {
            for (size_t k = 0; k < ncols; k++) {
                float sum = 0.0f;
                for (size_t ki = 0; ki < wx; ki++) {
                    for (size_t kj = 0; kj < wy; kj++) {
                        size_t jj = j - ki + wx / 2;
                        size_t kk = k - kj + wy / 2;
                        if (jj < nrows && kk < ncols) {
                            sum += input[{i, jj, kk}] * kernel(ki, kj);
                        }
                    }
                }
                output[{i, j, k}] = sum;
            }
        }
    }

    return output;
}

int main(int argc, char **argv) {

    size_t N = 511;
    if (argc == 2) { N = (size_t)atoi(argv[1]); }

    Array x_true(1, N, N);
    Kernel kernel;

    create_ground_truth(x_true);
    // blur + noise
    auto tmp = Array(x_true.dims());
    tmp.fill_random();
    tmp += x_true;
    auto x_blur = forward(tmp, kernel);

    auto yT = adjoint(x_blur, kernel);
    Array x0 = Array(x_true.dims());
    x0.fill(0.f);

    auto backfwd = [&kernel](const Array &x) {
        auto tmp = forward(x, kernel);
        return adjoint(tmp, kernel);
    };
    float tol = 1.e-05;

    auto x_opt = tomocam::opt::cgsolver<float>(backfwd, yT, 50, tol);

    tomocam::tiff::write("x_opt.tif", x_opt);
    return 0;
}
