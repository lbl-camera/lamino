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

#ifndef FFTUTILS__H
#define FFTUTILS__H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <execution>
#include <stdexcept>
#include <vector>

#include "array.h"

namespace tomocam::fft {

    template <typename T>
    std::vector<T> fftfreq(size_t n) {
        std::vector<T> freqs(n);
        size_t n2 = n % 2 == 0 ? n / 2 : n / 2 + 1;
        for (size_t i = 0; i < n; ++i) {
            if (i < n2) {
                freqs[i] = static_cast<T>(i) / static_cast<T>(n);
            } else {
                freqs[i] = static_cast<T>(i - n) / static_cast<T>(n);
            }
        }
        return freqs;
    }

    template <typename T>
    Array<T> fftshift1(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);

        size_t d = dims.n3 / 2;

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < dims.n1; ++i) {
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t k = 0; k < dims.n3; ++k) {
                    size_t k2 = (k + d) % dims.n3;
                    output[{i, j, k}] = input[{i, j, k2}];
                }
            }
        }
        return output;
    }

    template <typename T>
    Array<T> fftshift2(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);

        size_t d1 = dims.n2 / 2;
        size_t d2 = dims.n3 / 2;

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < dims.n1; ++i) {
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t k = 0; k < dims.n3; ++k) {
                    size_t j2 = (j + d1) % dims.n2;
                    size_t k2 = (k + d2) % dims.n3;
                    output[{i, j, k}] = input[{i, j2, k2}];
                }
            }
        }
        return output;
    }

    template <typename T>
    Array<T> fftshift3(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);

        size_t d1 = dims.n1 / 2;
        size_t d2 = dims.n2 / 2;
        size_t d3 = dims.n3 / 2;

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < dims.n1; ++i) {
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t k = 0; k < dims.n3; ++k) {
                    size_t i2 = (i + d1) % dims.n1;
                    size_t j2 = (j + d2) % dims.n2;
                    size_t k2 = (k + d3) % dims.n3;
                    output[{i, j, k}] = input[{i2, j2, k2}];
                }
            }
        }
        return output;
    }

    template <typename T>
    Array<T> ifftshift1(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);

        size_t d = dims.n3 / 2;
        if (dims.n3 % 2 == 1) { d += 1; }
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < dims.n1; ++i) {
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t k = 0; k < dims.n3; ++k) {
                    size_t k2 = (k + d) % dims.n3;
                    output[{i, j, k}] = input[{i, j, k2}];
                }
            }
        }
        return output;
    }

    template <typename T>
    Array<T> ifftshift2(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);

        size_t d1 = dims.n2 / 2;
        if (dims.n2 % 2 == 1) { d1 += 1; }
        size_t d2 = dims.n3 / 2;
        if (dims.n3 % 2 == 1) { d2 += 1; }
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < dims.n1; ++i) {
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t k = 0; k < dims.n3; ++k) {
                    size_t j2 = (j + d1) % dims.n2;
                    size_t k2 = (k + d2) % dims.n3;
                    output[{i, j, k}] = input[{i, j2, k2}];
                }
            }
        }

        return output;
    }

    template <typename T>
    Array<T> ifftshift3(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);

        size_t d1 = dims.n1 / 2;
        if (dims.n1 % 2 == 1) { d1 += 1; }
        size_t d2 = dims.n2 / 2;
        if (dims.n2 % 2 == 1) { d2 += 1; }
        size_t d3 = dims.n3 / 2;
        if (dims.n3 % 2 == 1) { d3 += 1; }
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < dims.n1; ++i) {
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t k = 0; k < dims.n3; ++k) {
                    size_t i2 = (i + d1) % dims.n1;
                    size_t j2 = (j + d2) % dims.n2;
                    size_t k2 = (k + d3) % dims.n3;
                    output[{i, j, k}] = input[{i2, j2, k2}];
                }
            }
        }
        return output;
    }
} // namespace tomocam::fft

#endif // FFTUTILS__H
