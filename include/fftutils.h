// clang-format off
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
 //clang-format on

#ifndef FFTUTILS__H
#define FFTUTILS__H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <execution>
#include <stdexcept>

#include "array.h"

namespace tomocam::fft {

    enum class Axes { one, two, three };

    template <typename T>
    Array<T> fftshift1(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);
        
        size_t nz = dims.z();
        size_t half_z = nz / 2;
        size_t rem_z = nz - half_z;
        
        for (size_t i = 0; i < dims.x(); ++i) {
            for (size_t j = 0; j < dims.y(); ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, half_z}],
                          &output[{i, j, rem_z}]);
                std::copy(&input[{i, j, half_z}], &input[{i, j, nz}],
                          &output[{i, j, 0}]);
            }
        }
        return output;
    }

    template <typename T>
    Array<T> fftshift2(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);
        
        size_t ny = dims.y();
        size_t nz = dims.z();
        size_t half_y = ny / 2;
        size_t half_z = nz / 2;
        size_t rem_y = ny - half_y;
        size_t rem_z = nz - half_z;
        
        for (size_t i = 0; i < dims.x(); ++i) {
            for (size_t j = 0; j < half_y; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, half_z}],
                          &output[{i, j + rem_y, rem_z}]);
                std::copy(&input[{i, j, half_z}], &input[{i, j, nz}],
                          &output[{i, j + rem_y, 0}]);
            }
            for (size_t j = half_y; j < ny; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, half_z}],
                          &output[{i, j - half_y, rem_z}]);
                std::copy(&input[{i, j, half_z}], &input[{i, j, nz}],
                          &output[{i, j - half_y, 0}]);
            }
        }
        return output;
    }

    template <typename T>
    Array<T> fftshift3(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);
        
        size_t nx = dims.x();
        size_t ny = dims.y();
        size_t nz = dims.z();
        size_t half_x = nx / 2;
        size_t half_y = ny / 2;
        size_t half_z = nz / 2;
        size_t rem_x = nx - half_x;
        size_t rem_y = ny - half_y;
        size_t rem_z = nz - half_z;
        
        for (size_t i = 0; i < half_x; ++i) {
            for (size_t j = 0; j < half_y; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, half_z}],
                          &output[{i + rem_x, j + rem_y, rem_z}]);
                std::copy(&input[{i, j, half_z}], &input[{i, j, nz}],
                          &output[{i + rem_x, j + rem_y, 0}]);
            }
            for (size_t j = half_y; j < ny; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, half_z}],
                          &output[{i + rem_x, j - half_y, rem_z}]);
                std::copy(&input[{i, j, half_z}], &input[{i, j, nz}],
                          &output[{i + rem_x, j - half_y, 0}]);
            }
        }
        for (size_t i = half_x; i < nx; ++i) {
            for (size_t j = 0; j < half_y; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, half_z}],
                          &output[{i - half_x, j + rem_y, rem_z}]);
                std::copy(&input[{i, j, half_z}], &input[{i, j, nz}],
                          &output[{i - half_x, j + rem_y, 0}]);
            }
            for (size_t j = half_y; j < ny; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, half_z}],
                          &output[{i - half_x, j - half_y, rem_z}]);
                std::copy(&input[{i, j, half_z}], &input[{i, j, nz}],
                          &output[{i - half_x, j - half_y, 0}]);
            }
        }
        return output;
    }

    template <typename T>
    Array<T> ifftshift1(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);
        
        size_t nz = dims.z();
        size_t half_z = nz / 2;
        size_t rem_z = nz - half_z;
        
        for (size_t i = 0; i < dims.x(); ++i) {
            for (size_t j = 0; j < dims.y(); ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, rem_z}],
                          &output[{i, j, half_z}]);
                std::copy(&input[{i, j, rem_z}], &input[{i, j, nz}],
                          &output[{i, j, 0}]);
            }
        }
        return output;
    }

    template <typename T>
    Array<T> ifftshift2(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);
        
        size_t ny = dims.y();
        size_t nz = dims.z();
        size_t half_y = ny / 2;
        size_t half_z = nz / 2;
        size_t rem_y = ny - half_y;
        size_t rem_z = nz - half_z;
        
        for (size_t i = 0; i < dims.x(); ++i) {
            for (size_t j = 0; j < rem_y; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, rem_z}],
                          &output[{i, j + half_y, half_z}]);
                std::copy(&input[{i, j, rem_z}], &input[{i, j, nz}],
                          &output[{i, j + half_y, 0}]);
            }
            for (size_t j = rem_y; j < ny; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, rem_z}],
                          &output[{i, j - rem_y, half_z}]);
                std::copy(&input[{i, j, rem_z}], &input[{i, j, nz}],
                          &output[{i, j - rem_y, 0}]);
            }
        }
        return output;
    }

    template <typename T>
    Array<T> ifftshift3(const Array<T> &input) {
        auto dims = input.dims();
        Array<T> output(dims);
        
        size_t nx = dims.x();
        size_t ny = dims.y();
        size_t nz = dims.z();
        size_t half_x = nx / 2;
        size_t half_y = ny / 2;
        size_t half_z = nz / 2;
        size_t rem_x = nx - half_x;
        size_t rem_y = ny - half_y;
        size_t rem_z = nz - half_z;
        
        for (size_t i = 0; i < rem_x; ++i) {
            for (size_t j = 0; j < rem_y; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, rem_z}],
                          &output[{i + half_x, j + half_y, half_z}]);
                std::copy(&input[{i, j, rem_z}], &input[{i, j, nz}],
                          &output[{i + half_x, j + half_y, 0}]);
            }
            for (size_t j = rem_y; j < ny; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, rem_z}],
                          &output[{i + half_x, j - rem_y, half_z}]);
                std::copy(&input[{i, j, rem_z}], &input[{i, j, nz}],
                          &output[{i + half_x, j - rem_y, 0}]);
            }
        }
        for (size_t i = rem_x; i < nx; ++i) {
            for (size_t j = 0; j < rem_y; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, rem_z}],
                          &output[{i - rem_x, j + half_y, half_z}]);
                std::copy(&input[{i, j, rem_z}], &input[{i, j, nz}],
                          &output[{i - rem_x, j + half_y, 0}]);
            }
            for (size_t j = rem_y; j < ny; ++j) {
                std::copy(&input[{i, j, 0}], &input[{i, j, rem_z}],
                          &output[{i - rem_x, j - rem_y, half_z}]);
                std::copy(&input[{i, j, rem_z}], &input[{i, j, nz}],
                          &output[{i - rem_x, j - rem_y, 0}]);
            }
        }
        return output;
    }

    template <typename T>
    Array<T> fftshift(const Array<T> &input, Axes axes) {
        if (axes == Axes::one) {
            return fftshift1(input);
        } else if (axes == Axes::two) {
            return fftshift2(input);
        } else if (axes == Axes::three) {
            return fftshift3(input);
        } else {
            throw std::runtime_error("unknown shift axes");
        }
    }

    template <typename T>
    Array<T> ifftshift(const Array<T> &input, Axes axes) {
        if (axes == Axes::one) {
            return ifftshift1(input);
        } else if (axes == Axes::two) {
            return ifftshift2(input);
        } else if (axes == Axes::three) {
            return ifftshift3(input);
        } else {
            throw std::runtime_error("unknown shift axes");
        }
    }

} // namespace tomocam::fft

#endif // FFTUTILS__H
