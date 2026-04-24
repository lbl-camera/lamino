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
#ifndef POLAR_GRID__H
#define POLAR_GRID__H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "array.h"
namespace tomocam {

    template <typename T>
    struct PolarGrid {
        size_t npts;
        std::vector<T> theta;
        Array<T> x;
        Array<T> y;
        Array<T> z;

        // default constructor
        PolarGrid() : npts(0) {}

        // constructor
        PolarGrid(const std::vector<T> &angles, size_t nrows, size_t ncols,
                  T gamma) {

            theta = angles;
            dims_t dims = dims_t{theta.size(), nrows, ncols};
            x = Array<T>(dims);
            y = Array<T>(dims);
            z = Array<T>(dims);
            npts = dims.size();

            // compute grid points
            T dX = (2 * M_PI) / static_cast<T>(ncols);
            T dY = (2 * M_PI) / static_cast<T>(nrows);

#pragma omp parallel for collapse(3)
            for (size_t i = 0; i < dims.n1; ++i) {
                for (size_t j = 0; j < dims.n2; ++j) {
                    for (size_t k = 0; k < dims.n3; ++k) {

                        T qX = k * dX - M_PI + dX / 2;
                        T qY = j * dY - M_PI + dY / 2;

                        // apply rotations
                        z[{i, j, k}] = -qX * std::sin(gamma) +
                                       qY * std::cos(theta[i]) * std::cos(gamma);
                        y[{i, j, k}] = -qY * std::sin(theta[i]);
                        x[{i, j, k}] = qX * std::cos(gamma) +
                                       qY * std::cos(theta[i]) * std::sin(gamma);
                    }
                }
            }
        }
        // delete copy constructor and assignment
        PolarGrid(const PolarGrid<T> &) = delete;
        PolarGrid<T> &operator=(const PolarGrid<T> &) = delete;

        // move constructor and assignment
        PolarGrid(PolarGrid<T> &&other) noexcept
            : npts(other.npts), theta(std::move(other.theta)), x(std::move(other.x)),
              y(std::move(other.y)), z(std::move(other.z)) {}

        PolarGrid<T> &operator=(PolarGrid<T> &&other) noexcept {
            if (this != &other) {
                npts = other.npts;
                theta = std::move(other.theta);
                x = std::move(other.x);
                y = std::move(other.y);
                z = std::move(other.z);
            }
            return *this;
        }

        // array dimensions for non-uniform points
        [[nodiscard]] dims_t dims() const { return x.dims(); }

        // size of the array
        [[nodiscard]] size_t size() const { return x.size(); }

        // theta values
        [[nodiscard]] const std::vector<T> &angles() const { return theta; }

        // number of angles
        [[nodiscard]] size_t nprojs() const { return theta.size(); }
    };

} // namespace tomocam

#endif // POLAR_GRID__H
