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
#ifndef POLAR_GRID_H
#define POLAR_GRID_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <iostream>
#include <ranges>
#include <vector>

#include "array.h"
#include "array_ops.h"

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
            npts = dims.n1 * dims.n2 * dims.n3;
            x = Array<T>(dims);
            y = Array<T>(dims);
            z = Array<T>(dims);
            // rotation matrix
            T cos_gamma = std::cos(gamma);
            T sin_gamma = std::sin(gamma);

            // compute grid points
            T L = 2 * M_PI;
            T dX = L / static_cast<T>(ncols);
            T dY = L / static_cast<T>(nrows);
            T L_half = L / 2;

#pragma omp parallel for collapse(3)
            for (size_t i = 0; i < dims.n1; ++i) {
                for (size_t j = 0; j < dims.n2; ++j) {
                    for (size_t k = 0; k < dims.n3; ++k) {

                        T qX = (k + 0.5) * dX - L_half;
                        T qY = (j + 0.5) * dY - L_half;

                        // apply rotations
                        x[{i, j, k}] = qX * std::cos(gamma) -
                                       qY * std::sin(gamma) * std::cos(theta[i]);
                        y[{i, j, k}] = qX * std::sin(gamma) +
                                       qY * std::cos(gamma) * std::cos(theta[i]);
                        z[{i, j, k}] = qY * std::sin(theta[i]);
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

        PolarGrid<T> clone() const {
            PolarGrid<T> out;
            out.npts = this->npts;
            out.theta = this->theta;
            out.x = this->x.clone();
            out.y = this->y.clone();
            out.z = this->z.clone();
            return out;
        }
        void print_limits() const {
            std::cout << "PolarGrid limits:\n";
            std::cout << std::format("qx: [{}, {}]\n", array::min(x), array::max(x));
            std::cout << std::format("qy: [{}, {}]\n", array::min(y), array::max(y));
            std::cout << std::format("qz: [{}, {}]\n", array::min(z), array::max(z));
        }

        // array dimensions for non-uniform points
        [[nodiscard]] dims_t dims() const { return x.dims(); }

        // size of the array
        [[nodiscard]] size_t size() const { return x.size(); }

        // theta values
        [[nodiscard]] const std::vector<T> &angles() const { return theta; }

        // get theta value for a given index
        [[nodiscard]] T angle(size_t i) const { return theta[i]; }

        // number of angles
        [[nodiscard]] size_t nprojs() const { return theta.size(); }
    };

} // namespace tomocam

#endif // POLAR_GRID_H
