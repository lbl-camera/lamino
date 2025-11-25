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

#include "array.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace tomocam {

    template <typename T>
    struct PolarGrid {
        size_t npts;
        Array<T> x;
        Array<T> y;
        Array<T> z;

        // array dimensions for non-uniform points
        [[nodiscard]] dims_t dims() const { return x.dims(); }

        // size of the array
        [[nodiscard]] size_t size() const { return x.size(); }

        // constructor
        PolarGrid(const std::vector<T> &theta, size_t nrows, size_t ncols) {

            dims_t dims = dims_t{theta.size(), nrows, ncols};
            x = Array<T>(dims);
            y = Array<T>(dims);
            z = Array<T>(dims);
            npts = dims.size();

            // rotation matrix
            T dz = (2 * M_PI) / static_cast<T>(nrows - 1);
            T dr = (2 * M_PI) / static_cast<T>(ncols - 1);

#pragma omp parallel for collapse(3)
            for (size_t i = 0; i < dims.n1; ++i) {
                for (size_t j = 0; j < dims.n2; ++j) {
                    for (size_t k = 0; k < dims.n3; ++k) {
                        T radius = j * dr - M_PI;
                        y[{i, j, k}] = radius * std::cos(theta[i]);
                        z[{i, j, k}] = radius * std::sin(theta[i]);
                        x[{i, j, k}] = k * dz - M_PI;
                    }
                }
            }
        }

        PolarGrid<T> rotate(T angle) const {

            PolarGrid<T> out = *this;

            T cos_t = std::cos(angle);
            T sin_t = std::sin(angle);

            auto dims = this->dims();
            for (size_t i = 0; i < x.size(); i++) {
                out.x[i] = x[i] * cos_t - y[i] * sin_t;
                out.y[i] = x[i] * sin_t + y[i] * cos_t;
            }
            return out;
        }
    };

} // namespace tomocam

#endif // POLAR_GRID__H
