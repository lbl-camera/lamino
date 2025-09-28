
#include "array.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#ifndef POLAR_GRID__H
#define POLAR_GRID__H
namespace tomocam {

    template <typename T>
    struct PolarGrid {
        size_t npts;
        Array<T> x;
        Array<T> y;
        Array<T> z;

        // array dimensions for non-uniform points
        [[nodiscard]] dims_t dims() const { return x.dims(); }

        PolarGrid(const std::vector<T> theta, uint64_t nrows, uint64_t ncols,
            T rotation) {
            dims_t dims = dims_t{theta.size(), nrows, ncols};
            x = Array<T>(dims);
            y = Array<T>(dims);
            z = Array<T>(dims);
            npts = dims.size();

            // rotation matrix
            T cos_t = std::cos(rotation);
            T sin_t = std::sin(rotation);
            T dx = (2 * M_PI) / static_cast<T>(dims.z() - 1);
            T dr = (2 * M_PI) / static_cast<T>(dims.y() - 1);

#pragma omp parallel for collapse(3)
            for (uint64_t i = 0; i < dims.x(); ++i) {
                for (uint64_t j = 0; j < dims.y(); ++j) {
                    for (uint64_t k = 0; k < dims.z(); ++k) {
                        T r_x = k * dx - M_PI;
                        T radius = j * dr - M_PI;
                        T r_y = radius * std::sin(theta[i]);
                        x[i, j, k] = r_x * cos_t - r_y * sin_t;
                        y[i, j, k] = r_x * sin_t + r_y * cos_t;
                        z[i, j, k] = radius * std::cos(theta[i]);
                    }
                }
            }
        }
    };

} // namespace tomocam

#endif // POLAR_GRID__H
