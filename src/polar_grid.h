
#include "array.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#ifndef POLAR_GRID__H
#define POLAR_GRID__H
namespace tomocam {

    template <typename T> struct PolarGrid {
        size_t npts;
        Array<T> x;
        Array<T> y;
        Array<T> z;

        // array dimensions for non-uniform points
        dims_t dims() const { return x.dims(); }

        PolarGrid(const std::vector<T> theta, uint64_t nrows, uint64_t ncols) {
            dims_t dims = dims_t{theta.size(), nrows, ncols};
            x = Array<T>(dims);
            y = Array<T>(dims);
            z = Array<T>(dims);
            npts = dims.size();

#pragma omp parallel for collapse(3)
            for (uint64_t i = 0; i < dims.x(); ++i) {
                for (uint64_t j = 0; j < dims.y(); ++j) {
                    for (uint64_t k = 0; k < dims.z(); ++k) {
                        x[i, j, k] = -M_PI + i * (2 * M_PI) / static_cast<T>(dims.z() - 1);
                        T r = -M_PI + k * (2 * M_PI) / static_cast<T>(dims.y() - 1);
                        y[i, j, k] = r * std::cos(theta[i]);
                        z[i, j, k] = r * std::sin(theta[i]);
                    }
                }
            }
        }
    };

} // namespace tomocam

#endif // POLAR_GRID__H
