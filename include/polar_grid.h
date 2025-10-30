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
            T dx = (2 * M_PI) / static_cast<T>(dims.z() - 1);
            T dr = (2 * M_PI) / static_cast<T>(dims.y() - 1);

    #pragma omp parallel for collapse(3)
            for (size_t i = 0; i < dims.x(); ++i) {
                for (size_t j = 0; j < dims.y(); ++j) {
                    for (size_t k = 0; k < dims.z(); ++k) {
                        T radius = j * dr - M_PI;
                        x[{i, j, k}] = k * dx - M_PI;
                        y[{i, j, k}] = radius * sin(theta[i]);
                        z[{i, j, k}] = radius * std::cos(theta[i]);
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
