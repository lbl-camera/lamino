
#include "array.h"
#include <cmath>
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

        PolarGrid(const std::vector<T> theta, int nrows, int ncols) {
            dims_t dims = {theta.size(), nrows, ncols};
            x(dims);
            y(dims);
            z(dims);
            npts = dims.size();

#pragma omp parallel for collapse(3)
            for (int i = 0; i < dims.x; ++i) {
                for (int j = 0; j < dims.y; ++j) {
                    for (int k = 0; k < dims.z; ++k) {
                        x[i, j, k] = 2 * M_PI / (dims.z - 1) * (k - dims.z / 2);
                        y[i, j, k] = 2 * M_PI / (dims.y - 1) *
                                     (j - dims.y / 2) * std::sin(theta[i]);
                        z[i, j, k] = 2 * M_PI / (dims.y - 1) *
                                     (j - dims.y / 2) * std::cos(theta[i]);
                    }
                }
            }
        }
    };

} // namespace tomocam

#endif // POLAR_GRID__H
