
#include "array.h"
#include "bregman.h"

namespace tomocam::opt {

    // Forward differences with a zero (Neumann) difference at the last index of
    // each axis:  (D u)[i] = u[i + 1] - u[i] for i < n - 1,  (D u)[n - 1] = 0.
    // divergence() is the exact adjoint D^T:
    //   (D^T p)[i] = p[i - 1] - p[i], with p[-1] = p[n - 1] = 0,
    // so laplacian = D^T D is positive semidefinite.
    // Axis 0, 1, 2 differentiate along slices, rows and columns respectively.

    template <typename T>
    std::array<Array<T>, 3> grad_u(const Array<T> &u) {
        const size_t ns = u.nslices(), nr = u.nrows(), nc = u.ncols();
        std::array<Array<T>, 3> grad;
        for (auto &g : grad) { g = Array<T>::zeros(u.dims()); }
#pragma omp parallel for collapse(3)
        for (size_t slc = 0; slc < ns; ++slc) {
            for (size_t row = 0; row < nr; ++row) {
                for (size_t col = 0; col < nc; ++col) {
                    const T c = u[{slc, row, col}];
                    if (slc + 1 < ns)
                        grad[0][{slc, row, col}] = u[{slc + 1, row, col}] - c;
                    if (row + 1 < nr)
                        grad[1][{slc, row, col}] = u[{slc, row + 1, col}] - c;
                    if (col + 1 < nc)
                        grad[2][{slc, row, col}] = u[{slc, row, col + 1}] - c;
                }
            }
        }
        return grad;
    }
    // explicit instantiation for float and double
    template std::array<Array<float>, 3> grad_u(const Array<float> &u);
    template std::array<Array<double>, 3> grad_u(const Array<double> &u);

    // adjoint of grad_u (negative divergence)
    template <typename T>
    Array<T> divergence(const std::array<Array<T>, 3> &du) {
        const size_t ns = du[0].nslices(), nr = du[0].nrows(), nc = du[0].ncols();
        Array<T> div = Array<T>::zeros(du[0].dims());
#pragma omp parallel for collapse(3)
        for (size_t slc = 0; slc < ns; ++slc) {
            for (size_t row = 0; row < nr; ++row) {
                for (size_t col = 0; col < nc; ++col) {
                    T val = 0;
                    if (slc + 1 < ns) val -= du[0][{slc, row, col}];
                    if (slc > 0) val += du[0][{slc - 1, row, col}];
                    if (row + 1 < nr) val -= du[1][{slc, row, col}];
                    if (row > 0) val += du[1][{slc, row - 1, col}];
                    if (col + 1 < nc) val -= du[2][{slc, row, col}];
                    if (col > 0) val += du[2][{slc, row, col - 1}];
                    div[{slc, row, col}] = val;
                }
            }
        }
        return div;
    }
    // explicit instantiation for float and double
    template Array<float> divergence(const std::array<Array<float>, 3> &du);
    template Array<double> divergence(const std::array<Array<double>, 3> &du);

} // namespace tomocam::opt
