#include <vector>

#include "array.h"
#include "dtypes.h"
#include "polar_grid.h"

using std::vector;

#ifndef TOMOCAM__H
    #define TOMOCAM__H

namespace tomocam {
    template <typename T>
    Array<T> forward(const Array<T> &, const PolarGrid<T> &, T);

    template <typename T>
    Array<T> backward(const Array<T> &, const PolarGrid<T> &, T, dims_t, bool);

    namespace opt {

        /**
         * Conjugate-Gradient method to solve for the linear system (Ax = yT)
         * for the Split-Bregman method.
         *
         * @param A Function that computes the left-hand side of the linear
         * system (Ax).
         * @param yT Right-hand side of the linear system.
         * @param max_iter Maximum number of iterations.
         * @param tol Tolerance for convergence.
         * @return Solution vector x.
         */
        template <typename T>
        Array<T> cgsolver(std::function<Array<T>(const Array<T> &)>,
            const Array<T> &, size_t, T);

        /**
         * Split-Bregman method for TV-regularized tomographic reconstruction
         * Solves: min_x ||R(x) - y||^2 + λ ||grad(x)||_1
         *
         * The problem is split into:
         * - min_x ||R(x) - y||^2 + (mu/2)||grad(x) - d + b||^2
         * - min_d λ||d||_1 + (mu/2)||grad(x) - d + b||^2
         *
         * @param function composed of A x = R^T R x => g = A x - yT
         * @param yT: Adjoint_R on y
         * @param lambda: TV regularization parameter
         * @param mu: Penalty parameter for Split-Bregman (default: 1.0)
         * @param outer_iter: Number of outer Bregman iterations (default: 50)
         * @param inner_iter: Number of inner iterations for x-update (default:
         * 10)
         * @param tol: Convergence tolerance (default: 1e-4)
         */
        template <typename T>
        Array<T> split_bregman(std::function<Array<T>(const Array<T> &)>,
            const Array<T> &, T, T, size_t, size_t, T);
    } // namespace opt

} // namespace tomocam

#endif // TOMOCAM__H
