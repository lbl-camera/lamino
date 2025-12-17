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

#ifndef TOMOCAM__H
#define TOMOCAM__H

#include <vector>

#include "array.h"
#include "dtypes.h"
#include "optimize.h"
#include "polar_grid.h"
#include "projection.h"
#include "tiff.h"

namespace tomocam {
    /**
     * @brief Adjoint of the Radon operator for polar grid data.
     * @param projections The input polar grid data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @param dims The dimensions of the output volume.
     * @param gamma orientation of the polar grid (default is 0).
     * @return The adjoint projected data as an Array.
     */
    template <typename T>
    Array<T> backproj(const Array<T> &projections, const PolarGrid<T> &grid,
                      const dims_t &dims, T gamma = T(0)) {
        std::string filter_type = "";
        return backward(projections, grid, dims, gamma, false, filter_type);
    }
    /**
     * @brief Filtered-backprojection of polar grid data to reconstruct a 3D volume.
     * @param projections The input polar grid data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @param dims The dimensions of the output volume.
     * @param gamma orientation of the polar grid (default is 0).
     * @pram filter_type The type of filter to apply (default is "ram-lak").
     * @return The reconstructed volume data as an Array.
     */
    template <typename T>
    Array<T> fbp(const Array<T> &projections, const PolarGrid<T> &grid,
                 const dims_t &dims, const std::string &filter_type = "ramp",
                 T gamma = T(0)) {
        return backward(projections, grid, dims, gamma, true, filter_type);
    }
    /**
     * @brief Function equivalent of system matrix y = A*x
     *
     * @param x Input volume data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @return The equivalent of A^T(A*x) result as an Array.
     */
    template <typename T>
    Array<T> sysmat(const Array<T> &x, const PolarGrid<T> &grid);

    /**
     * @brief Computes the gradient of the objective function for iterative
     * reconstruction.
     *
     * @param x Current estimate of the volume data as an Array.
     * @param b backprojection data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @return The gradient as an Array.
     */
    template <typename T>
    Array<T> gradient(const Array<T> &x, const Array<T> &b, const PolarGrid<T> &grid,
                      T gamma);

    /**
     * @brief Computes the residual between the projected data and the measured data.
     *
     * @param x Current estimate of the volume data as an Array.
     * @param b backprojection data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @param yTy Precomputed inner product of the measured data.
     * @return The computed residual as a scalar value of type T.
     */
    template <typename T>
    T residual(const Array<T> &x, const Array<T> &b, const PolarGrid<T> &grid,
               T yTy);

    /**
     * @brief Performs Model-Based Iterative Reconstruction (MBIR) for vector
     * tomographic imaging.
     *
     * Reconstructs 3D volumetric data from projection images using a regularized
     * optimization approach with QGGMRF (q-Generalized Gaussian Markov Random Field)
     * penalty. The algorithm operates in a padded domain (sqrt(2) factor) to avoid
     * aliasing artifacts during backprojection, then crops the result to the desired
     * dimensions. The reconstruction uses Nesterov's Accelerated Gradient (NAG)
     * method for optimization.
     *
     * The QGGMRF penalty provides edge-preserving regularization, with the parameter
     * p controlling the degree of edge preservation.
     *
     * @tparam T Floating-point type (float or double).
     * @param proj Input projection data as an Array containing measurements from
     * multiple viewing angles.
     * @param angles Vector of projection angles corresponding to each projection in
     * radians (or degrees depending on convention).
     * @param gamma Laminography angle parameter controlling the tilt geometry of the
     * imaging system. For standard tomography, gamma = 0.
     * @param recon_dims Desired output dimensions (n1, n2, n3) for the reconstructed
     * volume. The algorithm internally uses padded dimensions.
     * @param max_iter Maximum number of NAG optimization iterations.
     * @param sigma Regularization strength parameter for QGGMRF penalty. Larger
     * values increase smoothing (default is 1000).
     * @param p Power parameter for QGGMRF penalty controlling edge preservation
     * characteristics (default is 1.2).
     * @param tol Convergence tolerance for the loss function (residual norm).
     * Algorithm stops if relative change in loss falls below this threshold.
     * @param xtol Convergence tolerance for solution updates. Algorithm stops if
     * relative change in reconstruction falls below this threshold.
     a @return std::array<Array<T>, 3> containing three reconstructed volumes
     * representing the three spatial components of magnetization or multi-channel
     * data. Each component is an Array with dimensions specified by recon_dims.
     */
    template <typename T>
    std::array<Array<T>, 3> MBIR(const Array<T> &proj, const std::vector<T> &angles,
                                 T gamma, const dims_t &recon_dims, size_t max_iter,
                                 T sigma, T p, T tol, T xtol);

} // namespace tomocam

#endif // TOMOCAM__H
