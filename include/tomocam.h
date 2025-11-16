// clang-format off
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
 //clang-format on

#ifndef TOMOCAM__H
#define TOMOCAM__H

#include <vector>

#include "array.h"
#include "dtypes.h"
#include "optimize.h"
#include "polar_grid.h"

namespace tomocam {
    /**
     * @brief Performs a forward projection of 3D volume data onto a polar grid.
     * @param volume The input 3D volume data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @param gamma orientation of the polar grid (default is 0).
     * @return The projected data as an Array.
     */
    template <typename T>
    Array<T> forward(const Array<T> &, const PolarGrid<T> &, T gamma = T(0));

    /**
     * @brief Performs a backward projection from polar grid data to reconstruct a 3D
     * volume.
     * @param projections The input polar grid data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @param dims The dimensions of the output volume.
     * @param gamma orientation of the polar grid (default is 0).
     * @param filter Whether to apply a filter during backprojection (default is
     * false).
     * @return The backprojected polar data as an Array
     */
    template <typename T>
    Array<T> backward(const Array<T> &, const PolarGrid<T> &, const dims_t &,
                      T gamma = T(0), bool filter = false);

    /**
     * @brief Function equivalent of system matrix y = A*x
     * @param x Input volume data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @return The equivalent of A^T(A*x) result as an Array.
     */
    template <typename T>
    Array<T> sysmat(const Array<T> &, const PolarGrid<T> &);

    /**
     * @brief Computes the gradient of the objective function for iterative
     * reconstruction.
     * @param x Current estimate of the volume data as an Array.
     * @param b backprojection data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @return The gradient as an Array.
     */
    template <typename T>
    Array<T> gradient(const Array<T> &, const Array<T> &, const PolarGrid<T> &);

    /**
     * @brief Computes the residual between the projected data and the measured data.
     * @param x Current estimate of the volume data as an Array.
     * @param b backprojection data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @param yTy Precomputed inner product of the measured data.
     * @return The computed residual as a scalar value of type T.
     */
    template <typename T>
    T residual(const Array<T> &, const Array<T> &, const PolarGrid<T> &, T);

    /**
     * @brief Performs Model-Based Iterative Reconstruction (MBIR) of volume data.
     * @param theta std::vector containing the projection angles.
     * @param recon_dims Dimensions of the output reconstructed volume.
     * @param max_iter Maximum number of iterations for the optimization.
     * @param sigma Parameter for the QGGMRF penalty function.
     * @param p Parameter for the QGGMRF penalty function.
     * @return The reconstructed volume data as an Array.
     */
    template <typename T>
    Array<T> MBIR(const Array<T> &, std::vector<T>, dims_t, size_t, T, T, T, T);

} // namespace tomocam

#endif // TOMOCAM__H
