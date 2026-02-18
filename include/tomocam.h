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
#include "config.h"
#include "dtypes.h"
#include "padding.h"
#include "polar_grid.h"
#include "projection.h"
#include "tiff.h"
#include "timer.h"
#include "write_vti.h"

namespace tomocam {
    /// Restrict T to floating-point types
    template <typename T>
    concept Float = std::is_floating_point<T>::value;

    /**
     * @brief Composition of backprojection and projection operators defined as a
     * system matrix.
     *
     * @param x 3D vector field represented as an array of three components.
     * @param grid The polar grid defining the projection geometry.
     * @param gamma Sample orientation in plane normal to beam direction.
     * @return  \f$ A x = R^T R x \f$
     */
    template <typename Float>
    std::array<Array<Float>, 3> sysmat(const std::array<Array<Float>, 3> &x,
                                       const PolarGrid<Float> &grid, Float gamma);

    /**
     * @brief Computes the gradient of the objective function for iterative
     * reconstruction.
     *
     * @param m Current estimate of the 3-d vector field
     * @param pT backprojection data as an Array.
     * @param grid The polar grid defining the projection geometry.
     * @param gamma Sample orientation in plane normal to beam direction.
     * @return The gradient as an Array.
     */
    template <typename Float>
    std::array<Array<Float>, 3> gradient(const std::array<Array<Float>, 3> &m,
                                         const std::array<Array<Float>, 3> &pT,
                                         const PolarGrid<Float> &grid, Float gamma);

    /**
     * @brief Computes the residual between the projected data and the measured data.
     *
     * @param m Current estimate of the 3-d vector field
     * @param pT backprojected data as split into its components.
     * @param grid The polar grid defining the projection geometry.
     * @param pTp Precomputed inner product of the measured data.
     * @param gamma Sample orientation in plane normal to beam direction.
     * @return The computed residual as a scalar value of type T.
     */
    template <typename Float>
    Float residual(const std::array<Array<Float>, 3> &m,
                   const std::array<Array<Float>, 3> &pT,
                   const PolarGrid<Float> &grid, Float pTp, Float gamma);

    /**
     * @brief Performs Model-Based Iterative Reconstruction (MBIR) for vector
     * tomographic imaging.
     *
     * Reconstructs 3D volumetric data from set of projection images using a
     * regularized optimization. There are two options for regularization: 1).
     * Split-Bregram method with l1-regularization, 2). q-generalized Gaussian Markov
     * Random Field (QGGMRF) edge-preserving The algorithm operates in a padded
     * domain (sqrt(2) factor) to avoid aliasing artifacts during backprojection,
     * then crops the result to the desired dimensions. The reconstruction uses
     * Nesterov's Accelerated Gradient (NAG) method for optimization.
     *
     * @tparam T Floating-point type (float or double).
     * @param datasets A vector of tuples, where each tuple contains:
     *   1). An Array representing projection data
     *   2). A vector of representing the projection angles for the dataset
     *   3). A orientation angle (gamma), a proxy for angle between the sample
     * magnetization and the beam polarization direction.
     * @param recon_dims The dimensions of the reconstructed volume as a tuple
     * (nslice, nrow, ncol).
     * @param params Reconstruction parameters including regularization type, and
     * optimization parameters.
     * @return A array of three components representing the reconstructed 3D vector
     * field, cropped to the specified dimensions.
     */
    template <typename Float>
    std::array<Array<Float>, 3>
    MBIR(const std::vector<std::tuple<Array<Float>, std::vector<Float>, Float>>
             &datasets,
         const dims_t &recon_dims, const ReconParams &params);

} // namespace tomocam

#endif // TOMOCAM__H
