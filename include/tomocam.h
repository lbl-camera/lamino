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

#ifndef TOMOCAM_H
#define TOMOCAM_H

#include <tuple>
#include <vector>

#include "array.h"
#include "recon_params.h"

namespace tomocam {

    /**
     * @brief Performs Model-Based Iterative Reconstruction (MBIR) for vector
     * tomography using unconstrained Conjugate Gradient (CG) optimization.
     *
     * @tparam T Floating-point type (float or double).
     * @param datasets A vector of tuples, where each tuple contains:
     *   1). A tomocam::Array representing projection data
     *   2). An std::vector of representing the projection angles for the dataset
     *   3). Orientation angle (gamma), rotation around the beam axis
     * @param recon_dims The dimensions of the reconstructed volume as a dims_t
     * object.
     * @param params Reconstruction parameters including regularization type, and
     * optimization parameters.
     * @return A array of three components representing the reconstructed 3D vector
     * field, cropped to the specified dimensions.
     */
    template <typename T>
    std::array<Array<T>, 3>
    MBIR(const std::vector<std::tuple<Array<T>, std::vector<T>, T>> &datasets,
         const ReconParams &params);

} // namespace tomocam

#endif // TOMOCAM_H
