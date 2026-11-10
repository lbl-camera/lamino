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

#ifndef TOMOCAM_GPU_H
#define TOMOCAM_GPU_H

#include <string>
#include <vector>

#include "gpu/device_array.h"
#include "gpu/polar_grid.h"
#include "gpu/projection.h"
#include "recon_params.h"

namespace tomocam::gpu {

    // reuse the CPU host-side struct verbatim: it's plain host memory
    // (Array<T>, std::vector<T>, scalars), and MBIR's callers already build
    // one of these directly (see src/recon.cpp).
    template <typename T>
    using Dataset_t = tomocam::Dataset_t<T>;

    /**
     * @brief Host wrapper for GPU MBIR: takes projections and angles on the host,
     *        moves data to GPU, calls the GPU MBIR function, and returns the
     *        reconstructed volume back on the host.
     * @param dataset  Dataset containing projections, angles, and gamma reference.
     * @param params   Reconstruction parameters (e.g. number of iterations,
     * regularization weight, etc.).
     * @param mask     Binary support (dims = params.recon_dims) for the
     *                 const_mag regularizers, see support_mask() in mask.h
     * @return         Reconstructed vector components
     * similar).
     */
    template <typename T>
    std::array<Array<T>, 3> MBIR(const std::vector<Dataset_t<T>> &datasets,
                                 const ReconParams &params, const Array<T> &mask);

} // namespace tomocam::gpu

#endif // TOMOCAM_GPU_H
