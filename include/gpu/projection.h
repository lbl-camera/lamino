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

#ifndef TOMOCAM_GPU_PROJECTION_H
#define TOMOCAM_GPU_PROJECTION_H

#include <array>

#include "dtypes.h"
#include "gpu/device_array.h"
#include "gpu/polar_grid.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu {

    /**
     * @brief GPU forward projection: volume -> projections.
     * @param magnetization 3D magnetization components as std::array of 3
     * DeviceArrays.
     * @param pg      GPU polar grid.
     * @return        Device array of projections with shape (ntheta, nrows, ncols).
     */
    template <typename T>
    DeviceArray<T> forward(const VecArray<T> &magnetization, const PolarGrid<T> &pg,
                           T gamma);

    /**
     * @brief GPU adjoint projection: projections -> magnetization.
     * @param proj        Device array of projections.
     * @param pg          GPU polar grid.
     * @param recon_dims  Dimensions of the output
     * @return            array of 3 DeviceArrays with shape recon_dims,
     * corresponding to the 3 magnetization components.
     */
    template <typename T>
    VecArray<T> adjoint(const DeviceArray<T> &proj, const PolarGrid<T> &pg,
                        const dims_t &recon_dims, T gamma);

    /**
     * @brief GPU system matrix: (A^T A) x.
     * @param x   DeviceArray[3] with shape (nz, ny, nx).
     * @param pg  GPU polar grid.
     * @param gamma  Orientation in plane normal to beam direction
     * @return    DeviceArray[3] of shape (nz, ny, nx).
     */
    template <typename T>
    VecArray<T> sysmat(const VecArray<T> &x, const PolarGrid<T> &pg, T gamma);
} // namespace tomocam::gpu

#endif // TOMOCAM_GPU_PROJECTION_H
