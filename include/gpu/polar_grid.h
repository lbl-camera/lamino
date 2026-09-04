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

#ifndef TOMOCAM_GPU_POLAR_GRID_H
#define TOMOCAM_GPU_POLAR_GRID_H

#include <cstddef>
#include <vector>

#include <thrust/device_vector.h>

#include "dtypes.h"
#include "gpu/device_array.h"
#include "gpu/device_ptr.h"
#include "gpu/utils.h"

namespace tomocam::gpu {

    /// Mirrors the interface of tomocam::PolarGrid but stores coordinate
    /// arrays in device memory as DeviceArray<T>.
    ///
    /// Construction uploads the computed (x, y, z) grid points directly
    /// to the GPU via a CUDA kernel. Per-projection theta/gamma/beta are
    /// kept host-side (mirroring CPU PolarGrid) since Toeplitz PSF weight
    /// computation (beam_dir_vector per row) is done on the host; gamma/beta
    /// are additionally uploaded so the grid-construction and projection
    /// kernels can read them per-row on the device.
    ///
    /// @tparam T Floating-point element type (float or double)
    template <typename T>
    struct PolarGrid {
        size_t npts;
        std::vector<T> theta;   // host, length nprojs
        std::vector<T> gammas;  // host, length nprojs (broadcast-constant per dataset)
        std::vector<T> betas;   // host, length nprojs
        DeviceArray<T> x;
        DeviceArray<T> y;
        DeviceArray<T> z;
        DeviceArray<T> w; // 1 inside [-pi,pi]^3, 0 outside (masks aliased q-points)
        thrust::device_vector<T> angles;
        thrust::device_vector<T> d_gammas;
        thrust::device_vector<T> d_betas;

        /// Returns dimensions (nangles, nrows, ncols) of the coordinate arrays.
        [[nodiscard]] dims_t dims() const { return x.dims(); }

        /// Returns total number of non-uniform grid points.
        [[nodiscard]] size_t size() const { return x.size(); }

        /// Returns theta value for a given projection index.
        [[nodiscard]] T angle(size_t i) const { return theta[i]; }

        /// Returns gamma value for a given projection index.
        [[nodiscard]] T gamma(size_t i) const { return gammas[i]; }

        /// Returns beta value for a given projection index.
        [[nodiscard]] T beta(size_t i) const { return betas[i]; }

        /// Returns number of projections.
        [[nodiscard]] size_t nprojs() const { return theta.size(); }

        /// Constructs the polar grid on the GPU.
        ///
        /// @param theta  Host-side vector of projection angles (radians)
        /// @param nrows  Number of radial samples
        /// @param ncols  Number of axial samples
        /// @param gamma  Out-of-plane tilt angle (radians), broadcast to all angles
        /// @param beta   Out-of-plane tilt angle (radians), broadcast to all angles
        PolarGrid(const std::vector<T> &theta, size_t nrows, size_t ncols, T gamma,
                 T beta = T(0));

        // default constructor
        PolarGrid()
            : npts(0), x(), y(), z(), w(), angles(), d_gammas(), d_betas() {}
    };
} // namespace tomocam::gpu

#endif // TOMOCAM_GPU_POLAR_GRID_H
