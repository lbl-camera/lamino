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

#ifndef TOMOCAM_GPU_PHASE_SHIFT_H
#define TOMOCAM_GPU_PHASE_SHIFT_H

#include <cuda/std/complex>
#include <thrust/device_vector.h>

#include "gpu/device_array.h"

namespace tomocam::gpu {

    template <typename T>
    using Complex = cuda::std::complex<T>;

    // Applies a per-projection center-of-rotation phase shift in-place:
    // data[i,j,k] *= exp(-i*(qx*dx[i] + qy*dy[i])), mirroring
    // tomocam::fft::phase_shift2d (include/fft.h). dx/dy have length nprojs
    // (= data.dims().n1), one entry per projection.
    template <typename T>
    void phase_shift2d(DeviceArray<Complex<T>> &data,
                       const thrust::device_vector<T> &dx,
                       const thrust::device_vector<T> &dy);

} // namespace tomocam::gpu

#endif // TOMOCAM_GPU_PHASE_SHIFT_H
