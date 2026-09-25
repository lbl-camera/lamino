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
#ifndef GPU_PADDING_H
#define GPU_PADDING_H

#include "dtypes.h"
#include "gpu/device_array.h"

namespace tomocam::gpu {

    enum class PadType { LEFT, RIGHT, SYMMETRIC };

    template <typename T>
    size_t n_pad(size_t n, T factor) {
        size_t n2 = static_cast<size_t>(n * factor);
        return 2 * ((n2 - n) / 2);
    }

    // Odd padded size shared by volume and detector grid; mirrors CPU
    // tomocam::padded_dim (include/padding.h), see there for why it is odd.
    template <typename T>
    size_t padded_dim(size_t n, T factor) {
        size_t N = n + n_pad(n, factor);
        return (N % 2 == 0) ? N - 1 : N;
    }

    /**
     * @tparam T
     * @param input
     * @param padding
     * @param type
     * @return A new DeviceArray with the specified padding applied to the input
     * array. The type of padding is determined by the PadType enum.
     *
     */
    template <typename T>
    DeviceArray<T> pad2d(const DeviceArray<T> &input, float factor, PadType type);

    /**
     * @brief zero-pad `input` into a `new_dims`-sized array, placing `input`
     * starting at `offset` along each axis (e.g. offset={0,0,0} = zero
     * right-pad). Mirrors CPU's pad3d(arr, new_dims, offset) (include/padding.h).
     */
    template <typename T>
    DeviceArray<T> pad3d(const DeviceArray<T> &input, dims_t new_dims,
                         dims_t offset);

    /**
     * @brief zero-pad `input` by `factor` in every dimension, per `type`.
     */
    template <typename T>
    DeviceArray<T> pad3d(const DeviceArray<T> &input, float factor, PadType type);

    /**
     * @tparam T
     * @param input
     * @param output dimensions
     * @param type
     * @return A new DeviceArray with the specified padding applied to the input
     * array. The type of padding is determined by the PadType enum.
     *
     */
    template <typename T>
    DeviceArray<T> crop3d(const DeviceArray<T> &input, dims_t out_dims,
                          PadType type);

    /**
     * @brief crop a `new_dims`-sized block out of `input`, reading starting at
     * `offset` along each axis. Mirrors CPU's crop3d(arr, new_dims, offset)
     * (include/padding.h).
     */
    template <typename T>
    DeviceArray<T> crop3d(const DeviceArray<T> &input, dims_t new_dims,
                          dims_t offset);

} // namespace tomocam::gpu

#endif // GPU_PADDING_H
