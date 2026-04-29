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

#ifndef PROJECTION__H
#define PROJECTION__H

#include <array>
#include <cmath>

#include "array.h"
#include "dtypes.h"
#include "polar_grid.h"

namespace tomocam {

    // Returns the beam direction vector for a laminography geometry,
    // giving the projection coefficient for each magnetization component.
    // Components: {sin(gamma)*sin(alpha), cos(alpha), cos(gamma)*sin(alpha)}.
    template <typename T>
    inline std::array<T, 3> beam_dir_vector(T gamma, T alpha) {
        return {std::sin(gamma) * std::sin(alpha),
                std::cos(alpha),
                std::cos(gamma) * std::sin(alpha)};
    }

    /**
     * @brief Performs a forward projection of 3D magnetization data onto a polar
     * grid.
     * @param magnetization The input 3D magnetization components as std::array of 3
     * Arrays.
     * @param pg The polar grid defining the projection geometry.
     * @param gamma orientation of the polar grid.
     * @return The projected data as an Array.
     */
    template <typename T>
    Array<T> forward(const std::array<Array<T>, 3> &magnetization,
                     const PolarGrid<T> &pg, const T &gamma);

    /**
     * @brief Adjoint of the forward projection operator for polar grid data.
     * @param proj The input polar grid data as an Array.
     * @param pg The polar grid defining the projection geometry.
     * @param recon_dims The dimensions of the output volume.
     * @param gamma orientation of the polar grid.
     * @return The adjoint projected data as std::array of 3 Arrays.
     */
    template <typename T>
    std::array<Array<T>, 3> adjoint(const Array<T> &proj, const PolarGrid<T> &pg,
                                    const dims_t &recon_dims, T gamma);

} // namespace tomocam

#endif // PROJECTION__H
