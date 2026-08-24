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
#include <vector>

#include "array.h"
#include "dtypes.h"
#include "polar_grid.h"
#include "rotation.h"

namespace tomocam {

    // Third row of RotationTranspose(theta, gamma, beta) = R^T[2].
    // R = Rx(theta)*Ry(beta)*Rz(gamma); gamma does not appear in this row.
    template <typename T>
    inline std::array<T, 3> beam_dir_vector(T theta, T beta) {
        return {std::sin(beta),
                -std::sin(theta) * std::cos(beta),
                 std::cos(theta) * std::cos(beta)};
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
                     const PolarGrid<T> &pg, T beta);

    // adjoint with explicit beta; optional per-projection COR shifts
    template <typename T>
    std::array<Array<T>, 3>
    adjoint(const Array<T> &proj, const PolarGrid<T> &pg, const dims_t &recon_dims,
            T beta, const std::vector<std::array<T, 2>> &shifts = {});

    // A^T A x with explicit beta
    template <typename T>
    std::array<Array<T>, 3> sysmat(const std::array<Array<T>, 3> &x,
                                   const PolarGrid<T> &grid, T beta);

    // Overloads for unified PolarGrid — beta read per-angle from pg.beta(j)

    template <typename T>
    std::array<Array<T>, 3>
    adjoint(const Array<T> &proj, const PolarGrid<T> &pg, const dims_t &recon_dims);

    template <typename T>
    std::array<Array<T>, 3> sysmat(const std::array<Array<T>, 3> &x,
                                   const PolarGrid<T> &grid);

} // namespace tomocam

#endif // PROJECTION__H
