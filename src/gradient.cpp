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
#include <complex>
#include <execution>
#include <format>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "nufft.h"
#include "polar_grid.h"
#include "projection.h"
#include "tomocam.h"

namespace tomocam {
    template <typename T>
    std::array<Array<T>, 3> sysmat(const std::array<Array<T>, 3> &x,
                                   const PolarGrid<T> &grid, T beta) {
        using complex_t = std::complex<T>;

        T scale = static_cast<T>(grid.dims().n2 * grid.dims().n3);

        std::array<Array<complex_t>, 3> c_components;
        for (size_t i = 0; i < 3; ++i) {
            auto x_cmplx = array::to_complex(x[i]);
            c_components[i] = Array<complex_t>::zeros(grid.dims());
            nufft::nufft3d2(c_components[i], x_cmplx, grid);
        }

        std::array<Array<complex_t>, 3> result_components;
        for (size_t i = 0; i < 3; ++i) {
            result_components[i] = Array<complex_t>::zeros(grid.dims());
        }

        for (size_t j = 0; j < grid.nprojs(); ++j) {
            auto coeff = beam_dir_vector(grid.angle(j), beta);

            // Matrix multiplication: result = coeff.T * coeff * c_components
            // This is outer product of coeff with itself, applied to c_components
            for (size_t i = 0; i < 3; ++i) {
                auto result_slice = result_components[i].slice(j, j + 1);
                for (size_t k = 0; k < 3; ++k) {
                    T weight = coeff[i] * coeff[k];
                    auto c_slice = c_components[k].slice(j, j + 1);
                    std::transform(
                        std::execution::par_unseq, c_slice.begin(), c_slice.end(),
                        result_slice.begin(), result_slice.begin(),
                        [weight](const complex_t &c_val, complex_t &r_val) {
                            return r_val + weight * c_val;
                        });
                }
            }
        }

        // Step 3: Apply nufft3d1 to each resulting component
        std::array<Array<T>, 3> output;
        for (size_t i = 0; i < 3; ++i) {
            auto out_cmplx = Array<complex_t>(x[i].dims());
            nufft::nufft3d1(result_components[i], out_cmplx, grid);
            output[i] = array::to_real(out_cmplx) / scale;
        }
        return output;
    }
    template std::array<Array<float>, 3> sysmat(const std::array<Array<float>, 3> &x,
                                                const PolarGrid<float> &grid,
                                                float beta);
    template std::array<Array<double>, 3>
    sysmat(const std::array<Array<double>, 3> &x, const PolarGrid<double> &grid,
           double beta);

    // Overload: beta read per-angle from grid.beta(j)
    template <typename T>
    std::array<Array<T>, 3> sysmat(const std::array<Array<T>, 3> &x,
                                   const PolarGrid<T> &grid) {
        using complex_t = std::complex<T>;

        T scale = static_cast<T>(grid.dims().n2 * grid.dims().n3);

        std::array<Array<complex_t>, 3> c_components;
        for (size_t i = 0; i < 3; ++i) {
            auto x_cmplx = array::to_complex(x[i]);
            c_components[i] = Array<complex_t>::zeros(grid.dims());
            nufft::nufft3d2(c_components[i], x_cmplx, grid);
        }

        std::array<Array<complex_t>, 3> result_components;
        for (size_t i = 0; i < 3; ++i) {
            result_components[i] = Array<complex_t>::zeros(grid.dims());
        }

        for (size_t j = 0; j < grid.nprojs(); ++j) {
            auto coeff = beam_dir_vector(grid.angle(j), grid.beta(j));

            for (size_t i = 0; i < 3; ++i) {
                auto result_slice = result_components[i].slice(j, j + 1);
                for (size_t k = 0; k < 3; ++k) {
                    T weight = coeff[i] * coeff[k];
                    auto c_slice = c_components[k].slice(j, j + 1);
                    std::transform(
                        std::execution::par_unseq, c_slice.begin(), c_slice.end(),
                        result_slice.begin(), result_slice.begin(),
                        [weight](const complex_t &c_val, complex_t &r_val) {
                            return r_val + weight * c_val;
                        });
                }
            }
        }

        std::array<Array<T>, 3> output;
        for (size_t i = 0; i < 3; ++i) {
            auto out_cmplx = Array<complex_t>(x[i].dims());
            nufft::nufft3d1(result_components[i], out_cmplx, grid);
            output[i] = array::to_real(out_cmplx) / scale;
        }
        return output;
    }
    template std::array<Array<float>, 3>
    sysmat(const std::array<Array<float>, 3> &x, const PolarGrid<float> &grid);
    template std::array<Array<double>, 3>
    sysmat(const std::array<Array<double>, 3> &x, const PolarGrid<double> &grid);

} // namespace tomocam
