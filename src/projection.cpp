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

#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <execution>
#include <functional>
#include <stdexcept>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "fft.h"
#include "fftutils.h"
#include "filter.h"
#include "nufft.h"
#include "padding.h"
#include "polar_grid.h"
#include "tomocam.h"

namespace tomocam {

    template <typename T>
    Array<T> forward(const std::array<Array<T>, 3> &magnetization,
                     const PolarGrid<T> &pg, const T &gamma) {

        // nomralization factor
        auto dims = pg.dims();
        T scale = static_cast<T>(dims.n2 * dims.n3);

        // initialize output array
        using complex_t = std::complex<T>;
        auto proj = Array<complex_t>::zeros(pg.dims());

        // projection vector
        std::array<T, 3> c_gamma = {std::cos(gamma), std::sin(gamma), 1.0};
        std::array<std::function<T(T)>, 3> c_alpha = {
            [](T alpha) { return std::sin(alpha); },
            [](T alpha) { return -std::sin(alpha); },
            [](T alpha) { return std::cos(alpha); }};

        // loop over magnetization components
        for (size_t i = 0; i < 3; ++i) {
            // cast array to complex
            auto m_cmplx = array::to_complex(magnetization[i]);
            auto c_cmplx = Array<complex_t>::zeros(pg.dims());

            // call NUFFT3d type-2
            nufft::nufft3d2<T>(c_cmplx, m_cmplx, pg);

            // add to projection scaled by coeff
            for (size_t j = 0; j < pg.nprojs(); ++j) {
                auto slice = c_cmplx.slice(j, j + 1);
                T coeff = c_gamma[i] * c_alpha[i](pg.angles()[j]);
                std::for_each(std::execution::par_unseq, slice.begin(), slice.end(),
                              [coeff](complex_t &val) { val *= coeff; });
            }
            proj += c_cmplx;
        }

        // ifft
        proj = fft::fftshift2(proj);
        proj = fft::ifft2(proj);
        proj = fft::ifftshift2(proj);
        return array::to_real<T>(proj) / scale;
    }
    // Explicit instantiation forward
    template Array<float>
    forward<float>(const std::array<Array<float>, 3> &magnetization,
                   const PolarGrid<float> &pg, const float &gamma);
    template Array<double>
    forward<double>(const std::array<Array<double>, 3> &magnetization,
                    const PolarGrid<double> &pg, const double &gamma);

    /*
     * Backprojection operation
     */
    template <typename T>
    std::array<Array<T>, 3> adjoint(const Array<T> &proj, const PolarGrid<T> &pg,
                                    const dims_t &recon_dims, T gamma) {

        // cast to complex
        auto c_cmplx = array::to_complex(proj);

        // 2-D Fourier transforms
        c_cmplx = fft::fftshift2(c_cmplx);
        c_cmplx = fft::fft2(c_cmplx);
        c_cmplx = fft::ifftshift2(c_cmplx);

        std::array<T, 3> c_gamma = {std::cos(gamma), std::sin(gamma), 1.0};
        std::array<std::function<T(T)>, 3> c_alpha = {
            [](T alpha) { return std::sin(alpha); },
            [](T alpha) { return -std::sin(alpha); },
            [](T alpha) { return std::cos(alpha); }};

        // nufft - each component separately
        std::array<Array<T>, 3> m_components;
        using complex_t = std::complex<T>;

        for (size_t i = 0; i < 3; ++i) {
            auto c_cmplx_copy = c_cmplx.clone();

            // scale by coefficients for this component
            for (size_t j = 0; j < pg.nprojs(); ++j) {
                T coeff = c_gamma[i] * c_alpha[i](pg.angles()[j]);
                auto slice = c_cmplx_copy.slice(j, j + 1);
                std::for_each(std::execution::par_unseq, slice.begin(), slice.end(),
                              [coeff](complex_t &val) { val *= coeff; });
            }

            // apply NUFFT for this component
            Array<complex_t> m_cmplx(recon_dims);
            nufft::nufft3d1<T>(c_cmplx_copy, m_cmplx, pg);
            m_components[i] = std::move(array::to_real<T>(m_cmplx));
        }

        return m_components;
    }
    // Explicit instantiation adjoint
    template std::array<Array<float>, 3> adjoint<float>(const Array<float> &proj,
                                                        const PolarGrid<float> &pg,
                                                        const dims_t &recon_dims,
                                                        float gamma);
    template std::array<Array<double>, 3>
    adjoint<double>(const Array<double> &proj, const PolarGrid<double> &pg,
                    const dims_t &recon_dims, double gamma);

} // namespace tomocam
