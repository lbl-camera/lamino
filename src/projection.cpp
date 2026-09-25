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
#include "projection.h"

#include "tomocam.h"

namespace tomocam {

    template <typename T>
    Array<T> forward(const std::array<Array<T>, 3> &magnetization,
                     const PolarGrid<T> &pg, T gamma, T beta) {

        auto dims = magnetization[0].dims();
        T scale = static_cast<T>(dims.n1 * dims.n2 * dims.n3);

        using complex_t = std::complex<T>;
        auto proj = Array<complex_t>::zeros(pg.dims());

        for (size_t i = 0; i < 3; ++i) {
            auto m_cmplx = array::to_complex(magnetization[i]);
            auto c_cmplx = Array<complex_t>::zeros(pg.dims());
            nufft::nufft3d2<T>(c_cmplx, m_cmplx, pg);

            // discard aliased points (outside [-pi,pi]) before ifft2
            std::transform(std::execution::par_unseq, c_cmplx.begin(),
                           c_cmplx.end(), pg.w.begin(), c_cmplx.begin(),
                           [](complex_t c, T m) { return c * m; });

            for (size_t j = 0; j < pg.nprojs(); ++j) {
                auto slice = c_cmplx.slice(j, j + 1);
                T coeff = beam_dir_vector(pg.angle(j), gamma, beta)[i];
                std::for_each(std::execution::par_unseq, slice.begin(), slice.end(),
                              [coeff](complex_t &val) { val *= coeff; });
            }
            proj += c_cmplx;
        }

        proj = fft::fftshift2(proj);
        proj = fft::ifft2(proj);
        proj = fft::ifftshift2(proj);
        return array::to_real<T>(proj) / scale;
    }
    template Array<float>
    forward<float>(const std::array<Array<float>, 3> &magnetization,
                   const PolarGrid<float> &pg, float gamma, float beta);
    template Array<double>
    forward<double>(const std::array<Array<double>, 3> &magnetization,
                    const PolarGrid<double> &pg, double gamma, double beta);

    template <typename T>
    std::array<Array<T>, 3> adjoint(const Array<T> &proj, const PolarGrid<T> &pg,
                                    const dims_t &recon_dims, T gamma, T beta,
                                    const std::vector<std::array<T, 2>> &shifts) {

        auto c_cmplx = array::to_complex(proj);
        T scale = static_cast<T>(recon_dims.n1 * recon_dims.n2 * recon_dims.n3);

        c_cmplx = fft::fftshift2(c_cmplx);
        c_cmplx = fft::fft2(c_cmplx);
        c_cmplx = fft::ifftshift2(c_cmplx);

        if (!shifts.empty()) fft::phase_shift2d(c_cmplx, shifts);

        // discard aliased points (outside [-pi,pi]) before backprojecting
        std::transform(std::execution::par_unseq, c_cmplx.begin(), c_cmplx.end(),
                       pg.w.begin(), c_cmplx.begin(),
                       [](std::complex<T> c, T m) { return c * m; });

        std::array<Array<T>, 3> m_components;
        using complex_t = std::complex<T>;

        for (size_t i = 0; i < 3; ++i) {
            auto c_cmplx_copy = c_cmplx.clone();

            for (size_t j = 0; j < pg.nprojs(); ++j) {
                T coeff = beam_dir_vector(pg.angle(j), gamma, beta)[i];
                auto slice = c_cmplx_copy.slice(j, j + 1);
                std::for_each(std::execution::par_unseq, slice.begin(), slice.end(),
                              [coeff](complex_t &val) { val *= coeff; });
            }
            Array<complex_t> m_cmplx(recon_dims);
            nufft::nufft3d1<T>(c_cmplx_copy, m_cmplx, pg);
            m_components[i] = array::to_real<T>(m_cmplx) / scale;
        }
        return m_components;
    }

    // explicit instantiations for float and double
    template std::array<Array<float>, 3>
    adjoint<float>(const Array<float> &proj, const PolarGrid<float> &pg,
                   const dims_t &recon_dims, float gamma, float beta,
                   const std::vector<std::array<float, 2>> &shifts);
    template std::array<Array<double>, 3>
    adjoint<double>(const Array<double> &proj, const PolarGrid<double> &pg,
                    const dims_t &recon_dims, double gamma, double beta,
                    const std::vector<std::array<double, 2>> &shifts);

    // Overload: read gamma and beta read per-angle from pg
    template <typename T>
    std::array<Array<T>, 3> adjoint(const Array<T> &proj, const PolarGrid<T> &pg,
                                    const dims_t &recon_dims,
                                    const std::vector<std::array<T, 2>> &shifts) {

        auto c_cmplx = array::to_complex(proj);
        T scale = static_cast<T>(recon_dims.n1 * recon_dims.n2 * recon_dims.n3);

        c_cmplx = fft::fftshift2(c_cmplx);
        c_cmplx = fft::fft2(c_cmplx);
        c_cmplx = fft::ifftshift2(c_cmplx);

        if (!shifts.empty()) fft::phase_shift2d(c_cmplx, shifts);

        // discard aliased points (outside [-pi,pi]) before backprojecting
        std::transform(std::execution::par_unseq, c_cmplx.begin(), c_cmplx.end(),
                       pg.w.begin(), c_cmplx.begin(),
                       [](std::complex<T> c, T m) { return c * m; });

        std::array<Array<T>, 3> m_components;
        using complex_t = std::complex<T>;

        for (size_t i = 0; i < 3; ++i) {
            auto c_cmplx_copy = c_cmplx.clone();

            for (size_t j = 0; j < pg.nprojs(); ++j) {
                T coeff = beam_dir_vector(pg.angle(j), pg.gamma(j), pg.beta(j))[i];
                auto slice = c_cmplx_copy.slice(j, j + 1);
                std::for_each(std::execution::par_unseq, slice.begin(), slice.end(),
                              [coeff](complex_t &val) { val *= coeff; });
            }

            Array<complex_t> m_cmplx(recon_dims);
            nufft::nufft3d1<T>(c_cmplx_copy, m_cmplx, pg);
            m_components[i] = array::to_real<T>(m_cmplx) / scale;
        }

        return m_components;
    }
    template std::array<Array<float>, 3>
    adjoint<float>(const Array<float> &proj, const PolarGrid<float> &pg,
                   const dims_t &recon_dims,
                   const std::vector<std::array<float, 2>> &shifts);
    template std::array<Array<double>, 3>
    adjoint<double>(const Array<double> &proj, const PolarGrid<double> &pg,
                    const dims_t &recon_dims,
                    const std::vector<std::array<double, 2>> &shifts);

} // namespace tomocam
