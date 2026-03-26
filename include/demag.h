
#ifndef DEMAG_H
#define DEMAG_H

#include <array>
#include <omp.h>

#include "array.h"
#include "fft.h"
#include "fftutils.h"

namespace tomocam::opt {
    template <typename T>
    std::array<Array<T>, 3> demag(const std::array<Array<T>, 3> &m) {

        auto dims = m[0].dims();
        auto kx = fft::fftfreq<T>(dims.n1);
        auto ky = fft::fftfreq<T>(dims.n2);
        auto kz = fft::rfftfreq<T>(dims.n3);

        std::array<Array<std::complex<T>>, 3> m_hat;
        for (int i = 0; i < 3; ++i) { m_hat[i] = fft::fft3_r2c(m[i]); }
        T norm_factor = 1.0 / static_cast<T>(dims.n1 * dims.n2 * dims.n3);

        // Pre-compute k² values to avoid redundant multiplications
        std::vector<T> kx_sq(dims.n1), ky_sq(dims.n2), kz_sq(dims.n3 / 2 + 1);
        for (size_t i = 0; i < dims.n1; ++i) { kx_sq[i] = kx[i] * kx[i]; }
        for (size_t j = 0; j < dims.n2; ++j) { ky_sq[j] = ky[j] * ky[j]; }
        for (size_t k = 0; k < dims.n3 / 2 + 1; ++k) { kz_sq[k] = kz[k] * kz[k]; }

        // Calculate the demagnetization field in Fourier space
        std::array<Array<std::complex<T>>, 3> demag_hat;
        dims_t fourier_dims = {dims.n1, dims.n2, dims.n3 / 2 + 1};
        for (int i = 0; i < 3; ++i) {
            demag_hat[i] = Array<std::complex<T>>(fourier_dims);
        }

        size_t nk = dims.n3 / 2 + 1;
// Parallelize the outer two loops using OpenMP
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < dims.n1; ++i) {
            for (size_t j = 0; j < dims.n2; ++j) {
                for (size_t k = 0; k < nk; ++k) {
                    size_t flat_idx = (i * dims.n2 + j) * nk + k;

                    // Compute k magnitude squared using pre-computed values
                    T k_mag_sq = kx_sq[i] + ky_sq[j] + kz_sq[k];

                    // Avoid division by zero at k=0
                    if (i == 0 && j == 0 && k == 0) { k_mag_sq = 1.0; }

                    // Pre-compute reciprocal to replace 3 divisions with 1 division
                    // + 3 multiplications
                    T k_inv = 1.0 / k_mag_sq;

                    // Cache array accesses using flat indexing
                    const auto &m0 = m_hat[0][flat_idx];
                    const auto &m1 = m_hat[1][flat_idx];
                    const auto &m2 = m_hat[2][flat_idx];

                    std::complex<T> k_dot_m_hat =
                        kx[i] * m0 + ky[j] * m1 + kz[k] * m2;

                    demag_hat[0][flat_idx] = -k_dot_m_hat * kx[i] * k_inv;
                    demag_hat[1][flat_idx] = -k_dot_m_hat * ky[j] * k_inv;
                    demag_hat[2][flat_idx] = -k_dot_m_hat * kz[k] * k_inv;
                }
            }
        }

        // Set DC component to zero
        for (size_t comp = 0; comp < 3; ++comp) { demag_hat[comp][{0, 0, 0}] = 0.0; }

        std::array<Array<T>, 3> H;
        for (size_t comp = 0; comp < 3; ++comp) {
            H[comp] = fft::fft3_c2r(demag_hat[comp], dims) * norm_factor;
        }
        return H;
    }
} // namespace tomocam::opt
#endif
