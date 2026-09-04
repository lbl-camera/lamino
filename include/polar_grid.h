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
#ifndef POLAR_GRID_H
#define POLAR_GRID_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <iostream>
#include <ranges>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "rotation.h"

namespace tomocam {

    // A non-owning view into one row of a ragged (variable-length-per-row)
    // flat buffer. Mirrors Slice<T>'s role for Array<T>.
    template <typename T>
    class RaggedSlice {
        T *data_;
        size_t n_;

      public:
        RaggedSlice(T *data, size_t n) : data_(data), n_(n) {}

        [[nodiscard]] size_t size() const { return n_; }
        T *begin() { return data_; }
        T *end() { return data_ + n_; }
        const T *begin() const { return data_; }
        const T *end() const { return data_ + n_; }
    };

    // Describes how a flat, per-nonuniform-point buffer is partitioned into
    // rows (one row per projection). Every row currently has the same
    // length (nrows*ncols) — out-of-[-pi,pi] points are masked via
    // PolarGrid::w rather than dropped — but this is kept generic (rather
    // than a plain dims_t) so any co-sized flat buffer (PolarGrid::x/y/z/w,
    // or an unrelated buffer like ToeplitzVectorOp's per-angle weights) can
    // be sliced against it uniformly.
    struct RaggedShape {
        std::vector<size_t> offsets; // size nrows()+1; row i is [offsets[i], offsets[i+1])

        [[nodiscard]] size_t nrows() const {
            return offsets.empty() ? 0 : offsets.size() - 1;
        }
        [[nodiscard]] size_t size() const {
            return offsets.empty() ? 0 : offsets.back();
        }
        [[nodiscard]] size_t row_begin(size_t i) const { return offsets[i]; }
        [[nodiscard]] size_t row_size(size_t i) const {
            return offsets[i + 1] - offsets[i];
        }

        // view into row i of any flat buffer sharing this row structure
        template <typename U>
        [[nodiscard]] RaggedSlice<U> slice(U *data, size_t i) const {
            return RaggedSlice<U>(data + offsets[i], row_size(i));
        }
        template <typename U>
        [[nodiscard]] RaggedSlice<const U> slice(const U *data, size_t i) const {
            return RaggedSlice<const U>(data + offsets[i], row_size(i));
        }
    };

    template <typename T>
    struct PolarGrid {
        size_t npts;
        std::vector<T> theta;
        std::vector<T> gammas;
        std::vector<T> betas;
        Array<T> x;
        Array<T> y;
        Array<T> z;
        Array<T> w; // 1 inside [-π,π]^3, 0 outside (masks aliased q-points)
        RaggedShape rows; // per-projection row structure of x/y/z/w

        // default constructor
        PolarGrid() : npts(0) {}

        // constructor — single gamma and beta for all angles
        PolarGrid(const std::vector<T> &angles, size_t nrows, size_t ncols,
                  T gamma, T beta = T(0)) {

            theta = angles;
            gammas = std::vector<T>(angles.size(), gamma);
            betas  = std::vector<T>(angles.size(), beta);
            dims_t dims = dims_t{theta.size(), nrows, ncols};
            npts = dims.n1 * dims.n2 * dims.n3;
            x = Array<T>(dims);
            y = Array<T>(dims);
            z = Array<T>(dims);
            w = Array<T>(dims);

            rows.offsets.resize(dims.n1 + 1);
            for (size_t i = 0; i <= dims.n1; ++i) {
                rows.offsets[i] = i * dims.n2 * dims.n3;
            }

            T L = 2 * M_PI;
            T dX = L / static_cast<T>(ncols);
            T dY = L / static_cast<T>(nrows);
            T L_half = L / 2;

#pragma omp parallel for collapse(3)
            for (size_t i = 0; i < dims.n1; ++i) {
                for (size_t j = 0; j < dims.n2; ++j) {
                    for (size_t k = 0; k < dims.n3; ++k) {

                        /*
                         * Offset 0.5 gives qX[N-1-k] = -qX[k] for ODD N, ensuring
                         * Hermitian-symmetric NUFFT output for real inputs.
                         * Even N is not supported: FINUFFT's asymmetric mode range
                         * [-N/2, N/2-1] breaks Hermitian symmetry at the Nyquist.
                         */
                        T qX = (k + 0.5) * dX - L_half;
                        T qY = (j + 0.5) * dY - L_half;

                        auto R = RotationTranspose(theta[i], gamma, beta);
                        auto q = matvec(R, {qX, qY, T(0)});
                        x[{i, j, k}] = q[0];
                        y[{i, j, k}] = q[1];
                        z[{i, j, k}] = q[2];
                        w[{i, j, k}] = (std::abs(q[0]) <= L_half &&
                                        std::abs(q[1]) <= L_half &&
                                        std::abs(q[2]) <= L_half) ? T(1) : T(0);
                    }
                }
            }
        }

        // constructor — one (gamma, beta) per dataset, concatenates all non-uniform points
        PolarGrid(
            const std::vector<std::tuple<std::vector<T>, T, T>> &angle_gamma_beta,
            size_t nrows, size_t ncols) {

            for (auto &[angles, gamma, beta] : angle_gamma_beta) {
                theta.insert(theta.end(), angles.begin(), angles.end());
                gammas.insert(gammas.end(), angles.size(), gamma);
                betas.insert(betas.end(), angles.size(), beta);
            }

            dims_t dims = dims_t{theta.size(), nrows, ncols};
            npts = dims.n1 * dims.n2 * dims.n3;
            x = Array<T>(dims);
            y = Array<T>(dims);
            z = Array<T>(dims);
            w = Array<T>(dims);

            rows.offsets.resize(dims.n1 + 1);
            for (size_t i = 0; i <= dims.n1; ++i) {
                rows.offsets[i] = i * dims.n2 * dims.n3;
            }

            T L = 2 * M_PI;
            T dX = L / static_cast<T>(ncols);
            T dY = L / static_cast<T>(nrows);
            T L_half = L / 2;

#pragma omp parallel for collapse(3)
            for (size_t i = 0; i < dims.n1; ++i) {
                for (size_t j = 0; j < dims.n2; ++j) {
                    for (size_t k = 0; k < dims.n3; ++k) {
                        T qX = (k + 0.5) * dX - L_half;
                        T qY = (j + 0.5) * dY - L_half;
                        auto R = RotationTranspose(theta[i], gammas[i], betas[i]);
                        auto q = matvec(R, {qX, qY, T(0)});
                        x[{i, j, k}] = q[0];
                        y[{i, j, k}] = q[1];
                        z[{i, j, k}] = q[2];
                        w[{i, j, k}] = (std::abs(q[0]) <= L_half &&
                                        std::abs(q[1]) <= L_half &&
                                        std::abs(q[2]) <= L_half) ? T(1) : T(0);
                    }
                }
            }
        }

        // delete copy constructor and assignment
        PolarGrid(const PolarGrid<T> &) = delete;
        PolarGrid<T> &operator=(const PolarGrid<T> &) = delete;

        // move constructor and assignment
        PolarGrid(PolarGrid<T> &&other) noexcept
            : npts(other.npts), theta(std::move(other.theta)),
              gammas(std::move(other.gammas)), betas(std::move(other.betas)),
              x(std::move(other.x)), y(std::move(other.y)), z(std::move(other.z)),
              w(std::move(other.w)), rows(std::move(other.rows)) {}

        PolarGrid<T> &operator=(PolarGrid<T> &&other) noexcept {
            if (this != &other) {
                npts = other.npts;
                theta = std::move(other.theta);
                gammas = std::move(other.gammas);
                betas = std::move(other.betas);
                x = std::move(other.x);
                y = std::move(other.y);
                z = std::move(other.z);
                w = std::move(other.w);
                rows = std::move(other.rows);
            }
            return *this;
        }

        PolarGrid<T> clone() const {
            PolarGrid<T> out;
            out.npts = this->npts;
            out.theta = this->theta;
            out.gammas = this->gammas;
            out.betas = this->betas;
            out.x = this->x.clone();
            out.y = this->y.clone();
            out.z = this->z.clone();
            out.w = this->w.clone();
            out.rows = this->rows;
            return out;
        }
        void print_limits() const {
            std::cout << "PolarGrid limits:\n";
            std::cout << std::format("qx: [{}, {}]\n", array::min(x), array::max(x));
            std::cout << std::format("qy: [{}, {}]\n", array::min(y), array::max(y));
            std::cout << std::format("qz: [{}, {}]\n", array::min(z), array::max(z));
        }

        // array dimensions for non-uniform points
        [[nodiscard]] dims_t dims() const { return x.dims(); }

        // size of the array
        [[nodiscard]] size_t size() const { return x.size(); }

        // theta values
        [[nodiscard]] const std::vector<T> &angles() const { return theta; }

        // get theta value for a given index
        [[nodiscard]] T angle(size_t i) const { return theta[i]; }

        // get gamma value for a given projection index
        [[nodiscard]] T gamma(size_t i) const { return gammas[i]; }

        // get beta value for a given projection index
        [[nodiscard]] T beta(size_t i) const { return betas[i]; }

        // number of angles
        [[nodiscard]] size_t nprojs() const { return theta.size(); }
    };

} // namespace tomocam

#endif // POLAR_GRID_H
