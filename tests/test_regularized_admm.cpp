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

#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "array.h"
#include "bregman.h"
#include "mask.h"
#include "optimize.h"
#include "tomocam.h"

using namespace tomocam;
using VecArray = opt::VecArray<double>;

static double rms_error(const VecArray &a, const VecArray &b) {
    auto d = opt::operator-(a, b);
    return std::sqrt(opt::dot(d, d) / static_cast<double>(a[0].size()));
}

static bool check(bool ok, const char *msg) {
    std::cout << (ok ? "PASS: " : "FAIL: ") << msg << "\n";
    return ok;
}

// 1 inside the block of size roi centred in dims (as crop3d places it), else 0
static Array<double> block_mask(const dims_t &dims, const dims_t &roi) {
    Array<double> m = Array<double>::zeros(dims);
    dims_t lo = (dims - roi) / 2;
    for (size_t k = 0; k < roi.n1; ++k)
        for (size_t i = 0; i < roi.n2; ++i)
            for (size_t j = 0; j < roi.n3; ++j)
                m[{k + lo.n1, i + lo.n2, j + lo.n3}] = 1.0;
    return m;
}

int main() {
    bool ok = true;
    dims_t dims = {4, 16, 16};
    size_t n = dims.size();
    Array<double> all_mask = Array<double>::ones(dims);

    // smooth unit-magnitude field: slowly rotating direction
    VecArray truth = opt::createVecArray<double>(dims);
    for (size_t k = 0; k < dims.n1; ++k)
        for (size_t i = 0; i < dims.n2; ++i)
            for (size_t j = 0; j < dims.n3; ++j) {
                double th =
                    0.15 * static_cast<double>(i) + 0.05 * static_cast<double>(k);
                double ph = 0.10 * static_cast<double>(j);
                truth[0][{k, i, j}] = std::sin(th) * std::cos(ph);
                truth[1][{k, i, j}] = std::sin(th) * std::sin(ph);
                truth[2][{k, i, j}] = std::cos(th);
            }

    // 1. projection helper
    {
        std::mt19937 gen(1);
        std::normal_distribution<double> nd(0.0, 1.0);
        VecArray v = opt::createVecArray<double>(dims);
        for (size_t c = 0; c < 3; ++c)
            for (auto &e : v[c]) e = nd(gen);
        // one zero vector
        v[0][{0, 0, 0}] = v[1][{0, 0, 0}] = v[2][{0, 0, 0}] = 0.0;
        const double M0 = 2.5;
        auto z = opt::project_magnitude(v, M0);
        double max_dev = 0.0, max_dir = 0.0;
        for (size_t p = 1; p < n; ++p) {
            double vx = v[0].begin()[p], vy = v[1].begin()[p], vz = v[2].begin()[p];
            double zx = z[0].begin()[p], zy = z[1].begin()[p], zz = z[2].begin()[p];
            double zn = std::sqrt(zx * zx + zy * zy + zz * zz);
            double vn = std::sqrt(vx * vx + vy * vy + vz * vz);
            max_dev = std::max(max_dev, std::abs(zn - M0));
            max_dir = std::max(max_dir, std::abs(zx * vn - vx * zn) / (zn * vn));
        }
        ok &= check(max_dev < 1e-12, "projection gives |z| = M0");
        ok &= check(max_dir < 1e-12, "projection preserves direction");
        ok &= check(z[0].begin()[0] == 0.0 && z[1].begin()[0] == 0.0 &&
                        z[2].begin()[0] == 0.0,
                    "projection leaves zero vectors at zero");
        ok &= check(
            std::abs(opt::mean_magnitude(z) - M0 * static_cast<double>(n - 1) /
                                                  static_cast<double>(n)) < 1e-12,
            "mean_magnitude of projected field = M0 (one zero voxel)");
    }

    // 2. constant magnitude + L2 on a denoising problem (A = identity); the
    // magnitude is not given, it has to be recovered from the data
    {
        const double M_true = 2.0;
        std::mt19937 gen(7);
        std::normal_distribution<double> nd(0.0, 0.3);
        VecArray clean = opt::clone(truth);
        for (size_t c = 0; c < 3; ++c) clean[c] *= M_true;
        VecArray y = opt::clone(clean);
        for (size_t c = 0; c < 3; ++c)
            for (auto &e : y[c]) e += nd(gen);

        opt::Function<double> A = [](const VecArray &u) { return opt::clone(u); };
        auto x0 = opt::clone(y);

        auto z = opt::constmag_l2<double>(A, y, x0, 0.05, 1.0, 200, 50, 1e-10, 1e-8,
                                          all_mask);

        // uniform magnitude at every voxel
        double M0 = opt::mean_magnitude(z);
        double max_dev = 0.0;
        for (size_t p = 0; p < n; ++p) {
            double a = z[0].begin()[p], b = z[1].begin()[p], c = z[2].begin()[p];
            max_dev =
                std::max(max_dev, std::abs(std::sqrt(a * a + b * b + c * c) - M0));
        }
        ok &= check(max_dev < 1e-8 * M0, "result has uniform |z| at every voxel");
        ok &= check(std::abs(M0 - M_true) < 0.05 * M_true,
                    "magnitude is recovered from the data");

        double err_noisy = rms_error(y, clean);
        double err_rec = rms_error(z, clean);
        std::cout << "  rms error noisy = " << err_noisy
                  << ", const-mag = " << err_rec << ", M0 = " << M0 << "\n";
        ok &= check(err_rec < err_noisy, "const-mag + L2 beats noisy data");
    }

    // 3. constraint restricted to a centred block: the padding stays free
    {
        dims_t roi = {2, 8, 8};
        dims_t lo = {1, 4, 4};
        auto inside = [&](size_t k, size_t i, size_t j) {
            return k >= lo.n1 && k < lo.n1 + roi.n1 && i >= lo.n2 &&
                   i < lo.n2 + roi.n2 && j >= lo.n3 && j < lo.n3 + roi.n3;
        };

        // unit magnitude inside the block, 5x larger in the padding
        VecArray y = opt::clone(truth);
        for (size_t k = 0; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j)
                    if (!inside(k, i, j))
                        for (size_t c = 0; c < 3; ++c) y[c][{k, i, j}] *= 5.0;

        opt::Function<double> A = [](const VecArray &u) { return opt::clone(u); };
        auto x0 = opt::clone(y);

        // no smoothing: the solution is y itself in the padding
        auto z = opt::constmag_l2<double>(A, y, x0, 0.0, 1.0, 500, 50, 1e-12, 1e-10,
                                          block_mask(dims, roi));

        double sum_in = 0.0, max_dev_in = 0.0, max_err_out = 0.0;
        std::vector<double> mags;
        for (size_t k = 0; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j) {
                    double a = z[0][{k, i, j}], b = z[1][{k, i, j}],
                           c = z[2][{k, i, j}];
                    if (inside(k, i, j)) {
                        mags.push_back(std::sqrt(a * a + b * b + c * c));
                        sum_in += mags.back();
                    } else {
                        for (size_t d = 0; d < 3; ++d)
                            max_err_out =
                                std::max(max_err_out, std::abs(z[d][{k, i, j}] -
                                                               y[d][{k, i, j}]));
                    }
                }
        double M0 = sum_in / static_cast<double>(mags.size());
        for (double m : mags) max_dev_in = std::max(max_dev_in, std::abs(m - M0));
        std::cout << "  M0 = " << M0 << ", max dev inside = " << max_dev_in
                  << ", max error outside = " << max_err_out << "\n";
        ok &= check(mags.size() == roi.size(), "block has the expected size");
        ok &= check(max_dev_in < 1e-8, "uniform |z| inside the block");
        ok &=
            check(std::abs(M0 - 1.0) < 1e-6, "M0 is estimated from the block only");
        ok &= check(max_err_out < 1e-6, "padding is left unconstrained");
    }

    // 3b. non-rectangular support (a disc hole) and support_mask()
    {
        // 2D file mask: 1 everywhere except a disc of radius 3 at the centre
        Array<double> file_mask(1, dims.n2, dims.n3);
        for (size_t i = 0; i < dims.n2; ++i)
            for (size_t j = 0; j < dims.n3; ++j) {
                double di = double(i) - 8.0, dj = double(j) - 8.0;
                file_mask[{0, i, j}] = (di * di + dj * dj > 9.0) ? 1.0 : 0.0;
            }
        auto mask = support_mask(file_mask, dims);
        bool extruded = mask.dims() == dims;
        for (size_t k = 1; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j)
                    extruded &= mask[{k, i, j}] == mask[{0, i, j}];
        ok &= check(extruded, "support_mask extrudes a 2D mask over n1");

        bool threw = false;
        try {
            support_mask(Array<double>(1, dims.n2 - 1, dims.n3), dims);
        } catch (const std::invalid_argument &) { threw = true; }
        ok &= check(threw, "support_mask rejects a mismatching mask");

        // 3D mask window: centred in all three axes, free (0) outside it
        {
            const std::array<size_t, 3> md = {dims.n1 - 2, dims.n2 - 4, dims.n3 - 4};
            Array<double> m3d(md[0], md[1], md[2]);
            for (size_t k = 0; k < md[0]; ++k)
                for (size_t i = 0; i < md[1]; ++i)
                    for (size_t j = 0; j < md[2]; ++j)
                        m3d[{k, i, j}] = (k == 0) ? 0.0 : 1.0;
            auto w = support_mask(m3d, dims, md);
            bool win = w.dims() == dims;
            for (size_t k = 0; k < dims.n1; ++k)
                for (size_t i = 0; i < dims.n2; ++i)
                    for (size_t j = 0; j < dims.n3; ++j) {
                        bool in = k >= 1 && k < 1 + md[0] && i >= 2 &&
                                  i < 2 + md[1] && j >= 2 && j < 2 + md[2];
                        double expect = (in && k != 1) ? 1.0 : 0.0;
                        win &= w[{k, i, j}] == expect;
                    }
            ok &= check(win, "support_mask centres a 3D mask in all axes");

            auto ones = support_mask(Array<double>(), dims, md);
            double sum = 0;
            for (size_t i = 0; i < ones.size(); ++i) sum += ones.begin()[i];
            ok &= check(sum == double(md[0] * md[1] * md[2]),
                        "mask_dims without a file gives a centred window of ones");
        }

        // 5x larger magnitude inside the hole; no smoothing, so z = y there
        VecArray y = opt::clone(truth);
        for (size_t k = 0; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j)
                    if (mask[{k, i, j}] < 0.5)
                        for (size_t c = 0; c < 3; ++c) y[c][{k, i, j}] *= 5.0;
        opt::Function<double> A = [](const VecArray &u) { return opt::clone(u); };
        auto z =
            opt::constmag_l2<double>(A, y, y, 0.0, 1.0, 500, 50, 1e-12, 1e-10, mask);
        double sum = 0.0, cnt = 0.0, err_out = 0.0;
        for (size_t p = 0; p < n; ++p) {
            double a = z[0].begin()[p], b = z[1].begin()[p], c = z[2].begin()[p];
            if (mask.begin()[p] > 0.5) {
                sum += std::sqrt(a * a + b * b + c * c);
                cnt += 1.0;
            } else {
                for (size_t d = 0; d < 3; ++d)
                    err_out = std::max(err_out,
                                       std::abs(z[d].begin()[p] - y[d].begin()[p]));
            }
        }
        double M0 = sum / cnt, max_dev = 0.0;
        for (size_t p = 0; p < n; ++p) {
            if (mask.begin()[p] < 0.5) continue;
            double a = z[0].begin()[p], b = z[1].begin()[p], c = z[2].begin()[p];
            max_dev =
                std::max(max_dev, std::abs(std::sqrt(a * a + b * b + c * c) - M0));
        }
        ok &= check(max_dev < 1e-8, "uniform |z| inside a non-rectangular mask");
        ok &= check(std::abs(M0 - 1.0) < 1e-6, "M0 estimated from the mask only");
        ok &= check(err_out < 1e-6, "voxels outside the mask are free");
    }

    // 3c. zero_outside: x and z are exactly 0 outside the mask (masked CG)
    {
        Array<double> file_mask(1, dims.n2, dims.n3);
        for (size_t i = 0; i < dims.n2; ++i)
            for (size_t j = 0; j < dims.n3; ++j)
                file_mask[{0, i, j}] =
                    (i >= 3 && i < 13 && j >= 3 && j < 13) ? 1.0 : 0.0;
        auto mask = support_mask(file_mask, dims);
        opt::Function<double> A = [](const VecArray &u) { return opt::clone(u); };
        auto x0 = opt::clone(truth);

        opt::RegWeights<double> w{true, 0.0, 0.05};
        opt::MagOptions<double> mo;
        mo.zero_outside = true;
        auto z = opt::regularized_admm<double>(A, truth, x0, w, 1.0, 20, 20, 1e-12,
                                               0.0, mask, nullptr, mo);
        auto zl = opt::constmag_l2<double>(A, truth, x0, 0.05, 1.0, 20, 20, 1e-12,
                                           0.0, mask, nullptr, true);
        double out = 0.0, in = 0.0;
        for (size_t c = 0; c < 3; ++c)
            for (size_t p = 0; p < n; ++p) {
                if (mask.begin()[p] < 0.5)
                    out = std::max(out, std::max(std::abs(z[c].begin()[p]),
                                                 std::abs(zl[c].begin()[p])));
                else
                    in = std::max(in, std::abs(z[c].begin()[p]));
            }
        // x0 is nonzero outside the mask: the multiplier b must still be 0 there
        opt::AdmmHistory<double> hist;
        opt::regularized_admm<double>(A, truth, x0, w, 1.0, 5, 20, 1e-12, 0.0, mask,
                                      nullptr, mo, &hist);
        double b_out = 0.0;
        for (size_t c = 0; c < 3; ++c)
            for (size_t p = 0; p < n; ++p)
                if (mask.begin()[p] < 0.5)
                    b_out = std::max(b_out, std::abs(hist.b_z[c].begin()[p]));
        ok &= check(b_out == 0.0, "zero_outside: x0 is projected, b = 0 outside");
        ok &= check(out == 0.0, "zero_outside: z = 0 outside the mask");
        ok &= check(in > 0.0, "zero_outside: z nonzero inside the mask");

        // masked CG alone: solution of P A P x = P y is P y for A = I
        auto x = opt::cgsolver<double>(A, truth, x0, 50, 1e-14, 1e-14, mask);
        double dev = 0.0;
        for (size_t c = 0; c < 3; ++c)
            for (size_t p = 0; p < n; ++p)
                dev = std::max(dev, std::abs(x[c].begin()[p] -
                                             truth[c].begin()[p] * mask.begin()[p]));
        ok &= check(dev < 1e-10, "masked CG solves the restricted system");
    }

    // 4. L2 operator is SPD: <u, (lambda ∇^T∇) u> >= 0
    {
        double q = 0.0;
        for (size_t c = 0; c < 3; ++c) {
            auto lap = opt::laplacian<double>(truth[c]);
            q += array::dot(truth[c], lap);
        }
        ok &= check(q >= 0.0, "<u, ∇^T∇ u> >= 0");
    }

    // 4b. divergence is the exact adjoint of grad_u: <grad u, p> == <u, D^T p>
    {
        std::mt19937 g(5);
        std::normal_distribution<double> nd(0.0, 1.0);
        auto rnd = [&]() {
            Array<double> a(truth[0].dims());
            for (auto &v : a) v = nd(g);
            return a;
        };
        auto u = rnd();
        std::array<Array<double>, 3> p{rnd(), rnd(), rnd()};
        auto gu = opt::grad_u<double>(u);
        double lhs = 0;
        for (size_t i = 0; i < 3; ++i) lhs += array::dot(gu[i], p[i]);
        double rhs = array::dot(u, opt::divergence<double>(p));
        ok &= check(std::abs(lhs - rhs) < 1e-9 * std::abs(lhs),
                    "divergence is the adjoint of grad_u");
    }

    // 5. regularized_admm: each term can run independently or together
    {
        const double M_true = 2.0;
        std::mt19937 gen(11);
        std::normal_distribution<double> nd(0.0, 0.3);
        VecArray clean = opt::clone(truth);
        for (size_t c = 0; c < 3; ++c) clean[c] *= M_true;
        VecArray y = opt::clone(clean);
        for (size_t c = 0; c < 3; ++c)
            for (auto &e : y[c]) e += nd(gen);
        opt::Function<double> A = [](const VecArray &u) { return opt::clone(u); };
        auto x0 = opt::clone(y);

        // mag + L2 reproduces constmag_l2
        auto ref = opt::constmag_l2<double>(A, y, x0, 0.05, 1.0, 20, 20, 1e-10, 0.0,
                                            all_mask);
        opt::RegWeights<double> w_ml{true, 0.0, 0.05};
        auto z_ml = opt::regularized_admm<double>(A, y, x0, w_ml, 1.0, 20, 20, 1e-10,
                                                  0.0, all_mask);
        ok &= check(rms_error(ref, z_ml) < 1e-10, "mag + L2 matches constmag_l2");

        // TV only reproduces split_bregman
        auto sb = opt::split_bregman<double>(A, y, x0, 0.1, 1.0, 20, 20, 1e-10, 0.0,
                                             Array<double>());
        opt::RegWeights<double> w_tv{false, 0.1, 0.0};
        auto x_tv = opt::regularized_admm<double>(A, y, x0, w_tv, 1.0, 20, 20, 1e-10,
                                                  0.0, all_mask);
        ok &= check(rms_error(sb, x_tv) < 1e-8, "TV only matches split_bregman");

        // all three on: uniform magnitude, closer to truth than the data
        opt::RegWeights<double> w_all{true, 0.05, 0.02};
        auto z_all = opt::regularized_admm<double>(A, y, x0, w_all, 1.0, 100, 20,
                                                   1e-10, 1e-8, all_mask);
        double M0 = opt::mean_magnitude(z_all);
        double max_dev = 0.0;
        for (size_t p = 0; p < n; ++p) {
            double a = z_all[0].begin()[p], b = z_all[1].begin()[p],
                   c = z_all[2].begin()[p];
            max_dev =
                std::max(max_dev, std::abs(std::sqrt(a * a + b * b + c * c) - M0));
        }
        ok &= check(max_dev < 1e-8 * M0, "mag + TV + L2 has uniform |z|");
        ok &= check(rms_error(z_all, clean) < rms_error(y, clean),
                    "mag + TV + L2 beats noisy data");

        // L2 only is a plain Tikhonov solve (no constraint): runs and reduces error
        opt::RegWeights<double> w_l2{false, 0.0, 0.5};
        auto x_l2 = opt::regularized_admm<double>(A, y, x0, w_l2, 1.0, 5, 50, 1e-10,
                                                  1e-8, all_mask);
        ok &= check(rms_error(x_l2, clean) < rms_error(y, clean),
                    "L2 only beats noisy data");

        // nothing enabled is an error
        bool threw = false;
        try {
            opt::regularized_admm<double>(A, y, x0, opt::RegWeights<double>{}, 1.0,
                                          1, 1, 1e-10, 1e-8, all_mask);
        } catch (const std::invalid_argument &) { threw = true; }
        ok &= check(threw, "no active term throws");
    }

    // 4. project_sphere, rho continuation, history, warm start
    {
        const double M_true = 2.0;
        std::mt19937 gen(11);
        std::normal_distribution<double> nd(0.0, 0.3);
        VecArray clean = opt::clone(truth);
        for (size_t c = 0; c < 3; ++c) clean[c] *= M_true;
        VecArray y = opt::clone(clean);
        for (size_t c = 0; c < 3; ++c)
            for (auto &e : y[c]) e += nd(gen);
        opt::Function<double> A = [](const VecArray &u) { return opt::clone(u); };

        // project_sphere: idempotent, zero vectors keep the previous value
        VecArray prev = opt::createVecArray<double>(dims);
        for (size_t c = 0; c < 3; ++c)
            for (auto &e : prev[c]) e = (c == 1) ? M_true : 0.0;
        VecArray v = opt::clone(y);
        v[0][{0, 0, 0}] = v[1][{0, 0, 0}] = v[2][{0, 0, 0}] = 0.0;
        auto p1 = opt::project_sphere<double>(v, M_true, prev, 1e-12);
        auto p2 = opt::project_sphere<double>(p1, M_true, prev, 1e-12);
        ok &= check(rms_error(p1, p2) < 1e-14, "project_sphere is idempotent");
        ok &= check(p1[0].begin()[0] == 0.0 && p1[1].begin()[0] == M_true &&
                        p1[2].begin()[0] == 0.0,
                    "project_sphere keeps previous z at a zero vector");

        // continuation: rho only grows, |z| = Ms, magnitude residuals shrink
        opt::RegWeights<double> w{true, 0.0, 0.05};
        opt::MagOptions<double> mo;
        mo.rho = 1.0;
        mo.rho_growth = 1.05;
        mo.rho_max = 50.0;
        opt::AdmmHistory<double> hist;
        auto z = opt::regularized_admm<double>(A, y, y, w, 1.0, 300, 50, 1e-12, 0.0,
                                               all_mask, nullptr, mo, &hist);
        bool monotone = true;
        for (size_t i = 1; i < hist.iters.size(); ++i)
            monotone &= hist.iters[i].rho >= hist.iters[i - 1].rho;
        ok &= check(monotone && hist.rho <= 50.0 && hist.rho > 1.0,
                    "rho is non-decreasing and capped by rho_max");
        double Ms = hist.Ms, max_dev = 0.0;
        for (size_t p = 0; p < n; ++p) {
            double a = z[0].begin()[p], b = z[1].begin()[p], c = z[2].begin()[p];
            max_dev =
                std::max(max_dev, std::abs(std::sqrt(a * a + b * b + c * c) - Ms));
        }
        ok &= check(max_dev < 1e-6, "max | |z| - Ms | < 1e-6 with continuation");
        ok &= check(hist.iters.back().primal_mag < 1e-4,
                    "primal magnitude residual < 1e-4");
        ok &= check(hist.iters.back().dual_mag < 1e-4,
                    "dual magnitude residual < 1e-4");
        ok &= check(rms_error(z, clean) < rms_error(y, clean),
                    "continuation result beats noisy data");

        // rescaling keeps the unscaled multiplier rho * b_z fixed
        opt::AdmmHistory<double> h_grow, h_ref;
        opt::MagOptions<double> m_grow = mo, m_ref = mo;
        m_grow.rho_growth = 2.0;
        m_grow.rho_max = 0.0;
        m_ref.rho_growth = 1.0;
        opt::regularized_admm<double>(A, y, y, w, 1.0, 1, 50, 1e-12, 0.0, all_mask,
                                      nullptr, m_grow, &h_grow);
        opt::regularized_admm<double>(A, y, y, w, 1.0, 1, 50, 1e-12, 0.0, all_mask,
                                      nullptr, m_ref, &h_ref);
        VecArray lhs = opt::clone(h_grow.b_z), rhs = opt::clone(h_ref.b_z);
        for (size_t c = 0; c < 3; ++c) {
            lhs[c] *= h_grow.rho;
            rhs[c] *= h_ref.rho;
        }
        ok &= check(h_grow.rho == 2.0 * mo.rho && rms_error(lhs, rhs) < 1e-12,
                    "rho rescaling leaves rho * b_z unchanged");

        // zero_outside: z = 0 outside the roi
        dims_t roi = {2, 8, 8};
        opt::MagOptions<double> mz;
        mz.zero_outside = true;
        auto zz = opt::regularized_admm<double>(A, y, y, w, 1.0, 20, 50, 1e-12, 0.0,
                                                block_mask(dims, roi), nullptr, mz);
        bool zero_out = true;
        for (size_t k = 0; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j) {
                    bool in =
                        k >= 1 && k < 3 && i >= 4 && i < 12 && j >= 4 && j < 12;
                    if (!in)
                        for (size_t c = 0; c < 3; ++c)
                            zero_out &= zz[c][{k, i, j}] == 0.0;
                }
        ok &= check(zero_out, "zero_outside gives z = 0 outside the roi");

        // initialisation sensitivity (report); warm start must not be worse
        // than the best random start. Objective is evaluated at z with A = I.
        auto objective = [&](const VecArray &zv) {
            double data = 0.0, g2 = 0.0;
            auto r = opt::operator-(zv, y);
            data = 0.5 * opt::dot(r, r);
            for (size_t c = 0; c < 3; ++c) {
                auto g = opt::grad_u(zv[c]);
                for (auto &gi : g)
                    for (auto e = gi.begin(); e != gi.end(); ++e) g2 += *e * *e;
            }
            return data + 0.5 * w.lambda_l2 * g2;
        };
        opt::MagOptions<double> mc = mo;
        auto run = [&](const VecArray &start) {
            return opt::regularized_admm<double>(A, y, start, w, 1.0, 300, 50, 1e-12,
                                                 0.0, all_mask, nullptr, mc);
        };
        double best_random = 1e300;
        for (unsigned seed = 0; seed < 3; ++seed) {
            std::mt19937 g2(100 + seed);
            std::normal_distribution<double> u(0.0, 1.0);
            VecArray start = opt::createVecArray<double>(dims);
            for (size_t c = 0; c < 3; ++c)
                for (auto &e : start[c]) e = u(g2);
            auto zr = run(start);
            double f = objective(zr);
            best_random = std::min(best_random, f);
            std::cout << "  random start " << seed << ": objective " << f
                      << ", rms error " << rms_error(zr, clean) << "\n";
        }
        opt::RegWeights<double> w_unc{false, 0.0, w.lambda_l2};
        auto warm0 = opt::regularized_admm<double>(A, y, y, w_unc, 1.0, 20, 50,
                                                   1e-12, 0.0, all_mask);
        auto zw = run(warm0);
        double fw = objective(zw);
        std::cout << "  warm start: objective " << fw << ", rms error "
                  << rms_error(zw, clean) << "\n";
        ok &= check(fw <= best_random * (1.0 + 1e-6),
                    "warm start objective <= best random start");
    }

    std::cout << (ok ? "ALL PASSED\n" : "SOME FAILED\n");
    return ok ? 0 : 1;
}
