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

// Check the GPU constant-magnitude solver (src/gpu/regularized_admm.cu) against the
// CPU solver (src/regularized_admm.cpp) and against the closed-form solution of the
// denoising problem (A = identity, no smoothing):
//   z = M0 y/|y| with M0 = mean|y| inside the constrained block, z = y outside.

#include <array>
#include <cmath>
#include <iostream>
#include <random>

#include "array.h"
#include "optimize.h"

#include "bregman.h"
#include "gpu/device_array.h"
#include "gpu/finitediff.h"
#include "gpu/gpu_opt.h"
#include "gpu/vec_array.h"

using namespace tomocam;
using CpuVec = opt::VecArray<double>;
using GpuVec = gpu::VecArray<double>;

static GpuVec to_device(const CpuVec &v) {
    return GpuVec{gpu::DeviceArray<double>(v[0]), gpu::DeviceArray<double>(v[1]),
                  gpu::DeviceArray<double>(v[2])};
}

static CpuVec to_host(const GpuVec &v) {
    return {v[0].to_host(), v[1].to_host(), v[2].to_host()};
}

static double max_abs_diff(const CpuVec &a, const CpuVec &b) {
    double d = 0.0;
    for (size_t c = 0; c < 3; ++c)
        for (size_t p = 0; p < a[c].size(); ++p)
            d = std::max(d, std::abs(a[c].begin()[p] - b[c].begin()[p]));
    return d;
}

static double rms_error(const CpuVec &a, const CpuVec &b) {
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
    dims_t dims = {8, 16, 16};
    size_t n = dims.size();
    const double M_true = 2.0;

    // smooth constant-magnitude field: slowly rotating direction
    CpuVec truth = opt::createVecArray<double>(dims);
    for (size_t k = 0; k < dims.n1; ++k)
        for (size_t i = 0; i < dims.n2; ++i)
            for (size_t j = 0; j < dims.n3; ++j) {
                double th =
                    0.15 * static_cast<double>(i) + 0.05 * static_cast<double>(k);
                double ph = 0.10 * static_cast<double>(j);
                truth[0][{k, i, j}] = M_true * std::sin(th) * std::cos(ph);
                truth[1][{k, i, j}] = M_true * std::sin(th) * std::sin(ph);
                truth[2][{k, i, j}] = M_true * std::cos(th);
            }

    std::mt19937 gen(7);
    std::normal_distribution<double> nd(0.0, 0.3);
    CpuVec y = opt::clone(truth);
    for (size_t c = 0; c < 3; ++c)
        for (auto &e : y[c]) e += nd(gen);

    opt::Function<double> A_cpu = [](const CpuVec &u) { return opt::clone(u); };
    gpu::opt::gpuFunction<double> A_gpu = [](const GpuVec &u) { return u.clone(); };

    // 1. projection helpers
    {
        CpuVec v = opt::clone(y);
        v[0].begin()[0] = v[1].begin()[0] = v[2].begin()[0] = 0.0;
        const double M0 = 2.5;
        auto v_dev = to_device(v);
        auto z = to_host(gpu::opt::project_magnitude(v_dev, M0));
        auto z_cpu = opt::project_magnitude(v, M0);
        ok &= check(max_abs_diff(z, z_cpu) < 1e-12, "projection matches CPU");
        ok &= check(z[0].begin()[0] == 0.0 && z[1].begin()[0] == 0.0 &&
                        z[2].begin()[0] == 0.0,
                    "projection leaves zero vectors at zero");
        ok &= check(std::abs(gpu::opt::mean_magnitude(v_dev) -
                             opt::mean_magnitude(v)) < 1e-12,
                    "mean_magnitude matches CPU");
    }

    // 2. no smoothing: CPU/GPU parity and closed-form solution, for the whole
    // volume and for a centred block
    for (dims_t roi : {dims, dims_t{4, 8, 8}}) {
        dims_t lo = (dims - roi) / 2;
        auto inside = [&](size_t k, size_t i, size_t j) {
            return k >= lo.n1 && k < lo.n1 + roi.n1 && i >= lo.n2 &&
                   i < lo.n2 + roi.n2 && j >= lo.n3 && j < lo.n3 + roi.n3;
        };
        auto mag = [&](const CpuVec &v, size_t k, size_t i, size_t j) {
            double a = v[0][{k, i, j}], b = v[1][{k, i, j}], c = v[2][{k, i, j}];
            return std::sqrt(a * a + b * b + c * c);
        };

        // closed form
        double M0 = 0.0;
        for (size_t k = 0; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j)
                    if (inside(k, i, j)) M0 += mag(y, k, i, j);
        M0 /= static_cast<double>(roi.size());
        CpuVec exact = opt::clone(y);
        for (size_t k = 0; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j)
                    if (inside(k, i, j)) {
                        double s = M0 / mag(y, k, i, j);
                        for (size_t c = 0; c < 3; ++c) exact[c][{k, i, j}] *= s;
                    }

        auto z_cpu = opt::constmag_l2<double>(A_cpu, y, y, 0.0, 1.0, 500, 50, 1e-12,
                                              1e-10, block_mask(dims, roi));
        auto y_dev = to_device(y);
        auto z_gpu = to_host(gpu::opt::constmag_l2<double>(
            A_gpu, y_dev, y_dev, 0.0, 1.0, 500, 50, 1e-12, 1e-10,
            gpu::DeviceArray<double>(block_mask(dims, roi))));

        double d_cpu = max_abs_diff(z_gpu, z_cpu);
        double d_exact = max_abs_diff(z_gpu, exact);
        std::cout << "  roi = " << roi.n1 << "x" << roi.n2 << "x" << roi.n3
                  << ": M0 = " << M0 << ", |gpu - cpu| = " << d_cpu
                  << ", |gpu - exact| = " << d_exact << "\n";
        ok &= check(d_cpu < 1e-8, "GPU matches CPU (lambda_l2 = 0)");
        ok &= check(d_exact < 1e-8, "GPU matches closed form (lambda_l2 = 0)");
    }

    // 3. with smoothing: uniform magnitude recovered from the data
    {
        auto y_dev = to_device(y);
        auto z_dev = gpu::opt::constmag_l2<double>(
            A_gpu, y_dev, y_dev, 0.05, 1.0, 200, 50, 1e-10, 1e-8,
            gpu::DeviceArray<double>::ones(dims));
        double M0 = gpu::opt::mean_magnitude(z_dev);
        auto z = to_host(z_dev);
        double max_dev = 0.0;
        for (size_t p = 0; p < n; ++p) {
            double a = z[0].begin()[p], b = z[1].begin()[p], c = z[2].begin()[p];
            max_dev =
                std::max(max_dev, std::abs(std::sqrt(a * a + b * b + c * c) - M0));
        }
        ok &= check(max_dev < 1e-8 * M0, "result has uniform |z| at every voxel");
        ok &= check(std::abs(M0 - M_true) < 0.05 * M_true,
                    "magnitude is recovered from the data");
        double err_noisy = rms_error(y, truth);
        double err_rec = rms_error(z, truth);
        std::cout << "  rms error noisy = " << err_noisy
                  << ", const-mag = " << err_rec << ", M0 = " << M0 << "\n";
        ok &= check(err_rec < err_noisy, "const-mag + L2 beats noisy data");
    }

    // 4a. finite-difference operators: CPU and GPU use the same forward
    // difference / Neumann boundary, so the Laplacian matches on the whole
    // volume, and GPU divergence is the exact adjoint of GPU grad
    {
        auto lap_c = opt::laplacian<double>(y[0]);
        auto lap_g =
            gpu::opt::neg_laplacian(gpu::DeviceArray<double>(y[0])).to_host();
        double dmax = 0;
        for (size_t k = 0; k < dims.n1; ++k)
            for (size_t i = 0; i < dims.n2; ++i)
                for (size_t j = 0; j < dims.n3; ++j)
                    dmax = std::max(dmax,
                                    std::abs(lap_c[{k, i, j}] - lap_g[{k, i, j}]));
        std::cout << "  laplacian max diff = " << dmax << "\n";
        ok &= check(dmax < 1e-12, "GPU neg_laplacian matches CPU laplacian");

        // <grad u, p> == <u, grad^T p> for random u, p
        std::mt19937 g(5);
        std::normal_distribution<double> nd(0.0, 1.0);
        auto rnd = [&]() {
            Array<double> a(dims);
            for (auto &v : a) v = nd(g);
            return a;
        };
        auto u = rnd();
        gpu::VecArray<double> p;
        double lhs = 0;
        std::array<Array<double>, 3> ph;
        for (size_t i = 0; i < 3; ++i) {
            ph[i] = rnd();
            p[i] = gpu::DeviceArray<double>(ph[i]);
        }
        auto gu = gpu::opt::grad(gpu::DeviceArray<double>(u));
        for (size_t i = 0; i < 3; ++i) lhs += array::dot(gu[i].to_host(), ph[i]);
        double rhs = array::dot(u, gpu::opt::neg_divergence(p).to_host());
        std::cout << "  adjoint <Du,p> = " << lhs << ", <u,D'p> = " << rhs << "\n";
        ok &= check(std::abs(lhs - rhs) < 1e-9 * std::abs(lhs),
                    "GPU neg_divergence is the adjoint of grad");
    }

    // 4. regularized_admm: every term combination runs and improves on the data
    {
        struct Case {
            const char *name;
            bool mag;
            double tv, l2;
        };
        for (Case c :
             {Case{"mag + L2", true, 0.0, 0.05}, Case{"TV only", false, 0.1, 0.0},
              Case{"L2 only", false, 0.0, 0.5}, Case{"mag + TV", true, 0.05, 0.0},
              Case{"mag + TV + L2", true, 0.05, 0.02}}) {
            gpu::opt::RegWeights<double> wg{c.mag, c.tv, c.l2};
            auto y_dev = to_device(y);
            auto z_dev = gpu::opt::regularized_admm<double>(
                A_gpu, y_dev, y_dev, wg, 1.0, 50, 20, 1e-10, 1e-8,
                gpu::DeviceArray<double>::ones(dims));
            auto z = to_host(z_dev);
            double err = rms_error(z, truth), err_noisy = rms_error(y, truth);
            std::cout << "  " << c.name << ": rms " << err << " (noisy " << err_noisy
                      << ")\n";
            ok &= check(err < err_noisy, c.name);
            if (c.mag) {
                double M0 = gpu::opt::mean_magnitude(z_dev), max_dev = 0.0;
                for (size_t p = 0; p < n; ++p) {
                    double a = z[0].begin()[p], b = z[1].begin()[p],
                           cc = z[2].begin()[p];
                    max_dev = std::max(
                        max_dev, std::abs(std::sqrt(a * a + b * b + cc * cc) - M0));
                }
                ok &= check(max_dev < 1e-8 * M0, "  uniform |z|");
            }
        }
        // no active term throws
        bool threw = false;
        try {
            auto y_dev = to_device(y);
            gpu::opt::regularized_admm<double>(
                A_gpu, y_dev, y_dev, gpu::opt::RegWeights<double>{}, 1.0, 1, 1,
                1e-10, 1e-8, gpu::DeviceArray<double>::ones(dims));
        } catch (const std::invalid_argument &) { threw = true; }
        ok &= check(threw, "no active term throws");
    }

    // 5. rho, continuation and zero_outside agree with the CPU solver
    {
        opt::Function<double> A_cpu = [](const CpuVec &u) { return opt::clone(u); };
        struct Case {
            const char *name;
            dims_t roi;
            bool zero_outside;
            double rho, growth, rho_max;
        };
        for (Case c :
             {Case{"separate rho", dims, false, 3.0, 1.0, 0.0},
              Case{"rho continuation", dims, false, 1.0, 1.1, 20.0},
              Case{"zero_outside + roi", {2, 8, 8}, true, 2.0, 1.05, 0.0}}) {
            opt::RegWeights<double> wc{true, 0.05, 0.02};
            gpu::opt::RegWeights<double> wg{true, 0.05, 0.02};
            opt::MagOptions<double> mc;
            gpu::opt::MagOptions<double> mg;
            mc.rho = mg.rho = c.rho;
            mc.rho_growth = mg.rho_growth = c.growth;
            mc.rho_max = mg.rho_max = c.rho_max;
            mc.zero_outside = mg.zero_outside = c.zero_outside;
            opt::AdmmHistory<double> hc;
            gpu::opt::AdmmHistory<double> hg;
            auto zc = opt::regularized_admm<double>(
                A_cpu, y, y, wc, 1.0, 30, 20, 1e-12, 0.0, block_mask(dims, c.roi),
                nullptr, mc, &hc);
            auto y_dev = to_device(y);
            auto zg = to_host(gpu::opt::regularized_admm<double>(
                A_gpu, y_dev, y_dev, wg, 1.0, 30, 20, 1e-12, 0.0,
                gpu::DeviceArray<double>(block_mask(dims, c.roi)), nullptr, mg,
                &hg));
            std::cout << "  " << c.name
                      << ": max |cpu - gpu| = " << max_abs_diff(zc, zg) << ", rho "
                      << hc.rho << " / " << hg.rho << "\n";
            ok &= check(max_abs_diff(zc, zg) < 1e-6 && hc.rho == hg.rho &&
                            std::abs(hc.Ms - hg.Ms) < 1e-6,
                        c.name);
        }

        // project_sphere keeps the previous value at a zero vector
        CpuVec v = opt::clone(y), prev = opt::clone(y);
        v[0].begin()[0] = v[1].begin()[0] = v[2].begin()[0] = 0.0;
        auto pc = opt::project_sphere<double>(v, 2.0, prev, 1e-12);
        auto pg = to_host(gpu::opt::project_sphere<double>(to_device(v), 2.0,
                                                           to_device(prev), 1e-12));
        ok &= check(max_abs_diff(pc, pg) < 1e-12, "project_sphere matches CPU");
    }

    std::cout << (ok ? "ALL PASSED\n" : "SOME FAILED\n");
    return ok ? 0 : 1;
}
