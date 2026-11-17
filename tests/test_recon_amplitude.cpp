// Amplitude of the reconstruction.
//
// Projection data follow the convention of the `forward` program,
// y = forward(m) / projection_scale(...) (see padding.h). The solvers multiply the
// data by the matching projection_scale, so that
//
//  1. the Toeplitz normal operator equals adjoint(forward(x)) including its
//     scale, also when the polar grid is larger than the volume (V < G);
//  2. the solver's model reproduces the scaled data exactly when x = m is placed
//     in its padded volume; with 1. this makes x = m the solution of the normal
//     equations, i.e. the reconstruction has the amplitude of m;
//  3. MBIR2 end to end: the reconstruction is linear in the data (there is no
//     hidden normalisation by the data range; it used to divide by max(y)), and
//     its best data-fit scale is ~1. That check is loose (5%): the solver pads
//     n1, part of the signal spreads into the padding along z (limited tilt),
//     and the output is cropped to recon_dims.

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <random>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "padding.h"
#include "polar_grid.h"
#include "projection.h"
#include "recon_params.h"
#include "toeplitz.h"
#include "tomocam.h"

using namespace tomocam;
using vec3 = std::array<Array<double>, 3>;

static double vec_dot(const vec3 &a, const vec3 &b) {
    double s = 0;
    for (size_t i = 0; i < 3; ++i) s += array::dot<double>(a[i], b[i]);
    return s;
}

// two smooth blobs with different directions, well inside an n1 x n x n box
static vec3 blob_field(size_t n1, size_t n) {
    vec3 m;
    for (auto &c : m) c = Array<double>::zeros(dims_t{n1, n, n});
    const double c1 = 0.5 * (n1 - 1), c2 = 0.5 * (n - 1);
    const std::array<std::array<double, 6>, 2> blobs = {{
        {0.0, -0.15, 0.1, 0.8, -0.5, 0.3}, // offset (z, y, x)/n, direction
        {0.0, 0.15, -0.12, -0.2, 0.7, -0.4},
    }};
    for (size_t i = 0; i < n1; ++i)
        for (size_t j = 0; j < n; ++j)
            for (size_t k = 0; k < n; ++k)
                for (const auto &b : blobs) {
                    double dz = (i - c1) / n1 - b[0], dy = (j - c2) / n - b[1],
                           dx = (k - c2) / n - b[2];
                    double g =
                        std::exp(-(dz * dz / 0.04 + (dy * dy + dx * dx) / 0.01));
                    for (size_t c = 0; c < 3; ++c) m[c][{i, j, k}] += b[3 + c] * g;
                }
    return m;
}

// projections as written by the forward program (src/forward.cpp)
static Array<double> forward_data(const vec3 &m, const std::vector<double> &theta,
                                  double gamma, size_t det) {
    dims_t vol = m[0].dims();
    size_t g =
        padded_dim(std::max(vol.n2, det), static_cast<double>(DEFAULT_PAD_FACTOR));
    PolarGrid<double> pg(theta, g, g, gamma, 0.0);
    auto p = forward<double>(m, pg, gamma, 0.0);
    p /= projection_scale<double>(vol.n1, vol, g, g);
    return crop2d<double>(p, dims_t{theta.size(), det, det}, PadType::SYMMETRIC);
}

static std::vector<double> tilt_series() {
    std::vector<double> th;
    for (int a = -70; a <= 70; a += 5) th.push_back(a * M_PI / 180.0);
    return th;
}

// 1. Toeplitz A^T A vs adjoint(forward), scale included
static bool test_normal_operator_scale(size_t n, size_t G) {
    std::mt19937 gen(5);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    dims_t vol{7, n, n};
    vec3 x;
    for (auto &c : x) {
        c = Array<double>::zeros(vol);
        for (size_t i = 0; i < c.size(); ++i) c[i] = dis(gen);
    }
    auto th = tilt_series();
    PolarGrid<double> pg(th, G, G, 0.3, 0.0);
    auto ata = adjoint<double>(forward<double>(x, pg, 0.3, 0.0), pg, vol, 0.3, 0.0);
    cpu::ToeplitzVectorOp<double> op(pg, vol);
    auto toe = cpu::sysmat(x, op);
    double s = vec_dot(ata, toe) / vec_dot(toe, toe);
    double err = 0;
    for (size_t i = 0; i < 3; ++i) {
        auto d = ata[i] - toe[i];
        err += array::dot<double>(d, d);
    }
    err = std::sqrt(err / vec_dot(ata, ata));
    bool ok = std::abs(s - 1) < 1e-6 && err < 1e-6;
    std::cout << std::format(
        "Toeplitz vs adjoint(forward), V={} G={}: scale {:.8f}, "
        "err {:.2e} {}\n",
        n, G, s, err, ok ? "PASS" : "FAIL");
    return ok;
}

static vec3 reconstruct(const vec3 &m, size_t det, std::array<size_t, 3> rdims,
                        double amp) {
    auto th = tilt_series();
    std::vector<Dataset_t<double>> ds;
    for (double gamma : {0.0, M_PI / 4}) {
        Dataset_t<double> d;
        d.projs = forward_data(m, th, gamma, det);
        d.projs *= amp;
        d.angles = th;
        d.gamma = gamma;
        ds.push_back(std::move(d));
    }
    ReconParams p;
    p.recon_dims = rdims;
    p.mask_dims = rdims;
    p.maxIters = 400;
    p.tol = 1e-12f;
    p.xtol = 1e-12f;
    p.logMode = LogMode::SILENT;
    dims_t rd{rdims[0], rdims[1], rdims[2]};
    return MBIR2<double>(ds, p, Array<double>::ones(rd));
}

// best data-fit scale c of the reconstruction: forward(c x) ~ y
static double data_fit_scale(const vec3 &x, const vec3 &m, size_t det) {
    auto th = tilt_series();
    double num = 0, den = 0;
    for (double gamma : {0.0, M_PI / 4}) {
        auto y = forward_data(m, th, gamma, det);
        auto p = forward_data(x, th, gamma, det);
        num += array::dot<double>(p, y);
        den += array::dot<double>(p, p);
    }
    return num / den;
}

// 2. data convention vs the solver model, as set up in MBIR / MBIR2
static bool test_data_convention(size_t n1, size_t n, size_t det) {
    auto m = blob_field(n1, n);
    auto th = tilt_series();
    const double padfac = DEFAULT_PAD_FACTOR;
    size_t G = recon_grid_dim<double>(det, n, padfac);
    dims_t out{padded_dim(n1, padfac), recon_vol_dim(det, n), recon_vol_dim(det, n)};
    const double scale = projection_scale<double>(n1, out, G, G);
    vec3 me;
    dims_t off{(out.n1 - n1) / 2, (out.n2 - n) / 2, (out.n3 - n) / 2};
    for (size_t i = 0; i < 3; ++i) me[i] = pad3d<double>(m[i], out, off);
    double err = 0, ref = 0;
    for (double gamma : {0.0, M_PI / 4}) {
        auto y = forward_data(m, th, gamma, det);
        y *= scale;
        PolarGrid<double> pg(th, G, G, gamma, 0.0);
        auto model = crop2d<double>(forward<double>(me, pg, gamma, 0.0),
                                    dims_t{th.size(), det, det}, PadType::SYMMETRIC);
        auto d = model - y;
        err += array::dot<double>(d, d);
        ref += array::dot<double>(y, y);
    }
    err = std::sqrt(err / ref);
    // when the solver grid differs from the one `forward` used (detector smaller
    // than the volume), the periodic ifft2 wraps the blob tails differently:
    // ~1e-6, not a scale error (a wrong scale shows up as O(1))
    bool ok = err < 1e-5;
    std::cout << std::format("solver model of m vs scaled data, {}x{}^2, detector "
                             "{} (volume {}x{}^2, grid {}): err {:.2e} {}\n",
                             n1, n, det, out.n1, out.n2, G, err,
                             ok ? "PASS" : "FAIL");
    return ok;
}

// 3. end to end with MBIR2 (CPU)
static bool test_end_to_end(size_t n1, size_t n, size_t det) {
    auto m = blob_field(n1, n);
    std::array<size_t, 3> rd = {n1, n, n};
    auto x1 = reconstruct(m, det, rd, 1.0);
    auto x3 = reconstruct(m, det, rd, 3.0);
    double c = data_fit_scale(x1, m, det);
    double lin = 0;
    for (size_t i = 0; i < 3; ++i) {
        auto d = x3[i] - x1[i] * 3.0;
        lin += array::dot<double>(d, d);
    }
    lin = std::sqrt(lin / (9 * vec_dot(x1, x1)));
    size_t G = padded_dim(std::max(n, det), static_cast<double>(DEFAULT_PAD_FACTOR));
    // parallel reductions make CG differ slightly between runs (~1e-4); a
    // normalisation by the data range would give an error of 2/3
    bool ok = std::abs(c - 1) < 0.05 && lin < 1e-3;
    std::cout << std::format("end to end {}x{}^2, detector {} (grid {}): data-fit "
                             "scale {:.4f}, linearity err {:.1e} {}\n",
                             n1, n, det, G, c, lin, ok ? "PASS" : "FAIL");
    return ok;
}

int main() {
    bool ok = true;
    ok = test_normal_operator_scale(21, 29) && ok; // V < G (default layout)
    ok = test_normal_operator_scale(21, 21) && ok; // V = G
    ok = test_data_convention(15, 21, 21) && ok;   // V < G
    ok = test_data_convention(15, 21, 15) && ok;   // detector smaller than volume
    ok = test_data_convention(15, 21, 31) && ok;   // detector larger than volume
    ok = test_data_convention(8, 21, 21) && ok;    // even thickness
    // the sample is thin in z inside a thicker box, so little of the z spread
    // reaches the solver padding
    ok = test_end_to_end(15, 21, 21) && ok;
    std::cout << (ok ? "ALL PASS\n" : "SOME FAILED\n");
    return ok ? 0 : 1;
}
