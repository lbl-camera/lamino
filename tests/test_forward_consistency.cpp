// Consistency between the forward model used to simulate XMCD projections
// (forward.cpp -> tomocam::forward) and the model inverted by recon
// (adjoint + Toeplitz A^T A).
//
//  1. padded_dim() is always odd and pad2d/pad3d agree with it, so forward's
//     padded volume and recon's padded detector grid have identical sizes.
//  2. On an odd grid, forward() at theta = gamma = beta = 0 equals the exact
//     line integral of m_z along n1 (divided by n1, per convention).
//  3. On an odd grid, the Toeplitz normal operator equals adjoint(forward(x)).
//
// For even grids the half-sample-offset polar grid has no q = 0 sample, so (2)
// and (3) fail; those cases are printed for reference but not asserted.

#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <random>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "config.h"
#include "dtypes.h"
#include "padding.h"
#include "polar_grid.h"
#include "projection.h"
#include "toeplitz.h"

using namespace tomocam;
using vec3 = std::array<Array<double>, 3>;

static Array<double> random_array(const dims_t &dims, std::mt19937 &gen) {
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    Array<double> a(dims);
    for (size_t i = 0; i < a.size(); ++i) a[i] = dis(gen);
    return a;
}

static double vec_dot(const vec3 &a, const vec3 &b) {
    double s = 0;
    for (size_t i = 0; i < 3; ++i) s += array::dot<double>(a[i], b[i]);
    return s;
}

// random magnetization, zero outside the central (unpadded) support
static vec3 random_padded_volume(size_t n1, size_t n, std::mt19937 &gen) {
    vec3 m;
    for (size_t i = 0; i < 3; ++i) {
        m[i] = pad3d<double>(random_array(dims_t{n1, n, n}, gen), DEFAULT_PAD_FACTOR,
                             PadType::SYMMETRIC);
    }
    return m;
}

static bool test_padded_sizes() {
    bool ok = true;
    // start at 5: padded_dim throws for n = 2, 4 (no room for odd padding)
    for (size_t n = 5; n <= 2048; ++n) {
        size_t N = padded_dim(n, DEFAULT_PAD_FACTOR);
        Array<double> a2(dims_t{1, n, n});
        auto p2 = pad2d<double>(a2, DEFAULT_PAD_FACTOR, PadType::SYMMETRIC);
        bool good = (N % 2 == 1) && N >= n && p2.nrows() == N && p2.ncols() == N;
        if (n <= 64) {
            Array<double> a3(dims_t{n, n, n});
            auto p3 = pad3d<double>(a3, DEFAULT_PAD_FACTOR, PadType::SYMMETRIC);
            good = good && p3.nslices() == N && p3.nrows() == N && p3.ncols() == N;
        }
        if (!good) {
            std::cout << std::format("  padded_dim({}) = {} inconsistent\n", n, N);
            ok = false;
        }
    }
    std::cout << std::format("padded sizes odd and consistent: {}\n",
                             ok ? "PASS" : "FAIL");
    return ok;
}

// returns relative error of forward() vs the exact projection of m_z
static double straight_projection_error(size_t N) {
    std::mt19937 gen(7);
    dims_t dims{9, N, N};
    vec3 m;
    for (size_t i = 0; i < 3; ++i) m[i] = random_array(dims, gen);

    PolarGrid<double> pg(std::vector<double>{0.0}, N, N, 0.0, 0.0);
    auto p = forward<double>(m, pg, 0.0, 0.0);

    // FFTW ifft2 is unnormalized, so forward() = sum along n1 / n1
    double scale = static_cast<double>(dims.n1);
    double err = 0, ref = 0;
    for (size_t j = 0; j < N; ++j) {
        for (size_t k = 0; k < N; ++k) {
            double s = 0;
            for (size_t i = 0; i < dims.n1; ++i) s += m[2][{i, j, k}];
            s /= scale;
            double d = p[{0, j, k}] - s;
            err += d * d;
            ref += s * s;
        }
    }
    return std::sqrt(err / ref);
}

// relative error between Toeplitz A^T A x and adjoint(forward(x)), after
// fitting the (convention-dependent) global scale
static double normal_operator_error(size_t n, double gamma) {
    std::mt19937 gen(42);
    size_t N = padded_dim(n, DEFAULT_PAD_FACTOR);
    auto x = random_padded_volume(5, n, gen);
    dims_t vol = x[0].dims();

    std::vector<double> theta;
    for (int i = 0; i < 12; ++i) theta.push_back(-M_PI / 3 + i * M_PI / 18);
    PolarGrid<double> pg(theta, N, N, gamma, 0.0);

    auto ata =
        adjoint<double>(forward<double>(x, pg, gamma, 0.0), pg, vol, gamma, 0.0);
    cpu::ToeplitzVectorOp<double> op(pg, vol);
    auto toe = cpu::sysmat(x, op);

    double s = vec_dot(ata, toe) / vec_dot(toe, toe);
    double err = 0;
    for (size_t i = 0; i < 3; ++i) {
        auto d = ata[i] - toe[i] * s;
        err += array::dot<double>(d, d);
    }
    return std::sqrt(err / vec_dot(ata, ata));
}

// forward.cpp and recon share to_radians_if_degrees; check it on the cases
// the old per-program rules disagreed on
static bool test_degree_detection() {
    auto conv = [](std::vector<float> a) {
        to_radians_if_degrees(a);
        return a;
    };
    const float d2r = static_cast<float>(M_PI) / 180.0f;
    bool ok = true;
    // all-negative degrees (old recon rule used signed max -> left as radians)
    auto a = conv({-180.0f, -90.0f, -1.0f});
    ok = ok && std::abs(a[0] + 180.0f * d2r) < 1e-6f && std::abs(a[2] + d2r) < 1e-6f;
    // radians up to 2*pi (old forward rule |a| > pi -> treated as degrees)
    auto b = conv({0.0f, 3.5f, 6.0f});
    ok = ok && b[1] == 3.5f && b[2] == 6.0f;
    // ordinary symmetric degree range
    auto c = conv({-90.0f, 0.0f, 89.0f});
    ok = ok && std::abs(c[0] + 90.0f * d2r) < 1e-6f;
    std::cout << std::format("degree detection shared rule: {}\n",
                             ok ? "PASS" : "FAIL");
    return ok;
}

// MBIR2 backprojects all datasets at once on a unified grid; with per-dataset
// COR shifts stacked in the same order, it must equal the sum of per-dataset
// adjoints that each apply their own shifts
static bool test_unified_adjoint_shifts() {
    size_t N = 15;
    dims_t vol{5, N, N};
    std::vector<double> th1, th2;
    for (int i = 0; i < 6; ++i) th1.push_back(-M_PI / 4 + i * M_PI / 12);
    for (int i = 0; i < 4; ++i) th2.push_back(-M_PI / 6 + i * M_PI / 9);
    double g1 = 0.0, g2 = M_PI / 4;

    std::mt19937 gen(3);
    std::uniform_real_distribution<double> dis(-2.0, 2.0);
    std::vector<std::array<double, 2>> s1, s2;
    for (size_t i = 0; i < th1.size(); ++i) s1.push_back({dis(gen), dis(gen)});
    for (size_t i = 0; i < th2.size(); ++i) s2.push_back({dis(gen), dis(gen)});

    auto y1 = random_array(dims_t{th1.size(), N, N}, gen);
    auto y2 = random_array(dims_t{th2.size(), N, N}, gen);

    PolarGrid<double> pg1(th1, N, N, g1, 0.0), pg2(th2, N, N, g2, 0.0);
    auto ref1 = adjoint<double>(y1, pg1, vol, g1, 0.0, s1);
    auto ref2 = adjoint<double>(y2, pg2, vol, g2, 0.0, s2);

    PolarGrid<double> pg({{th1, g1, 0.0}, {th2, g2, 0.0}}, N, N);
    Array<double> y(dims_t{th1.size() + th2.size(), N, N});
    std::copy(y1.begin(), y1.end(), y.begin());
    std::copy(y2.begin(), y2.end(), y.begin() + y1.size());
    auto shifts = s1;
    shifts.insert(shifts.end(), s2.begin(), s2.end());
    auto uni = adjoint<double>(y, pg, vol, shifts);
    auto uni_noshift = adjoint<double>(y, pg, vol);

    double err = 0, ref = 0, diff_noshift = 0;
    for (size_t i = 0; i < 3; ++i) {
        auto r = ref1[i] + ref2[i];
        auto d = uni[i] - r;
        auto d0 = uni_noshift[i] - r;
        err += array::dot<double>(d, d);
        diff_noshift += array::dot<double>(d0, d0);
        ref += array::dot<double>(r, r);
    }
    err = std::sqrt(err / ref);
    diff_noshift = std::sqrt(diff_noshift / ref);
    bool ok = err < 1e-10 && diff_noshift > 1e-2;
    std::cout << std::format(
        "unified adjoint with shifts vs per-dataset: err={:.2e} (without shifts "
        "{:.2e}) {}\n",
        err, diff_noshift, ok ? "PASS" : "FAIL");
    return ok;
}

int main() {
    bool ok = test_padded_sizes();
    ok = test_degree_detection() && ok;
    ok = test_unified_adjoint_shifts() && ok;

    for (size_t N : {21, 25, 20, 24}) {
        double e = straight_projection_error(N);
        bool odd = N % 2 == 1;
        bool pass = e < 1e-8;
        std::cout << std::format(
            "forward vs exact projection, N={:<3}: err={:.2e} {}\n", N, e,
            odd ? (pass ? "PASS" : "FAIL") : "(even, info)");
        if (odd) ok = ok && pass;
    }

    // n=14 pads to 19 (odd, via padded_dim); compare with the old even size 18
    for (double gamma : {0.0, M_PI / 4}) {
        double e = normal_operator_error(14, gamma);
        bool pass = e < 1e-6;
        std::cout << std::format(
            "Toeplitz A^T A vs adjoint(forward), gamma={:.2f}: err={:.2e} {}\n",
            gamma, e, pass ? "PASS" : "FAIL");
        ok = ok && pass;
    }

    std::cout << (ok ? "ALL PASSED\n" : "SOME TESTS FAILED\n");
    return ok ? 0 : 1;
}
