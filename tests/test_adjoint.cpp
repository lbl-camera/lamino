#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "projection.h"

using namespace tomocam;

static double vec_dot(const std::array<Array<double>, 3> &a,
                      const std::array<Array<double>, 3> &b) {
    double s = 0;
    for (size_t i = 0; i < 3; ++i) s += array::dot<double>(a[i], b[i]);
    return s;
}

static double rel_diff(double a, double b) { return std::abs(a - b) / std::abs(a); }

static bool run_tests(dims_t dims, size_t nangles, double gamma, double beta) {

    std::vector<double> theta(nangles);
    for (size_t i = 0; i < nangles; ++i)
        theta[i] = (static_cast<double>(i) - nangles / 2.0) * M_PI / 180.0;

    auto pg = PolarGrid<double>(theta, dims.n2, dims.n3, gamma, beta);

    std::array<Array<double>, 3> m, u;
    for (size_t i = 0; i < 3; ++i) {
        m[i] = Array<double>::random(dims);
        u[i] = Array<double>::random(dims);
    }

    bool pass = true;

    // Test 1: Adjoint identity  <Rm, y> = sum_i <m_i, (R^T y)_i>
    {
        auto Rm = forward(m, pg, gamma, beta);
        auto y = Array<double>::random(Rm.dims());
        double lhs = array::dot<double>(Rm, y);
        auto RTy = adjoint(y, pg, dims, gamma, beta);
        double rhs = 0;
        for (size_t i = 0; i < 3; ++i) rhs += array::dot<double>(m[i], RTy[i]);
        double rd = rel_diff(lhs, rhs);
        std::cout << std::format(
            "  Test 1 (Adjoint):    lhs={:.6e}  rhs={:.6e}  rel_diff={:.2e}", lhs,
            rhs, rd);
        if (rd < 1e-10) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL\n";
            pass = false;
        }
    }

    // Test 2: Quadratic identity  vec(m, Am) - 2*vec(m, RTy) + ||y||^2 = ||Rm-y||^2
    // Requires odd n2 and n3: the half-pixel frequency offset qX=(k+0.5)*dX-pi gives
    // qX[N-1-k] = -qX[k] exactly for odd N, so the NUFFT output is
    // Hermitian-symmetric for real inputs, making sysmat(m) = adjoint(forward(m)) to
    // floating-point precision.
    {
        auto Rm = forward(m, pg, gamma, beta);
        auto y = Array<double>::random(Rm.dims());
        auto RTy = adjoint(y, pg, dims, gamma, beta);
        auto Am = sysmat(m, pg, gamma, beta);
        double lhs =
            vec_dot(m, Am) - 2.0 * vec_dot(m, RTy) + array::dot<double>(y, y);
        auto res = Rm - y;
        double rhs = array::dot<double>(res, res);
        double rd = rel_diff(lhs, rhs);
        std::cout << std::format(
            "  Test 2 (Quadratic):  lhs={:.6e}  rhs={:.6e}  rel_diff={:.2e}", lhs,
            rhs, rd);
        if (rd < 1e-6) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL\n";
            pass = false;
        }
    }

    // Test 5: per-angle adjoint overload matches explicit overload
    {
        auto y = Array<double>::random(pg.dims());
        auto RTy_explicit = adjoint(y, pg, dims, gamma, beta);
        auto RTy_perangle = adjoint(y, pg, dims);
        double lhs = vec_dot(m, RTy_explicit);
        double rhs = vec_dot(m, RTy_perangle);
        double rd = rel_diff(lhs, rhs);
        std::cout << std::format(
            "  Test 5 (Adj overld): lhs={:.6e}  rhs={:.6e}  rel_diff={:.2e}", lhs,
            rhs, rd);
        if (rd < 1e-14) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL\n";
            pass = false;
        }
    }

    // Test 6: per-angle sysmat overload matches explicit overload
    {
        auto Am_explicit = sysmat(m, pg, gamma, beta);
        auto Am_perangle = sysmat(m, pg);
        double lhs = vec_dot(m, Am_explicit);
        double rhs = vec_dot(m, Am_perangle);
        double rd = rel_diff(lhs, rhs);
        std::cout << std::format(
            "  Test 6 (Sys overld): lhs={:.6e}  rhs={:.6e}  rel_diff={:.2e}", lhs,
            rhs, rd);
        if (rd < 1e-14) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL\n";
            pass = false;
        }
    }

    // Test 3: Sysmat symmetry  vec(u, Am) = vec(m, Au)
    {
        auto Am = sysmat(m, pg, 0.0, 0.0);
        auto Au = sysmat(u, pg, 0.0, 0.0);
        double lhs = vec_dot(u, Am);
        double rhs = vec_dot(m, Au);
        double rd = rel_diff(lhs, rhs);
        std::cout << std::format(
            "  Test 3 (Symmetry):   <u,Am>={:.6e}  <m,Au>={:.6e}  rel_diff={:.2e}",
            lhs, rhs, rd);
        if (rd < 1e-12) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL\n";
            pass = false;
        }
    }

    // Test 4: Sysmat PSD  vec(m, Am) >= 0
    {
        auto Am = sysmat(m, pg, 0.0, 0.0);
        double val = vec_dot(m, Am);
        std::cout << std::format("  Test 4 (PSD):        val={:.6e}", val);
        if (val >= 0) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL\n";
            pass = false;
        }
    }

    return pass;
}

int main() {
    // n2 and n3 must be ODD: the polar grid uses qX = (k+0.5)*dX - pi, which gives
    // qX[N-1-k] = -qX[k] exactly for odd N, so the NUFFT output is
    // Hermitian-symmetric for real inputs and sysmat(m) == adjoint(forward(m)). Even
    // N breaks this: FINUFFT's asymmetric mode range [-N/2, N/2-1] leaves the
    // Nyquist spatial mode (k=-N/2) unpaired, introducing a non-Hermitian
    // contribution.
    struct Geo {
        dims_t dims;
        size_t nangles;
        double gamma;
        double beta;
    };
    Geo geometries[] = {
        {dims_t{5, 15, 15}, 7, M_PI / 4, 0.0},
        {dims_t{8, 31, 31}, 13, M_PI / 4, 0.0},
    };

    bool all_pass = true;
    for (size_t g = 0; g < std::size(geometries); ++g) {
        auto &geo = geometries[g];
        std::cout << std::format(
            "Geometry {} ({}x{}x{}, {} angles, gamma=pi/4, beta=0):\n", g + 1,
            geo.dims.n1, geo.dims.n2, geo.dims.n3, geo.nangles);
        bool ok = run_tests(geo.dims, geo.nangles, geo.gamma, geo.beta);
        all_pass &= ok;
        std::cout << (ok ? "  --> PASS\n" : "  --> FAIL\n") << "\n";
    }

    std::cout << (all_pass ? "Overall: PASS\n" : "Overall: FAIL\n");
    return all_pass ? 0 : 1;
}
