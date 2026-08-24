#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "projection.h"

using namespace tomocam;

int main() {
    // n2=539 and n3=359 are both odd, so the NUFFT grid is Hermitian-symmetric
    // for real inputs and sysmat(m) == adjoint(forward(m)).
    dims_t dims = {21, 539, 359};
    std::array<Array<double>, 3> f;
    for (size_t i = 0; i < 3; i++) { f[i] = Array<double>::random(dims); }

    size_t ntheta = 71;
    double gamma = M_PI / 4.0;
    std::vector<double> theta(ntheta);
    for (size_t i = 0; i < ntheta; i++) { theta[i] = (i - 70.0) * M_PI / 180.0; }
    auto pg = PolarGrid(theta, dims.n2, dims.n3, gamma);

    bool pass = true;

    // Test 1: adjoint dot-product test  <Af, y> == <f, A^T y>
    // Verifies that adjoint() is the true mathematical adjoint of forward().
    {
        auto Af = forward(f, pg, gamma, 0.0);
        auto y = Array<double>::random(Af.dims());

        double lhs = array::dot<double>(Af, y);

        auto ATy = adjoint(y, pg, dims, 0.0);
        double rhs = 0.0;
        for (size_t i = 0; i < 3; ++i) { rhs += array::dot<double>(f[i], ATy[i]); }

        double rel_diff = std::abs(lhs - rhs) / std::abs(lhs);
        std::cout << std::format(
            "Adjoint test: <Af,y> = {:.6e}  <f,A^Ty> = {:.6e}  rel_diff = {:.2e}\n",
            lhs, rhs, rel_diff);

        if (rel_diff < 1e-10) {
            std::cout << "  PASS\n";
        } else {
            std::cout << "  FAIL\n";
            pass = false;
        }
    }

    // Test 2: Gradient equivalence  R^T(Rm - y)  vs  sysmat(m) - R^T y
    // Both express the gradient of 0.5*||Rm-y||^2 in two equivalent ways.
    // Requires odd n2, n3 so that sysmat(m) == adjoint(forward(m)).
    {
        auto Rm  = forward(f, pg, gamma, 0.0);
        auto y   = Array<double>::random(Rm.dims());
        auto res = Rm - y;
        auto grad_direct = adjoint(res, pg, dims, 0.0);

        auto Am  = sysmat(f, pg, 0.0);
        auto RTy = adjoint(y, pg, dims, 0.0);
        std::array<Array<double>, 3> grad_expanded;
        for (size_t i = 0; i < 3; ++i) grad_expanded[i] = Am[i] - RTy[i];

        for (size_t i = 0; i < 3; ++i) {
            auto &g1 = grad_direct[i];
            auto &g2 = grad_expanded[i];
            double scale = array::dot<double>(g1, g2) / array::dot<double>(g2, g2);
            auto diff = g1 - g2 * scale;
            double rel_diff = array::norm2<double>(diff) / array::norm2<double>(g1);
            std::cout << std::format(
                "Gradient component {}: scale = {:.6f}  rel_diff = {:.2e}\n",
                i, scale, rel_diff);

            if (rel_diff < 1e-9) {
                std::cout << "  PASS\n";
            } else {
                std::cout << "  FAIL\n";
                pass = false;
            }
        }
    }

    return pass ? 0 : 1;
}
