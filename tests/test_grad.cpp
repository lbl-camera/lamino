#include <array>
#include <fstream>
#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "nufft.h"
#include "timer.h"
#include "tomocam.h"

using namespace tomocam;
int main() {

    // Define random f
    dims_t dims = {21, 319, 419};
    std::array<Array<float>, 3> f;
    for (size_t i = 0; i < 3; i++) { f[i] = Array<float>::random(dims); }

    // Define theta
    size_t ntheta = 71;
    float gamma = 0.7853981634f; // pi/4
    std::vector<float> theta(ntheta, 0.0f);
    for (size_t i = 0; i < ntheta; i++) { theta[i] = (i - 70.f) * M_PI / 180.f; }
    auto pg = PolarGrid(theta, dims.n2, dims.n3, gamma);

    Timer timer;
    // direct method
    // gamma = 0.0f;                // for testing
    timer.start();
    auto fwd = forward(f, pg, gamma);
    auto grad_direct = adjoint(fwd, pg, dims, gamma);
    timer.stop();
    std::cout << std::format("Direct method time: {:.3f} s\n", timer.seconds());

    // Ajdoint method
    std::array<Array<float>, 3> yT;
    for (size_t i = 0; i < 3; i++) { yT[i] = Array<float>::zeros(fwd.dims()); }

    timer.start();
    auto grad_adj = gradient(f, yT, pg, gamma);
    timer.stop();
    std::cout << std::format("Adjoint method time: {:.3f} s\n", timer.seconds());

    // Compare
    for (size_t i = 0; i < 3; ++i) {
        auto &g1 = grad_direct[i];
        auto &g2 = grad_adj[i];
        auto scale = array::dot<float>(g1, g2) / array::dot<float>(g2, g2);
        auto diff = g1 - g2 * scale;
        std::cout << "Component " << i << ":\n";
        std::cout << "  Scale factor: " << scale << std::endl;
        std::cout << "  Relative difference: "
                  << array::norm2(diff) / array::norm2(g1) << std::endl;
    }

    return 0;
}
