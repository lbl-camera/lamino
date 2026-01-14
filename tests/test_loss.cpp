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
    for (size_t i = 0; i < 3; ++i) { f[i] = Array<float>::random(dims); }

    // Define theta
    size_t ntheta = 141;
    float gamma = 0.7853981634f; // pi/4
    std::vector<float> theta(ntheta, 0.0f);
    for (size_t i = 0; i < ntheta; i++) { theta[i] = (i - 70.f) * M_PI / 180.f; }

    auto pg = PolarGrid(theta, dims.n2, dims.n3, gamma);
    auto y = Array<float>::random(pg.x.dims());
    auto yTy = array::dot(y, y);

    Timer timer;

    // direct method
    timer.start();
    auto diff = forward(f, pg, gamma) - y;
    auto err = array::norm2(diff) / static_cast<float>(y.size());
    timer.stop();
    auto dt1 = timer.seconds();

    // Ajdoint method

    auto yT = adjoint(y, pg, dims, gamma);
    timer.start();
    auto err2 = residual(f, yT, pg, yTy, gamma) / static_cast<float>(y.size());
    timer.stop();
    auto dt2 = timer.seconds();
    std::cout << "Time comparision:\n";
    std::cout << std::format("Direct Method: {}, Adjoint Method: {}\n", dt1, dt2);
    std::cout << "Loss comparision:\n";
    std::cout << std::format("Direct Method Loss: {}, Adjoint Method Loss: {}\n",
                             err, err2);
    std::cout << std::format("Direct Method / Adjoint Method: {}\n", err / err2);
    std::cout << std::format("Adjoint Method / Direct Method: {}\n", err2 / err);
    return 0;
}
