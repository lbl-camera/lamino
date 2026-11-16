// Angle-convention regression test.
//
// R = Rx(theta) * Ry(beta) * Rz(gamma), every elementary rotation active and
// counter-clockwise (right-handed) positive. RotationTranspose(theta, gamma,
// beta) must equal the R^T written by alignlsq (alignlsq/geometry.py,
// validation.py: ROTATION_GOLDEN). Keep the golden values in sync with that file.
// beam_dir_vector() must be the third column of RotationTranspose (R^T e_z).

#include <array>
#include <cmath>
#include <cstdio>

#include "projection.h"
#include "rotation.h"

struct Golden {
    double theta, gamma, beta;
    std::array<double, 9> rt;
};

int main() {
    const Golden golden[] = {
        {0.7,
         0.4,
         0.15,
         {0.9107184718850753, 0.3865147523630429, 0.1455960547990535,
          -0.385045594092591, 0.6669767335077371, 0.6378729712326998,
          0.1494381324735992, -0.6369838170354759, 0.756253833972244}},
        {-1.1,
         1.3,
         -0.2,
         {0.2621666615466401, 0.4844284132022276, -0.8346243191136125,
          -0.9443511733327449, -0.04926691242669611, -0.3252285854057239,
          -0.1986693307950612, 0.8734425475223383, 0.4445543984476258}},
        {0.5,
         0.0,
         0.0,
         {1.0, 0.0, 0.0, 0.0, 0.8775825618903728, 0.479425538604203, 0.0,
          -0.479425538604203, 0.8775825618903728}},
    };

    bool ok = true;
    for (const auto &g : golden) {
        auto R = tomocam::RotationTranspose(g.theta, g.gamma, g.beta);
        double err = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                err = std::fmax(err, std::fabs(R[i][j] - g.rt[3 * i + j]));
        auto n = tomocam::beam_dir_vector(g.theta, g.gamma, g.beta);
        double nerr = 0;
        for (int i = 0; i < 3; ++i)
            nerr = std::fmax(nerr, std::fabs(n[i] - R[i][2]));
        bool pass = err < 1e-14 && nerr < 1e-14;
        std::printf("theta=%5.2f gamma=%5.2f beta=%5.2f  |R^T-golden|=%.1e  "
                    "|n-col3|=%.1e  %s\n",
                    g.theta, g.gamma, g.beta, err, nerr, pass ? "PASS" : "FAIL");
        ok &= pass;
    }
    std::puts(ok ? "Overall: PASS" : "Overall: FAIL");
    return ok ? 0 : 1;
}
