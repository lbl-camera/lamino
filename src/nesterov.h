/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 *National Laboratory (subject to receipt of any required approvals from the
 *U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 *IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 *the U.S. Government has been granted for itself and others acting on its
 *behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 *to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#include <cmath>

#include "array.h"

#ifndef NESTEROV_OPTIMIZE__H
#define NESTEROV_OPTIMIZE__H
namespace opt {
    class NAGoptimizer {
      private:
        float lipschitz_;
        Array<float> t_;
        Array<float> z_;

      public:
        NAGoptimizer(int num_pixels, int num_projs, const Array<float> & angles, float cen):
            t_(num_pixels, num_pixels), z_(num_pixels, num_pixels) {

            Array<float> x(num_pixels, num_pixels, 1);
            auto y = forward(x, angles, cen);
            auto g = backward(y, angles, cen);
            lipschitz_ = 1./g.max();
        }

        template <typename T>
        void update(Array<T> &recon, Array<T> &gradient) {
            #pragma omp parallel for
            for (int i = 0; i < recon.size(); i++) {
                float z_new = recon[i] - lipschitz_ * gradient[i];
                float t_new = 0.5 * (1 + std::sqrt(1 + 4 * t_[i] * t_[i]));
                float l = (1 - t_[i]) / t_new;
                recon[i] = (1 - l) * z_new + l * z_[i];
                t_[i] = t_new;
                z_[i] = z_new;
            }
        }
    };
} // namespace opt

#endif // NESTEROV_OPTIMIZE__H
