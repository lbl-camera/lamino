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

#ifndef MSD3_ADAM__H
#define MSD3_ADAM__H


namespace msd {

    float beta1 = 0.9;
    float beta2 = 0.999;
    float epsilon = 1.0E-08;
    float alpha = 1.0E-03;

    class AdamMinimizer {
        Array<float> beta1t_;
        Array<float> beta2t_;
        Array<float> mean_;
        Array<float> variance_;

      public:
        AdamMinimizer(int nrows, int ncols):
            beta1t_(nrows, ncols, beta1),
            beta2t_(nrows, ncols, beta2),
            mean_(nrows, ncols, 0),
            variance_(nrows, ncols, 0) {}

        Array<float> update(Array<float> grad) {
            mean_ = mean_ * beta1 + (1 - beta1) * grad;
            variance_ = variance_ * beta2 + (1 - beta2) * (grad * grad);
            auto m_hat = mean_ / (1 - beta1t_);
            auto v_hat = variance_ / (1. - beta2t_);
            beta1t_ *= beta1;
            beta2t_ *= beta2;
            return alpha * m_hat / (v_hat.sqrt() + epsilon);
        }
    };
} // namespace msd
#endif
