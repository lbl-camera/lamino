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
#include <iostream>

#include "array.h"
#include "tomocam.h"

#ifndef OPTIMIZE__H
#define OPTIMIZE__H

namespace tomocam::opt {

    template <typename T, template <typename> class Array, typename Gradient,
        typename Error>
    class Optimizer {
      private:
        Gradient gradient_;
        Error error_;
        int max_inner_iters_;

      public:
        // constructor
        Optimizer(Gradient gradient, Error error) :
            gradient_(gradient), error_(error), max_inner_iters_(20) {}

        Array<T> run(Array<T> sol, int max_iters, T step_size, T tol, T xtol) {

            // initialize
            Array<T> x = sol;
            Array<T> y = sol;
            T t = 1;
            T tnew = 1;
            T step0 = step_size;

            for (int iter = 0; iter < max_iters; iter++) {

                int inner_iter = 0;
                while (true) {

                    // update theta
                    T beta = tnew * (1 / t - 1);
                    tnew =
                        0.5 * (std::sqrt(std::pow(t, 4) + 4 * std::pow(t, 2)) -
                                  std::pow(t, 2));

                    // update y
                    y = sol + (sol - x) * beta;
                    auto g = gradient_(y);

                    // update x
                    sol = y - g * step_size;

                    // check if step size is small enough
                    T fx = error_(sol);
                    T fy = error_(y);
                    T gy = 0.5 * step_size * g.norm();
                    if (fx > (fy + gy)) step_size *= 0.9;
                    else {
                        step_size = step0;
                        t = tnew;
                        x = sol;
                        break;
                    }
                    inner_iter += 1;
                    if (inner_iter >= max_inner_iters_) {
                        std::runtime_error("meh!");
                    }
                }
                T e = error_(sol);
                if (e < tol) {
                    std::cout << "iter: " << iter << ", error: " << e
                              << std::endl;
                    break;
                }
                T xerr = (sol - x).norm();
                if (xerr < xtol) {
                    std::cout << "iter: " << iter << ", error: " << e
                              << std::endl;
                    break;
                }
                std::cout << "iter: " << iter << ", error: " << e
                          << ", xerr: " << xerr << std::endl;
            }
            return sol;
        }
    };

} // namespace tomocam::opt

#endif // TOMOCAM_OPTIMIZE__H
