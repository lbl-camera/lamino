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

#include "array.h"

namespace tomocam {

    const float weight[3][3][3] = {
        {{0.0302, 0.037, 0.0302}, {0.037, 0.0523, 0.037}, {0.0302, 0.037, 0.0302}},
        {{0.037, 0.0523, 0.037}, {0.0532, 0., 0.0523}, {0.037, 0.0523, 0.037}},
        {{0.0302, 0.037, 0.0302}, {0.037, 0.0523, 0.037}, {0.0302, 0.037, 0.0302}}
    };

    const double MRF_Q = 2;
    const double MRF_C = 0.001;

    template <typename T>
    T sign(T x) {
        if (x < 0.f) return -1;
        else return 1;
    }

    template <typename T>
    T d_potential_fcn(T delta, T sigma, T p) {

        T Q = static_cast<T>(MRF_Q);
        T C = static_cast<T>(MRF_C);

        T sigma_q   = std::pow(sigma, Q);
        T sigma_q_p = std::pow(sigma, Q - p);

        T temp1 = std::pow(std::abs(delta), Q - p) / sigma_q_p;
        T temp2 = std::pow(std::abs(delta), Q - 1);
        T temp3 = C + temp1;

        if (std::abs(delta) > 0.f) {
            return ((sign(delta) * temp2 / (temp3 * sigma_q)) * (MRF_Q - ((MRF_Q - p) * temp1) /
                (temp3)));
        } else {
            return MRF_Q / (MRF_SIGMA_Q * MRF_C);
        }
    }

    template <typename T>
    void calc_totalvar(Array<T> &input, Array<T> &output, T p, T sigma) {

        // dims
        auto dims = input.dims();

        // write out for clarity
        int nrows = dims.x;
        int ncols = dims.y;

        #pragma omp parallel for
        for (int i = 1; i < nrows - 1; i++)
            for (int j = 1; j < ncols - 1; j++) {
                float v = 0;
                float u = input(i, j);
                for (int x = -1; x < 2; x++)
                    for (int y = -1; y < 2; y++)
                        v += d_potential_fcn(input(i + x, j + y) - u, sigma, p);
                output(i, j) = v;
            }
    }

    // huber function
    template <typename T>
    inline T huber_fcn(T delta, T tau) {
        if (std::abs(delta) < tau)
            return 0.5f * delta * delta;
        else
            return tau * (std::abs(delta) - 0.5 * tau);
    }

    // calculate huber constraints on CPU
    template <typename T>
    void calc_huber_tv(Array<T> &input, Array<T> &output, T tau, T sigma) {

        // dims
        auto dims = input.dims();

        // write out for clarity
        int nrows = dims.x;
        int ncols = dims.y;

        #pragma omp parallel for
        for (int i = 1; i < nrows - 1; i++)
            for (int j = 1; j < ncols - 1; j++) {
                auto u = input(i, j);
                auto du = (input(i + 1, j) - input(i - 1, j) +
                    input(i, j + 1) - input(i, j - 1)) / 4;
                output(i, j) = sigma * huber_fcn(du, tau);
            }
    }
} // namespace tomocam
