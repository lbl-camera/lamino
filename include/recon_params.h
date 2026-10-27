/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals from the
 * U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 * IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 * the U.S. Government has been granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 * to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */
#ifndef TOMOCAM_RECON_PARAMS_H
#define TOMOCAM_RECON_PARAMS_H

#include <array>
#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "array.h"
#include "logger.h"
#include "padding.h"

namespace tomocam {

    template <typename T>
    struct Dataset_t {
        Array<T> projs;
        std::vector<T> angles;
        T gamma;
        T beta = T(0);
        std::vector<std::array<T, 2>> shifts;
    };

    enum class Regularizer { SPLIT_BREGMAN, UNCONSTRAINED, CONST_MAG, COMBINED };

    struct ReconParams {
        Regularizer regularizer = Regularizer::UNCONSTRAINED;
        std::array<size_t, 3> recon_dims = {0, 0, 0};
        std::string mask_file = ""; // binary support for const_mag; "" = all voxels
        // in-plane {n2, n3} of the mask file; parse_recon_params defaults it to
        // recon_dims[1:2] ({0, 0} is also treated as recon_dims[1:2])
        std::array<size_t, 2> mask_dims = {0, 0};
        size_t maxIters = 100;
        size_t innerIters = 1;
        float lambda = 0.1f;
        float mu = 10.0f;
        float lambda_l2 = 0.1f; // Tikhonov (L2 smoothing) weight, CONST_MAG/COMBINED
        bool const_mag = false; // COMBINED only: enforce |m| = M0
        // magnitude-constraint options (COMBINED); defaults keep rho = mu, no
        // continuation, Ms estimated every iteration, padding left free
        float rho = 0.0f; // 0: use mu
        float rho_growth = 1.0f;
        float rho_max = 0.0f; // 0: unbounded
        float Ms = 0.0f;
        bool estimate_Ms = true;
        bool zero_outside = false;
        float tol = 1e-5f;
        float xtol = 1e-5f;
        float PAD_FACTOR = DEFAULT_PAD_FACTOR;
        LogMode logMode = LogMode::STDOUT;
        std::string logfile = "";

        void print(std::ostream &os) const {
            std::string reg_str = "Unconstrained";
            if (regularizer == Regularizer::SPLIT_BREGMAN) reg_str = "Split-Bregman";
            if (regularizer == Regularizer::CONST_MAG)
                reg_str = "Constant magnitude + L2";
            if (regularizer == Regularizer::COMBINED) reg_str = "Combined";
            os << "Reconstruction Parameters:\n";
            os << "  max_outer_iters: " << maxIters << "\n";
            os << std::format("  recon_dims: [{}, {}, {}]\n", recon_dims[0],
                              recon_dims[1], recon_dims[2]);
            if (!mask_file.empty()) os << "  mask: " << mask_file << "\n";
            if (mask_dims[0] != recon_dims[1] || mask_dims[1] != recon_dims[2])
                os << std::format("  mask_dims: [{}, {}]\n", mask_dims[0],
                                  mask_dims[1]);
            os << "  tol: " << tol << "\n";
            os << "  xtol: " << xtol << "\n";
            os << "  regularizer: " << reg_str << "\n";
            if (regularizer == Regularizer::SPLIT_BREGMAN) {
                os << "    inner_iters: " << innerIters << "\n";
                os << "    lambda: " << lambda << "\n";
                os << "    mu: " << mu << "\n";
            }
            if (regularizer == Regularizer::CONST_MAG) {
                os << "    inner_iters: " << innerIters << "\n";
                os << "    lambda_l2: " << lambda_l2 << "\n";
                os << "    mu: " << mu << "\n";
            }
            if (regularizer == Regularizer::COMBINED) {
                os << "    inner_iters: " << innerIters << "\n";
                os << "    const_mag: " << (const_mag ? "on" : "off") << "\n";
                os << "    lambda_tv: " << lambda << (lambda > 0 ? "" : " (off)")
                   << "\n";
                os << "    lambda_l2: " << lambda_l2
                   << (lambda_l2 > 0 ? "" : " (off)") << "\n";
                os << "    mu: " << mu << "\n";
                if (const_mag) {
                    os << "    rho: " << rho << (rho > 0 ? "" : " (= mu)") << "\n";
                    os << "    rho_growth: " << rho_growth << "\n";
                    os << "    rho_max: " << rho_max
                       << (rho_max > 0 ? "" : " (none)") << "\n";
                    os << "    estimate_Ms: " << (estimate_Ms ? "on" : "off");
                    if (!estimate_Ms) os << " (Ms = " << Ms << ")";
                    os << "\n";
                    os << "    zero_outside: " << (zero_outside ? "on" : "off")
                       << "\n";
                }
            }
        }
    };
} // namespace tomocam

#endif // TOMOCAM_RECON_PARAMS_H
