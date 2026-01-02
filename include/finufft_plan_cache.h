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
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
 * Software to reproduce, distribute copies to the public, prepare derivative
 * works, and perform publicly and display publicly, and to permit other to do
 * so.
 *---------------------------------------------------------------------------------
 */

#ifndef FINUFFT_PLAN_CACHE__H
#define FINUFFT_PLAN_CACHE__H

#include "finufft_plan.h"
#include <mutex>

namespace tomocam::nufft {

    template <typename T>
    class FinufftPlanCache {
      private:
        FinufftPlanWrapper<T> type1_plan_;
        FinufftPlanWrapper<T> type2_plan_;
        std::once_flag type1_init_flag_;
        std::once_flag type2_init_flag_;

      public:
        FinufftPlanCache() = default;
        ~FinufftPlanCache() = default;

        FinufftPlanCache(const FinufftPlanCache &) = delete;
        FinufftPlanCache &operator=(const FinufftPlanCache &) = delete;
        FinufftPlanCache(FinufftPlanCache &&) = delete;
        FinufftPlanCache &operator=(FinufftPlanCache &&) = delete;

        FinufftPlanWrapper<T> &get_plan(int type, int dim, int64_t *n_modes,
                                        int iflag, T tol) {
            if (type == 1) {
                std::call_once(type1_init_flag_, [&]() {
                    type1_plan_.make_plan(1, dim, n_modes, iflag, tol);
                });
                return type1_plan_;
            } else if (type == 2) {
                std::call_once(type2_init_flag_, [&]() {
                    type2_plan_.make_plan(2, dim, n_modes, iflag, tol);
                });
                return type2_plan_;
            } else {
                throw std::invalid_argument(
                    "Only FINUFFT type 1 and 2 are supported");
            }
        }
    };

    namespace plans {
        template <typename T>
        inline FinufftPlanCache<T> cache;
    }

} // namespace tomocam::nufft

#endif // FINUFFT_PLAN_CACHE__H
