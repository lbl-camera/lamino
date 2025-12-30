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

#include <string>
#include <unordered_map>

#include "finufft_plan.h"

namespace tomocam::nufft {

    template <typename T>
    class FinufftPlanCache {
      private:
        std::unordered_map<std::string, FinufftPlanWrapper<T>> cache_;

        std::string make_key(int type, int dim, const int64_t *n_modes) {
            std::string key = std::to_string(type) + "_";
            for (int i = 0; i < dim - 1; ++i) {
                key += std::to_string(n_modes[i]) + "x";
            }
            key += std::to_string(n_modes[dim - 1]);
            return key;
        }

      public:
        FinufftPlanCache() = default;
        ~FinufftPlanCache() = default;

        FinufftPlanCache(const FinufftPlanCache &) = delete;
        FinufftPlanCache &operator=(const FinufftPlanCache &) = delete;
        FinufftPlanCache(FinufftPlanCache &&) = delete;
        FinufftPlanCache &operator=(FinufftPlanCache &&) = delete;

        FinufftPlanWrapper<T> &get_plan(int type, int dim, int64_t *n_modes,
                                        int iflag, T tol) {
            // create a hash key from the plan parameters
            std::string key = make_key(type, dim, n_modes);
            auto it = cache_.find(key);
            if (it == cache_.end()) {
                FinufftPlanWrapper<T> plan;
                plan.make_plan(type, dim, n_modes, iflag, tol);
                cache_[key] = std::move(plan);
            }
            return cache_[key];
        }
    };

    namespace plans {
        template <typename T>
        inline FinufftPlanCache<T> cache;
    }

} // namespace tomocam::nufft

#endif // FINUFFT_PLAN_CACHE__H
