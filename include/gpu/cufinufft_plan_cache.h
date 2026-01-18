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

#ifndef CUFINUFFT_PLAN_CACHE_H
#define CUFINUFFT_PLAN_CACHE_H

#include <array>
#include <mutex>
#include <stdexcept>

#include "cufinufft_plan.h"

namespace tomocam::nufft {

    template <typename T>
    class CufinufftPlanCache {
      private:
        CufinufftPlanWrapper<T> type1_plan_;
        CufinufftPlanWrapper<T> type2_plan_;
        std::once_flag type1_init_flag_;
        std::once_flag type2_init_flag_;

      public:
        CufinufftPlanCache() = default;
        ~CufinufftPlanCache() = default;

        CufinufftPlanCache(const CufinufftPlanCache &) = delete;
        CufinufftPlanCache &operator=(const CufinufftPlanCache &) = delete;
        CufinufftPlanCache(CufinufftPlanCache &&) = delete;
        CufinufftPlanCache &operator=(CufinufftPlanCache &&) = delete;

        CufinufftPlanWrapper<T> &get_plan(int type, int dim,
                                          std::array<int64_t, 3> n_modes, int iflag) {
            if (type == 1) {
                std::call_once(type1_init_flag_, [&]() {
                    type1_plan_.make_plan(1, dim, n_modes, iflag);
                });
                return type1_plan_;
            } else if (type == 2) {
                std::call_once(type2_init_flag_, [&]() {
                    type2_plan_.make_plan(2, dim, n_modes, iflag);
                });
                return type2_plan_;
            } else {
                throw std::invalid_argument(
                    "Only cufinufft type 1 and 2 are supported");
            }
        }
    };

    namespace plans {
        template <typename T>
        inline CufinufftPlanCache<T> cu_cache;
    }

} // namespace tomocam::nufft

#endif // CUFINUFFT_PLAN_CACHE_H
