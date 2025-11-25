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

#ifndef FFTPLANS__H
#define FFTPLANS__H

#include <fftw3.h>
#include <memory>
#include <string>
#include <unordered_map>

#include "fft.h"

namespace tomocam::fft {

    template <typename Real>
    class FFTWPlanHandle {
      private:
        std::shared_ptr<typename fftw_traits<Real>::plan> plan_ptr_;

      public:
        FFTWPlanHandle(typename fftw_traits<Real>::plan plan)
            : plan_ptr_(std::make_shared<typename fftw_traits<Real>::plan>(plan)) {}
        ~FFTWPlanHandle() {
            if constexpr (std::is_same<Real, double>::value) {
                fftw_destroy_plan(*plan_ptr_);
            } else {
                fftwf_destroy_plan(*plan_ptr_);
            }
        }
        typename fftw_traits<Real>::plan get() const { return *plan_ptr_; }
    };

    template <typename Real>
    class FFTPlanCache {
      private:
        std::unordered_map<std::string, FFTWPlanHandle<Real>> plan_cache_;

      public:
        FFTPlanCache() = default;
        ~FFTPlanCache() = default;

        FFTPlanCache(const FFTPlanCache &) = delete;
        FFTPlanCache &operator=(const FFTPlanCache &) = delete;
        FFTPlanCache(FFTPlanCache &&) = delete;
        FFTPlanCache &operator=(FFTPlanCache &&) = delete;

        std::string genKey(FFT_Type type, const dims_t &dims, bool is_double) {
            std::string key = std::to_string(static_cast<int>(type));
            key += "_" + std::to_string(dims.n1);
            key += "x" + std::to_string(dims.n2);
            key += "x" + std::to_string(dims.n3);
            key += is_double ? "_double" : "_float";
            return key;
        }

        template <typename Real, FFT_Type type>
        typename fftw_traits<Real>::plan get(const dims_t &dims) {
            std::string key = genKey(type, dims, std::is_same<Real, double>::value);
            auto it = plan_cache_.find(key);
            if (it != plan_cache_.end()) {
                return it->second.get();
            } else {
                // create new plan, add to cache and return it
                typename fftw_traits<Real>::plan new_plan;

                // Determine which plan creation function to use based on type
                if constexpr (FFT_Type::FFT2D || type == FFT_Type::IFFT2D) {
                    new_plan = fftw_traits<Real>::template plan_c2c<type>(dims);
                } else if constexpr (type == FFT_Type::RFFT2D) {
                    new_plan = fftw_traits<Real>::template plan_r2c<type>(dims);
                } else if constexpr (type == FFT_Type::IRFFT2D) {
                    new_plan = fftw_traits<Real>::template plan_c2r<type>(dims);
                } else {
                    throw std::invalid_argument("Unknown FFT_Type");
                }

                plan_cache_[key] = FFTWPlanHandle<Real>(new_plan);
                return new_plan;
            }
        }
    };

    namespace plans {
        inline FFTPlanCache cache;
    } // namespace plans

} // namespace tomocam::fft

#endif // FFTPLANS__H
