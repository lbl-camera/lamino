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

#ifndef FINUFFT_PLAN__H
#define FINUFFT_PLAN__H

#include <complex>
#include <finufft.h>
#include <type_traits>

#include "polar_grid.h"

namespace tomocam::nufft {

    template <typename T>
    struct FinufftTraits;

    template <>
    struct FinufftTraits<double> {
        using plan_type = finufft_plan;
        using complex_type = std::complex<double>;

        static int makeplan(int type, int dim, int64_t *n_modes, int iflag,
                            int ntrans, plan_type *plan, finufft_opts *opts) {
            constexpr double TOL = 1e-14;
            return finufft_makeplan(type, dim, n_modes, iflag, ntrans, TOL, plan,
                                    opts);
        }

        static int setpts(plan_type plan, int64_t npts, double *x, double *y,
                          double *z, int nk, double *s, double *t, double *u) {
            return finufft_setpts(plan, npts, x, y, z, nk, s, t, u);
        }

        static int execute(plan_type plan, complex_type *cz, complex_type *fz) {
            return finufft_execute(plan, cz, fz);
        }

        static void destroy(plan_type plan) { finufft_destroy(plan); }
    };

    template <>
    struct FinufftTraits<float> {
        using plan_type = finufftf_plan;
        using complex_type = std::complex<float>;

        static int makeplan(int type, int dim, int64_t *n_modes, int iflag,
                            int ntrans, plan_type *plan, finufft_opts *opts) {
            constexpr float TOL = 1.2e-06f;
            return finufftf_makeplan(type, dim, n_modes, iflag, ntrans, TOL, plan,
                                     opts);
        }

        static int setpts(plan_type plan, int64_t npts, float *x, float *y, float *z,
                          int nk, float *s, float *t, float *u) {
            return finufftf_setpts(plan, npts, x, y, z, nk, s, t, u);
        }

        static int execute(plan_type plan, complex_type *cz, complex_type *fz) {
            return finufftf_execute(plan, cz, fz);
        }

        static void destroy(plan_type plan) { finufftf_destroy(plan); }
    };

    template <typename T>
    class FinufftPlanWrapper {
      private:
        using Traits = FinufftTraits<T>;
        typename Traits::plan_type plan;
        bool initialized = false;

      public:
        FinufftPlanWrapper() = default;

        void make_plan(int type, int dim, std::array<int64_t, 3> n_modes,
                       int iflag) {

            if (dim != n_modes.size()) {
                throw std::runtime_error("FinufftPlanWrapper::make_plan: dim does "
                                         "not match size of n_modes");
            }
            finufft_opts opts;
            finufft_default_opts(&opts);
            opts.upsampfac = 2.0;
            int ierr =
                Traits::makeplan(type, dim, n_modes.data(), iflag, 1, &plan, &opts);
            if (ierr != 0) { throw std::runtime_error("Error in finufft_makeplan"); }
            initialized = true;
        }

        void set_points(const PolarGrid<T> &pg) {
            if (!initialized) {
                throw std::runtime_error(
                    "FinufftPlanWrapper::set_points called before make_plan");
            }
            T *x = (T *)pg.x.begin();
            T *y = (T *)pg.y.begin();
            T *z = (T *)pg.z.begin();
            int ierr =
                Traits::setpts(plan, pg.npts, x, y, z, 0, nullptr, nullptr, nullptr);
            if (ierr != 0) { throw std::runtime_error("Error in finufft_setpts"); }
        }

        int execute(std::complex<T> *cz, std::complex<T> *fz) {
            if (!initialized) {
                throw std::runtime_error(
                    "FinufftPlanWrapper::execute called before make_plan");
            }
            return Traits::execute(plan, (typename Traits::complex_type *)cz,
                                   (typename Traits::complex_type *)fz);
        }

        ~FinufftPlanWrapper() {
            if (initialized) { Traits::destroy(plan); }
        }

        FinufftPlanWrapper(const FinufftPlanWrapper &) = delete;
        FinufftPlanWrapper &operator=(const FinufftPlanWrapper &) = delete;

        bool valid() const { return initialized; }

        FinufftPlanWrapper(FinufftPlanWrapper &&other) noexcept
            : plan(other.plan), initialized(other.initialized) {
            other.initialized = false;
        }

        FinufftPlanWrapper &operator=(FinufftPlanWrapper &&other) noexcept {
            if (this != &other) {
                if (initialized) { Traits::destroy(plan); }
                plan = other.plan;
                initialized = other.initialized;
                other.initialized = false;
            }
            return *this;
        }
    };

} // namespace tomocam::nufft

#endif // FINUFFT_PLAN__H
