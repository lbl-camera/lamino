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

#ifndef CUFINUFFT_PLAN_H
#define CUFINUFFT_PLAN_H

#include <array>
#include <cuComplex.h>
#include <cufinufft.h>
#include <stdexcept>
#include <type_traits>

#include "polar_grid.h"

namespace tomocam::gpu::nufft {

    template <typename T>
    struct CufinufftTraits;

    template <>
    struct CufinufftTraits<double> {
        using plan_type = cufinufft_plan;
        using complex_type = cuDoubleComplex;

        static int makeplan(int type, int dim, int64_t *n_modes, int iflag,
                            int ntrans, plan_type *plan, cufinufft_opts *opts) {
            constexpr double TOL = 1e-14;
            return cufinufft_makeplan(type, dim, n_modes, iflag, ntrans, TOL, plan,
                                      opts);
        }

        static int setpts(plan_type plan, int64_t npts, double *x, double *y,
                          double *z, int nk, double *s, double *t, double *u) {
            return cufinufft_setpts(plan, npts, x, y, z, nk, s, t, u);
        }

        static int execute(plan_type plan, complex_type *cz, complex_type *fz) {
            return cufinufft_execute(plan, cz, fz);
        }

        static void destroy(plan_type plan) { cufinufft_destroy(plan); }
    };

    template <>
    struct CufinufftTraits<float> {
        using plan_type = cufinufftf_plan;
        using complex_type = cuFloatComplex;

        static int makeplan(int type, int dim, int64_t *n_modes, int iflag,
                            int ntrans, plan_type *plan, cufinufft_opts *opts) {
            constexpr float TOL = 1.2e-06f;
            return cufinufftf_makeplan(type, dim, n_modes, iflag, ntrans, TOL, plan,
                                       opts);
        }

        static int setpts(plan_type plan, int64_t npts, float *x, float *y, float *z,
                          int nk, float *s, float *t, float *u) {
            return cufinufftf_setpts(plan, npts, x, y, z, nk, s, t, u);
        }

        static int execute(plan_type plan, complex_type *cz, complex_type *fz) {
            return cufinufftf_execute(plan, cz, fz);
        }

        static void destroy(plan_type plan) { cufinufftf_destroy(plan); }
    };

    template <typename T>
    class CufinufftPlanWrapper {
      private:
        using Traits = CufinufftTraits<T>;
        typename Traits::plan_type plan;
        bool initialized = false;

      public:
        CufinufftPlanWrapper() = default;

        void make_plan(int type, int dim, std::array<int64_t, 3> n_modes,
                       int iflag) {

            if (dim != n_modes.size()) {
                throw std::runtime_error("CufinufftPlanWrapper::make_plan: dim does "
                                         "not match size of n_modes");
            }
            cufinufft_opts opts;
            cufinufft_default_opts(&opts);
            opts.upsampfac = 2.0;
            int ierr =
                Traits::makeplan(type, dim, n_modes.data(), iflag, 1, &plan, &opts);
            if (ierr != 0) {
                throw std::runtime_error("Error in cufinufft_makeplan");
            }
            initialized = true;
        }

        void set_points(const PolarGrid<T> &pg) {
            if (!initialized) {
                throw std::runtime_error(
                    "CufinufftPlanWrapper::set_points called before make_plan");
            }
            T *x = (T *)pg.x.begin();
            T *y = (T *)pg.y.begin();
            T *z = (T *)pg.z.begin();
            int ierr =
                Traits::setpts(plan, pg.npts, x, y, z, 0, nullptr, nullptr, nullptr);
            if (ierr != 0) { throw std::runtime_error("Error in cufinufft_setpts"); }
        }

        int execute(cuda::std::complex<T> *cz, cuda::std::complex<T> *fz) {
            if (!initialized) {
                throw std::runtime_error(
                    "CufinufftPlanWrapper::execute called before make_plan");
            }
            return Traits::execute(plan, (typename Traits::complex_type *)cz,
                                   (typename Traits::complex_type *)fz);
        }

        ~CufinufftPlanWrapper() {
            if (initialized) { Traits::destroy(plan); }
        }

        CufinufftPlanWrapper(const CufinufftPlanWrapper &) = delete;
        CufinufftPlanWrapper &operator=(const CufinufftPlanWrapper &) = delete;

        bool valid() const { return initialized; }

        CufinufftPlanWrapper(CufinufftPlanWrapper &&other) noexcept
            : plan(other.plan), initialized(other.initialized) {
            other.initialized = false;
        }

        CufinufftPlanWrapper &operator=(CufinufftPlanWrapper &&other) noexcept {
            if (this != &other) {
                if (initialized) { Traits::destroy(plan); }
                plan = other.plan;
                initialized = other.initialized;
                other.initialized = false;
            }
            return *this;
        }
    };

} // namespace tomocam::gpu::nufft

#endif // CUFINUFFT_PLAN__H
