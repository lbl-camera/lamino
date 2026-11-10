#include <tuple>

#include "gpu/toeplitz_host.h"
#include "polar_grid.h"
#include "toeplitz.h"

namespace tomocam::gpu {

    template <typename T>
    std::array<Array<std::complex<T>>, 6>
    build_psf_kernels_host(const std::vector<T> &theta, const std::vector<T> &gammas,
                           const std::vector<T> &betas, size_t nrows, size_t ncols,
                           const dims_t &recon_dims) {
        std::vector<std::tuple<std::vector<T>, T, T>> agb;
        agb.reserve(theta.size());
        for (size_t i = 0; i < theta.size(); ++i) {
            agb.push_back({std::vector<T>{theta[i]}, gammas[i], betas[i]});
        }
        tomocam::PolarGrid<T> grid(agb, nrows, ncols);
        cpu::ToeplitzVectorOp<T> op(grid, recon_dims);

        // same packing as ToeplitzVectorOp::idx
        std::array<Array<std::complex<T>>, 6> out;
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                int k = i == 0 ? j : (i == 1 ? j + 2 : 5);
                out[k] = op.get(i, j).kernel().clone();
            }
        }
        return out;
    }

    template std::array<Array<std::complex<float>>, 6>
    build_psf_kernels_host(const std::vector<float> &, const std::vector<float> &,
                           const std::vector<float> &, size_t, size_t,
                           const dims_t &);
    template std::array<Array<std::complex<double>>, 6>
    build_psf_kernels_host(const std::vector<double> &, const std::vector<double> &,
                           const std::vector<double> &, size_t, size_t,
                           const dims_t &);

} // namespace tomocam::gpu
