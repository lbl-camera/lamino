#ifndef TOMOCAM_GPU_TOEPLITZ_HOST_H
#define TOMOCAM_GPU_TOEPLITZ_HOST_H

#include <array>
#include <complex>
#include <vector>

#include "array.h"
#include "dtypes.h"

namespace tomocam::gpu {

    // Host (FINUFFT/FFTW) construction of the 6 unique Fourier-domain PSF
    // kernels, used when the device lacks the memory for cufinufft's fine
    // grid. Kernel order and layout are those of cpu::ToeplitzVectorOp and
    // gpu::ToeplitzVectorOp (dims n1 x n2 x (n3/2+1), R2C).
    // theta/gammas/betas hold one entry per projection; the polar grid has
    // nrows x ncols points per projection.
    template <typename T>
    std::array<Array<std::complex<T>>, 6>
    build_psf_kernels_host(const std::vector<T> &theta, const std::vector<T> &gammas,
                           const std::vector<T> &betas, size_t nrows, size_t ncols,
                           const dims_t &recon_dims);

} // namespace tomocam::gpu

#endif // TOMOCAM_GPU_TOEPLITZ_HOST_H
