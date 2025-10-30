#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "fft.h"
#include "fftutils.h"
#include "filter.h"
#include "nufft.h"
#include "padding.h"
#include "polar_grid.h"
#include "tomocam.h"
#include <cstdint>

constexpr double factor = 1.4142135624;

namespace tomocam {

    template <typename T>
    auto forward(const Array<T> &image, const PolarGrid<T> &pg,
        T gamma) -> Array<T> {

        auto F = pad3d(image, static_cast<T>(factor), PadType::SYMMETRIC);
        auto crop_size = F.dims() - image.dims();

        // cast double array to complex
        auto Ft = array::to_complex(F);
        Array<std::complex<T>> C(pg.dims());

        // call NUFFT3d type-2
        nufft::nufft3d2<T>(C, Ft, pg);

        // ifft
        auto C1 = fft::ifftshift(C, fft::Axes::two);
        auto C2 = fft::ifft2(C1);
        auto C3 = fft::fftshift(C2, fft::Axes::two);

        // don't crop projections, just the image
        crop_size.n1 = 0;
        auto C4 = crop2d(C3, crop_size, PadType::SYMMETRIC);
        return array::to_real<T>(C4);
    }

    template <typename T>
    auto backward(const Array<T> &proj, const PolarGrid<T> &pg, T gamma,
        dims_t recon_dims, bool use_filter) -> Array<T> {

        // cast to complex
        auto proj2 = pad2d(proj, static_cast<T>(factor), PadType::SYMMETRIC);
        auto C = array::to_complex(proj2);

        // 2-D Fourier transforms
        auto C1 = fft::ifftshift(C, fft::Axes::two);
        auto C2 = fft::fft2(C1);
        if (use_filter) { apply_filter(C2, "ram-lak"); }
        auto C3 = fft::fftshift(C2, fft::Axes::two);

        size_t n1 = static_cast<size_t>(factor * recon_dims.n1);
        if (n1 % 2 == 0) { n1 -= 1; }
        dims_t new_dims{n1, proj2.nrows(), proj2.ncols()};
        Array<std::complex<T>> F(new_dims);

        // nufft
        nufft::nufft3d1<T>(C3, F, pg);

        // crop image
        auto crop_size = new_dims - recon_dims;
        auto F2 = crop3d(F, crop_size, PadType::SYMMETRIC);
        return array::to_real<T>(F2);
    }

    // Explicit instantiation forward
    template Array<float> forward(const Array<float> &,
        const PolarGrid<float> &, float);
    template Array<double> forward(const Array<double> &,
        const PolarGrid<double> &, double);

    // Explicit instantiation backward
    template Array<float> backward(const Array<float> &,
        const PolarGrid<float> &, float, dims_t, bool);
    template Array<double> backward(const Array<double> &,
        const PolarGrid<double> &, double, dims_t, bool);

} // namespace tomocam
