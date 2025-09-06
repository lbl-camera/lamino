#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "fft.h"
#include "fftutils.h"
#include "nufft.h"
#include "padding.h"
#include "polar_grid.h"
#include "tomocam.h"

constexpr double tol = 1.0E-15;
constexpr double crop_factor = 0.4142135624;

namespace tomocam {

    template <typename T>
    Array<T> forward(const Array<T> &image, const PolarGrid<T> &pg) {

        int nslcs = image.nslices();
        int nrows = image.nrows();
        int ncols = image.ncols();
        dims_t pad_size;

        // pad 3d aray
        pad_size.x = 2 * (static_cast<int>(nrows * crop_factor) / 2);
        pad_size.y = 2 * (static_cast<int>(nrows * crop_factor) / 2);
        pad_size.z = 2 * (static_cast<int>(ncols * crop_factor) / 2);
        auto F = pad3d(image, pad_size, PadType::SYMMETRIC);

        // cast double array to complex
        auto Ft = to_complex(F);
        Array<std::complex<T>> C(pg.dims());

        // call NUFFT3d type-2
        nufft::nufft3d2<T>(C, Ft, pg);

        // ifft
        C = fft::ifftshift(C, fft::Axes::two);
        C = fft::ifft2(C);
        C = fft::fftshift(C, fft::Axes::two);

        // crop projection
        pad_size.x = 0;
        C = crop2d(C, pad_size, PadType::SYMMETRIC);
        return real<T>(C);
    }

    template <typename T>
    Array<T> backward(const Array<T> &proj, const PolarGrid<T> &pg) {

        // set switch for double or float
        bool is_double = std::is_same_v<T, double>;

        // dimensions
        dims_t dims = proj.dims();
        int nprojs = dims.x;
        int nrows = dims.y;
        int ncols = dims.z;

        // cast to complex
        auto C = to_complex(proj);

        // pad each projection by factor of sqrt(2)
        dims_t pad_size;
        pad_size.x = 0;
        pad_size.y = 2 * (static_cast<int>(nrows * crop_factor) / 2);
        pad_size.z = 2 * (static_cast<int>(ncols * crop_factor) / 2);
        C = pad2d(C, pad_size, PadType::SYMMETRIC);

        // 2-D Fourier transforms
        C = fft::ifftshift(C, fft::Axes::two);
        C = fft::fft2(C);
        C = fft::fftshift(C, fft::Axes::two);

        Array<std::complex<T>> F(dims_t{nrows, nrows, ncols});

        // nufft
        nufft::nufft3d1<T>(C, F, pg);

        // crop image
        pad_size.x = pad_size.y;
        F = crop3d(F, pad_size, PadType::SYMMETRIC);

        return to_real<T>(F);
    }

    // Explicit instantiation forward
    template Array<float> forward(const Array<float> &,
        const PolarGrid<float> &);
    template Array<double> forward(const Array<double> &,
        const PolarGrid<double> &);

    // Explicit instantiation backward
    template Array<float> backward(const Array<float> &,
        const PolarGrid<float> &);
    template Array<double> backward(const Array<double> &,
        const PolarGrid<double> &);

} // namespace tomocam
