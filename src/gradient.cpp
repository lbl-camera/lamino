#include <finufft.h>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "fft.h"
#include "nufft.h"
#include "padding.h"
#include "tomocam.h"

namespace tomocam {

    template <typename T>
    auto gradient(const Array<T> &f, const PolarGrid<T> &pg) -> Array<T> {

        // normlization factor
        T scale = std::pow(static_cast<T>(pg.size()), 3);

        // zero pad
        T factor = static_cast<T>(static_cast<T>(pg.dims().n3)) /
                   static_cast<T>(f.ncols());
        auto f2 = pad3d(f, factor, PadType::SYMMETRIC);

        // cast input to complex
        auto fz = array::to_complex(f2);
        Array<std::complex<T>> cz(pg.dims());

        // calculate gradient
        nufft::nufft3d2<T>(cz, fz, pg);
        nufft::nufft3d1<T>(cz, fz, pg);

        dims_t crop_size = f2.dims() - f.dims();
        crop3d(fz, crop_size, PadType::SYMMETRIC);
        fz /= scale;

        return array::to_real(fz);
    }

    // explicit instantiation
    template Array<float> gradient(const Array<float> &,
        const PolarGrid<float> &);
    template Array<double> gradient(const Array<double> &,
        const PolarGrid<double> &);
} // namespace tomocam
