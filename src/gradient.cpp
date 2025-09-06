#include <finufft.h>
#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "tomocam.h"
#include "dtypes.h"
#include "fft.h"
#include "nufft.h"
#include "padding.h"

template <typename Real_t>
Array<Real_t> gradient(const Array<Real_t> &f,
    const PolarGrid<Real_t> &pg) {

    // normlization factor
    Real_t scale = (Real_t) std::pow(pg.npixel, 3);

    // cast input to complex
    auto fz = to_complex(f);
    Array<complex<Real_t>> cz(pg.nprojs, pg.npixel);

    // calculate gradient
    nufft::nufft2d2<Real_t>(cz, fz, pg);
    nufft::nufft2d1<Real_t>(cz, fz, pg);

    Array<Real_t> g = real(fz) / scale;
    return g;
}

// explicit instantiation
template Array<float> gradient(const Array<float> &, const PolarGrid<float> &);
template Array<double> gradient(const Array<double> &, const PolarGrid<double> &);
