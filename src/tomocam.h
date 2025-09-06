#include <vector>

#include "array.h"
#include "dtypes.h"

using std::vector;

#ifndef TOMOCAM__H
#define TOMOCAM__H

namespace tomocam {
    template <typename T>
    Array<T> add_padding1d(const Array<T> &A, double oversample) {
        dims_t dims = A.dims();
        int jpad = (oversample - 1) * dims.y / 2;
        Array<T> B(dims.x, dims.y + 2 * jpad, 0);

#pragma omp parallel for
        for (int i = 0; i < dims.x; i++)
            for (int j = 0; j < dims.y; j++) B(i, j + jpad) = A(i, j);
        return B;
    }

    template <typename T>
    Array<T> remove_padding1d(const Array<T> &A, double oversample) {
        dims_t dims = A.dims();
        int jpad = dims.y / oversample / 2;
        Array<T> B(dims.x, dims.y - 2 * jpad);

#pragma omp parallel for
        for (int i = 0; i < B.dims().x; i++)
            for (int j = 0; j < B.dims().y; j++) { B(i, j) = A(i, j + jpad); }
        return B;
    }

    template <typename T>
    Array<T> fftconvolve(const Array<T> &, const Array<T> &);

    template <typename T>
    Array<T> forward(Array<T> &, const vector<T> &);

    template <typename T>
    Array<T> backward(Array<T> &, const vector<T> &);

    template <typename T>
    void ramp_filter(Array<T> &, T);

    template <typename Real_t>
    Array<Real_t> calc_psf(const PolarGrid<Real_t> &);

    template <typename Real_t>
    Array<Real_t> fftconvolve(const Array<Real_t> &, const Array<Real_t> &);

    template <typename Real_t>
    Array<Real_t> gradient(const Array<Real_t> &, const PolarGrid<Real_t> &);

    template <typename Real_t>
    Array<Real_t> gradient2(const Array<Real_t> &f, const Array<Real_t> &psf,
        const Array<Real_t> &gT) {
        Array<Real_t> g = fftconvolve<Real_t>(f, psf);
        g -= gT;
        return g;
    }

} // namespace tomocam

#endif // TOMOCAM__H
