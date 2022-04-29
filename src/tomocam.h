
#include "array.h"
#include "dtypes.h"

#ifndef TOMOCAM__H
#define TOMOCAM__H
template <typename T>
Array<T> add_padding1d(const Array<T> &A, float oversample) {
    dim2_t dims = A.dims();
    int jpad = (oversample-1) * dims.y / 2;
    Array<T> B(dims.x, dims.y + 2 * jpad, 0);

    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
            B(i, j+jpad) = A(i,j);
    return B;
}

template <typename T>
Array<T> remove_padding1d(const Array <T> &A, float oversample) {
    dim2_t dims = A.dims();
    int jpad = dims.y / oversample / 2;
    Array<T> B(dims.x, dims.y - 2 * jpad);

    #pragma omp parallel for
    for (int i = 0; i < B.dims().x; i++)  
        for (int j = 0; j < B.dims().y; j++) {
            B(i,j) = A(i, j+jpad);
		}
    return B;
}

template <typename T>
Array<T> add_padding2d(const Array<T> &A, float oversample) {
    dim2_t dims = A.dims();
    int ipad = (oversample-1) * dims.x / 2;
    int jpad = (oversample-1) * dims.y / 2;
    Array<T> B(dims.x + 2 * ipad, dims.y + 2 * jpad, 0);

    #pragma omp parallel for
    for (int i = 0; i < dims.x; i++)
        for (int j = 0; j < dims.y; j++)
            B(i+ipad, j+jpad) = A(i,j);
    return B;
}

template <typename T>
Array<T> remove_padding2d(const Array <T> &A, float oversample) {
    dim2_t dims = A.dims();
    int ipad = dims.x / oversample / 2;
    int jpad = dims.y / oversample / 2;
    Array<T> B(dims.x - 2*ipad, dims.y - 2*jpad);

    #pragma omp parallel for
    for (int i = 0; i < B.dims().x; i++)  
        for (int j = 0; j < B.dims().y; j++) {
            B(i,j) = A(i+ipad, j+jpad);
		}
    return B;
}


Array<float> forward(Array <float> &, const Array <float> &, float);
Array<float> backward(Array <float> &, const Array <float> &, float);

void ramp_filter(Array <complex_t> &, float);

Array<float> fftconvolve(const Array <float> &, const Array<float> &);
void normalize(Array<float> &, int);
void write_output(const char *, char *, size_t);
float calc_error(const Array<float> &, const Array<float> &);


#endif // TOMOCAM__H
