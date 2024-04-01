
#include "array.h"
#include "dtypes.h"

#ifndef TOMOCAM__H
#define TOMOCAM__H
template <typename T>
Array<T> add_padding1d(const Array<T> &A, double oversample) {
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
Array<T> remove_padding1d(const Array <T> &A, double oversample) {
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
Array<T> add_padding2d(const Array<T> &A, double oversample) {
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
Array<T> remove_padding2d(const Array <T> &A, double oversample) {
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

template <typename T>
Array<T> addPadding(const Array <T> &A, const dim2_t dims) {
    Array<T> out(dims.x, dims.y);
    memset(out.ptr(), 0, sizeof(T) * out.size());
    for (int i = 0; i < A.nrows(); i++)
        memcpy(out.row(i), A.row(i), sizeof(T) * A.ncols());
    return out;
}

template <typename T>
Array<T> removePadding(const Array <T> &A, const dim2_t dims) {
    Array<T> out(dims.x, dims.y);
    for (int i = 0; i < dims.x; i++)
        memcpy(out.row(i), A.row(i), sizeof(T) * dims.y);
    return out;
}

// debug functions
Array<complex_t> frwd(Array <double>);
Array<double> back(Array <complex_t>); 

Array<double> forward(Array <double> &, const Array <double> &, double);
Array<double> backward(Array <double> &, const Array <double> &, double);
Array<double> calc_psf(const Array<double> &, int);
void ramp_filter(Array <complex_t> &, double);

Array<double> gradient(const Array<double> &, const Array<double> &, const Array<complex_t> &);
Array<double> fftconvolve(const Array<double> &, const Array<complex_t> &);
Array<double> fftconvolve(Array<double>, Array<double>);
void normalize(Array<double> &, int);
void write_output(const char *, char *, size_t);
double calc_error(const Array<double> &, const Array<double> &);


#endif // TOMOCAM__H
