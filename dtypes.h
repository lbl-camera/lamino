
#ifndef DTYPES__H
#define DTYPES__H

#include <complex>

struct dim3_t {
    int x;
    int y;
    int z;
    dim3_t(int a, int b, int c): x(a), y(b), z(c) {}
};

template<class T>
struct coord_t {
    T x;
    T y;
};

typedef std::complex<double> complex_t;


#endif // DTYPES__H
