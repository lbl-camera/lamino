
#ifndef DTYPES__H
#define DTYPES__H

#include <complex>

struct dim2_t {
    int x;
    int y;
};

typedef std::complex <float >complex_t;
inline const complex_t I(0, 1);
#endif // DTYPES__H
