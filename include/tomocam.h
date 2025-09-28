#include <vector>

#include "array.h"
#include "dtypes.h"
#include "polar_grid.h"

using std::vector;

#ifndef TOMOCAM__H
#define TOMOCAM__H

namespace tomocam {
    template <typename T>
    Array<T> forward(const Array<T> &, const PolarGrid<T> &, T);

    template <typename T>
    Array<T> backward(const Array<T> &, const PolarGrid<T> &, T);

} // namespace tomocam

#endif // TOMOCAM__H
