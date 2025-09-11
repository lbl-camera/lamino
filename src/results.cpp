#include <fstream>
#include <iostream>

#include "array.h"

namespace tomocam {

    template <typename T>
    T calc_error(const Array<T> &array, const Array<T> &data) {

        T err = 0.f;
        for (int i = 0; i < array.size(); i++)
            err += std::pow(array[i] - data[i], 2);
        return std::sqrt(err);
    }
}
