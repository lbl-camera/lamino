
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <functional>
#include <ranges>
#include <stdexcept>

#include "array.h"

#ifndef FFTUTILS__H
    #define FFTUTILS__H

namespace tomocam {
    namespace fft {

        enum class Axes { one, two, three };

        inline size_t shiftIdx1(dims_t d, size_t idx) {

            auto [i, j, k] = d.unravel_idx(idx);

            // shift multi-index
            size_t k2 = (k + d.z() / 2) % d.z();

            // return flattend shifted index
            return d.flat_idx(i, j, k2);
        }

        inline size_t shiftIdx2(dims_t d, size_t idx) {

            auto [i, j, k] = d.unravel_idx(idx);

            // shift multi-index
            size_t j2 = (j + d.y() / 2) % d.y();
            size_t k2 = (k + d.z() / 2) % d.z();

            // return flattend shifted index
            return d.flat_idx(i, j2, k2);
        }

        inline size_t shiftIdx3(dims_t d, size_t idx) {

            auto [i, j, k] = d.unravel_idx(idx);

            // shift multi-index
            size_t i1 = (i + d.x() / 2) % d.x();
            size_t j1 = (j + d.y() / 2) % d.y();
            size_t k1 = (k + d.z() / 2) % d.z();

            // return flattend shifted index
            return d.flat_idx(i1, j1, k1);
        }

        inline std::function<size_t(dims_t, size_t)> select_shifter(Axes axis) {
            if (axis == Axes::one) {
                return std::function<size_t(dims_t, size_t)>(shiftIdx1);
            } else if (axis == Axes::two) {
                return std::function<size_t(dims_t, size_t)>(shiftIdx2);
            } else if (axis == Axes::three) {
                return std::function<size_t(dims_t, size_t)>(shiftIdx3);
            } else {
                throw std::runtime_error("unknown shift axes");
            }
            throw std::runtime_error("booo!");
        }

        template <typename T>
        Array<T> fftshift(const Array<T> &input, Axes axes) {

            auto dims = input.dims();
            Array<T> output(dims);

            auto shifter = select_shifter(axes);
            size_t nelems = input.size();
            for (size_t i0 = 0; i0 < nelems; ++i0) {
                auto i1 = shifter(dims, i0);
                output[i1] = input[i0];
            }
            return output;
        }

        template <typename T>
        Array<T> ifftshift(const Array<T> &input, Axes axes) {

            auto dims = input.dims();
            Array<T> output(dims);

            auto shifter = select_shifter(axes);
            size_t nelems = input.size();
            for (size_t i0 = 0; i0 < nelems; i0++) {
                auto i1 = shifter(dims, i0);
                output[i0] = input[i1];
            }
            return output;
        }
    } // namespace fft
} // namespace tomocam

#endif // FFTUTILS__H
