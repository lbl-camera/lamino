
#include <algorithm>
#include <cstddef>
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

            size_t i = idx / d.y / d.z;
            size_t j = (idx / d.z) % d.y;
            size_t k = idx % d.z;

            // shift multi-index
            size_t k2 = (k + d.z / 2) % d.z;

            // return flattend shifted index
            return (i * d.y + j) * d.z + k2;
        }

        inline size_t shiftIdx2(dims_t d, size_t idx) {

            size_t i = idx / d.y / d.z;
            size_t j = (idx / d.z) % d.y;
            size_t k = idx % d.z;

            // shift multi-index
            size_t j2 = (j + d.y / 2) % d.y;
            size_t k2 = (k + d.z / 2) % d.z;

            // return flattend shifted index
            return (i * d.y + j2) * d.z + k2;
        }

        inline size_t shiftIdx3(dims_t d, size_t idx) {

            size_t i = idx / d.y / d.z;
            size_t j = (idx / d.z) % d.y;
            size_t k = idx % d.z;

            // shift multi-index
            size_t i2 = (i + d.x / 2) % d.x;
            size_t j2 = (j + d.y / 2) % d.y;
            size_t k2 = (k + d.z / 2) % d.z;

            // return flattend shifted index
            return (i2 * d.y + j2) * d.z + k2;
        }

        inline std::function<size_t(dims_t, size_t)> select_shifter(Axes axis) {
            if (axis == Axes::one) {
                return std::function<size_t(dims_t, size_t)>(shiftIdx1);
            } else if (axis == Axes::two) {
                return std::function<size_t(dims_t, size_t)>(shiftIdx2);
            } else if (axis == Axes::three) {
                return std::function<size_t(dims_t, size_t)>(shiftIdx3);
            } else {
                std::runtime_error("unknown shift axes");
            }
        }

        template <typename T>
        Array<T> fftshift(const Array<T> &input, Axes axes) {

            auto dims = input.dims();
            Array<T> output(dims);

            auto shifter = select_shifter(axes);

            size_t nelems = input.size();
            std::for_each(std::execution::par_unseq,
                std::views::iota(size_t(0), nelems), [&](size_t idx) {
                    auto shifted_idx = shifter(dims, idx);
                    output[shifted_idx] = input[idx];
                });
            return output;
        }

        template <typename T>
        Array<T> ifftshift(const Array<T> &input, Axes axes) {

            auto dims = input.dims();
            Array<T> output(dims);

            auto shifter = select_shifter(axes);

            size_t nelems = input.size();
            std::for_each(std::execution::par_unseq,
                std::views::iota(size_t(0), nelems), [&](size_t idx) {
                    auto shifted_idx = shifter(dims, idx);
                    output[idx] = input[shifted_idx];
                });

            return output;
        }
    } // namespace fft
} // namespace tomocam

#endif // FFTUTILS__H
