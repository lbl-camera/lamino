
#include "array.h"
#include "dtypes.h"
#include <cmath>
#include <cstdint>
#include <system_error>

#ifndef PADDING__H
    #define PADDING__H

namespace tomocam {
    enum class PadType { LEFT, RIGHT, SYMMETRIC };

    template <typename T>
    auto pad1d(const Array<T> &arr, T factor, PadType pad_type) -> Array<T> {

        if ((int)std::round(factor) == 1) { return arr.clone(); }

        size_t n3 = static_cast<size_t>(factor * arr.ncols());
        if (n3 % 2 == 0) { n3 -= 1; }
        size_t pad_size = n3 - arr.ncols();

        // allocate return array
        dims_t dims = {arr.nslices(), arr.nrows(), n3};
        Array<T> arr2(dims);

        // assume symmetric padding by default
        size_t d = pad_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = 0;
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (size_t i = 0; i < dims.x(); i++)
            for (size_t j = 0; j < dims.y(); j++)
                for (size_t k = 0; k < dims.z(); k++)
                    arr2[{i, j, k + d}] = arr[{i, j, k}];
        return arr2;
    }

    template <typename T>
    auto crop1d(const Array<T> &arr, size_t crop_size,
        PadType pad_type) -> Array<T> {

        if (crop_size == 0) { return arr.clone(); }
        auto dims = arr.dims() - dims_t{0, 0, crop_size};
        Array<T> arr2(dims);

        // assume symmetric padding by default
        int d = crop_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = 0;
        } else if (pad_type == PadType::LEFT) {
            d = crop_size;
        }

        for (size_t i = 0; i < dims.x(); i++)
            for (size_t j = 0; j < dims.y(); j++)
                for (size_t k = 0; k < dims.z(); k++)
                    arr2[{i, j, k}] = arr[{i, j, k + d}];
        return arr2;
    }

    template <typename T>
    auto pad2d(const Array<T> &arr, T factor, PadType pad_type) -> Array<T> {

        if ((int)std::round(factor) == 1) { return arr.clone(); }
        auto n2 = static_cast<size_t>(factor * arr.nrows());
        if (n2 % 2 == 0) { n2 -= 1; }
        auto n3 = static_cast<size_t>(factor * arr.ncols());
        if (n3 % 2 == 0) { n3 -= 1; }

        // create and initialize return array
        dims_t dims{arr.nslices(), n2, n3};
        Array<T> arr2(dims);

        // assume symmetric padding by default
        auto pad_size = dims - arr.dims();
        dims_t d = pad_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = dims_t(0, 0, 0);
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (size_t i = 0; i < dims.x(); i++) {
            for (size_t j = 0; j < dims.y(); j++) {
                for (size_t k = 0; k < dims.z(); k++) {
                    arr2[{i, j + d.y(), k + d.z()}] = arr[{i, j, k}];
                }
            }
        }
        return arr2;
    }

    template <typename T>
    auto crop2d(const Array<T> &arr, dims_t crop_size,
        PadType pad_type) -> Array<T> {

        if ((crop_size.y() == 0) && (crop_size.z() == 0)) {
            return arr.clone();
        }

        auto dims = arr.dims();
        dims_t new_dims = dims - crop_size;
        Array<T> arr2(new_dims);

        // assume symmetric padding by default
        dims_t d = crop_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = dims_t(0, 0, 0);
        } else if (pad_type == PadType::LEFT) {
            d = crop_size;
        }

        for (size_t i = 0; i < new_dims.x(); i++) {
            for (size_t j = 0; j < new_dims.y(); j++) {
                for (size_t k = 0; k < new_dims.z(); k++) {
                    arr2[{i, j, k}] = arr[{i, j + d.y(), k + d.z()}];
                }
            }
        }
        return arr2;
    }

    template <typename T>
    auto pad3d(const Array<T> &arr, T factor, PadType pad_type) -> Array<T> {

        if ((int)std::round(factor) == 1) { return arr.clone(); }

        // calculate new dims
        auto n1 = static_cast<size_t>(static_cast<T>(arr.nslices() * factor));
        if (n1 % 2 == 0) { n1 -= 1; }
        auto n2 = static_cast<size_t>(static_cast<T>(arr.nrows() * factor));
        if (n2 % 2 == 0) { n2 -= 1; }
        auto n3 = static_cast<size_t>(static_cast<T>(arr.ncols() * factor));
        if (n3 % 2 == 0) { n3 -= 1; }

        // allocate return array
        dims_t dims{n1, n2, n3};
        Array<T> arr2(dims);

        auto pad_size = dims - arr.dims();

        // assume symmetric padding
        dims_t d = pad_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = dims_t(0, 0, 0);
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (size_t i = 0; i < dims.x(); i++) {
            for (size_t j = 0; j < dims.y(); j++) {
                for (size_t k = 0; k < dims.z(); k++) {
                    arr2[{i + d.x(), j + d.y(), k + d.z()}] = arr[{i, j, k}];
                }
            }
        }
        return arr2;
    }

    template <typename T>
    Array<T> crop3d(const Array<T> &arr, dims_t crop_size, PadType pad_type) {

        if ((crop_size.x() == 0) && (crop_size.y() == 0) &&
            (crop_size.z() == 0))
            return arr.clone();

        auto dims = arr.dims();
        dims_t new_dims = dims - crop_size;
        Array<T> arr2(new_dims);

        dims_t d = crop_size / 2;
        if (pad_type == PadType::LEFT) d = dims_t(0, 0, 0);
        else if (pad_type == PadType::RIGHT)
            d = crop_size;

        for (size_t i = 0; i < new_dims.x(); i++) {
            for (size_t j = 0; j < new_dims.y(); j++) {
                for (size_t k = 0; k < new_dims.z(); k++) {
                    arr2[{i, j, k}] = arr[{i + d.x(), j + d.y(), k + d.z()}];
                }
            }
        }
        return arr2;
    }
} // namespace tomocam
#endif // PADDING__H
