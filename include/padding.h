
#include "array.h"
#include "dtypes.h"
#include <cstdint>

#ifndef PADDING__H
#define PADDING__H

namespace tomocam {
    enum class PadType { LEFT, RIGHT, SYMMETRIC };

    template <typename T> Array<T> pad1d(const Array<T> &arr, uint64_t pad_size, PadType pad_type) {

        if (pad_size == 0) {
            return arr;
        }

        auto dims = arr.dims();
        Array<T> arr2(dims.x(), dims.y(), dims.z() + pad_size);

        // assume symmetric padding by default
        int d = pad_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = 0;
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (int i = 0; i < dims.x(); i++)
            for (int j = 0; j < dims.y(); j++)
                for (int k = 0; k < dims.z(); k++)
                    arr2[i, j, k + d] = arr[i, j, k];
        return arr2;
    }

    template <typename T>
    Array<T> crop1d(const Array<T> &arr, uint64_t crop_size, PadType pad_type) {

        if (crop_size == 0) {
            return arr;
        }
        auto dims = arr.dims() - dims_t{0, 0, crop_size};
        Array<T> arr2(dims);

        // assume symmetric padding by default
        int d = crop_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = 0;
        } else if (pad_type == PadType::LEFT) {
            d = crop_size;
        }

        for (int i = 0; i < dims.x(); i++)
            for (int j = 0; j < dims.y(); j++)
                for (int k = 0; k < dims.z(); k++)
                    arr2[i, j, k] = arr[i, j, k + d];
        return arr2;
    }

    template <typename T> Array<T> pad2d(const Array<T> &arr, dims_t pad_size, PadType pad_type) {

        if ((pad_size.y() == 0) && (pad_size.z() == 0)) {
            return arr;
        }

        // create and initialize return array
        auto dims = arr.dims();
        dims_t new_dims = dims + pad_size;
        Array<T> arr2(new_dims);
        arr2.fill(static_cast<T>(0));

        // assume symmetric padding by default
        dims_t d = pad_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = dims_t(0, 0, 0);
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (int i = 0; i < dims.x(); i++)
            for (int j = 0; j < dims.y(); j++)
                for (int k = 0; k < dims.z(); k++)
                    arr2[i, j + d.y(), k + d.z()] = arr[i, j, k];
        return arr2;
    }

    template <typename T> Array<T> crop2d(const Array<T> &arr, dims_t crop_size, PadType pad_type) {

        if ((crop_size.y() == 0) && (crop_size.z() == 0)) {
            return arr;
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

        for (int i = 0; i < new_dims.x(); i++) {
            for (int j = 0; j < new_dims.y(); j++)
                for (int k = 0; k < new_dims.z(); k++)
                    arr2[i, j, k] = arr[i, j + d.y(), k + d.z()];
        }
        return arr2;
    }

    template <typename T> Array<T> pad3d(const Array<T> &arr, dims_t pad_size, PadType pad_type) {

        if ((pad_size.x() == 0) && (pad_size.y() == 0) && (pad_size.z() == 0))
            return arr;

        auto dims = arr.dims();
        dims_t new_dims = dims + pad_size;
        Array<T> arr2(new_dims);

        // assume symmetric padding
        dims_t d = pad_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = dims_t(0, 0, 0);
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (int i = 0; i < dims.x(); i++)
            for (int j = 0; j < dims.y(); j++)
                for (int k = 0; k < dims.z(); k++)
                    arr2[i + d.x(), j + d.y(), k + d.z()] = arr[i, j, k];
        return arr2;
    }

    template <typename T> Array<T> crop3d(const Array<T> &arr, dims_t crop_size, PadType pad_type) {

        if ((crop_size.x() == 0) && (crop_size.y() == 0) && (crop_size.z() == 0))
            return arr;

        auto dims = arr.dims();
        dims_t new_dims = dims - crop_size;
        Array<T> arr2(new_dims);

        dims_t d = crop_size / 2;
        if (pad_type == PadType::LEFT)
            d = dims_t(0, 0, 0);
        else if (pad_type == PadType::RIGHT)
            d = crop_size;

        for (int i = 0; i < new_dims.x(); i++)
            for (int j = 0; j < new_dims.y(); j++)
                for (int k = 0; k < new_dims.z(); k++)
                    arr2[i, j, k] = arr[i + d.x(), j + d.y(), k + d.z()];
        return arr2;
    }
} // namespace tomocam
#endif // PADDING__H
