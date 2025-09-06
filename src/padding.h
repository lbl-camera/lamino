
#include "array.h"

#ifndef PADDING__H
#define PADDING__H

namespace tomocam {
    enum class PadType { LEFT, RIGHT, SYMMETRIC };

    template <typename T>
    Array<T> pad1d(const Array<T> &arr, int pad_size, PadType pad_type) {

        if (pad_size == 0) { return arr; }

        auto dims = arr.dims();
        Array<T> arr2(dims.x, dims.y + pad_size);

        // assume symmetric padding by default
        int d = pad_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = 0;
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (int i = 0; i < dims.x; i++) {
            for (int j = 0; j < dims.y; j++) new_arr(i, j + d) = arr(i, j);
        }
        return new_arr;
    }

    template <typename T>
    Array<T> crop1d(const Array<T> &arr, int crop_size, PadType pad_type) {

        if (crop_size == 0) { return arr; }
        auto dims = arr.dims();
        Array<T> new_arr(dims.x, dims.y - crop_size);

        // assume symmetric padding by default
        int d = crop_size / 2;
        if (pad_type == PadType::RIGHT) {
            d = 0;
        } else if (pad_type == PadType::LEFT) {
            d = crop_size;
        }

        for (int i = 0; i < dims.x; i++) {
            for (int j = 0; j < dims.y - crop_size; j++)
                new_arr(i, j) = arr(i, j + d);
        }
        return new_arr;
    }

    template <typename T>
    Array<T> pad2d(const Array<T> &arr, dims_t pad_size, PadType pad_type) {

        if ((pad_size.y == 0) && (pad_size.z == 0)) { return arr; }

        auto dims = arr.dims();
        auto new_dims = {dims.x, dims.y + pad_size.y, dims.z + pad_size.z};
        Array<T> new_arr(new_dims);
        new_arr = 0;

        // assume symmetric padding by default
        dims_t d = {pad_size.x, pad_size.y / 2, pad_size.z / 2};
        if (pad_type == PadType::RIGHT) {
            d = {0, 0, 0};
        } else if (pad_type == PadType::LEFT) {
            d = pad_size;
        }

        for (int i = 0; i < dims.x; i++)
            for (int j = 0; j < dims.y; j++)
                for (int k = 0; k < dims.z; k++)
                    new_arr[i, j + d.y, k + d.z] = arr[i, j, k];
        return new_arr;
    }

    template <typename T>
    Array<T> crop2d(const Array<T> &arr, dims_t crop_size, PadType pad_type) {
        if ((crop_size.y == 0) && (crop_size.z == 0)) { return arr; }

        auto dims = arr.dims();
        dims_t new_dims = {dims.x, dims.y - crop_size.y, dims.z - crop_size.z};
        Array<T> new_arr(new_dims);

        // assume symmetric padding by default
        dims_t d = {crop_size.x, crop_size.y / 2, crop_size.z / 2};
        if (pad_type == PadType::RIGHT) {
            d = {0, 0, 0};
        } else if (pad_type == PadType::LEFT) {
            d = crop_size;
        }

        for (int i = 0; i < new_dims.x; i++) {
            for (int j = 0; j < new_dims.y; j++)
                for (int k = 0; k < new_dims.z;)
                    new_arr[i, j, k] = arr[i, j + d.y, j + d.z];
        }
        return new_arr;
    }

    template <typename T>
    Array<T> pad3d(const Array<T> &arr, dims_t pad_size, PadType pad_type) {

        if ((pad_size.x == 0) && (pad_size.y == 0) && (pad_size.z == 0))
            return arr;

        auto dims = arr.dims();
        dims_t new_dims = { dims.x - pad_size.x }
    }

    template <typename T>
    Array<T> crop3d(const Array<T> &arr, dims_t crop_size, PadType pad_type) {

        if ((crop_size.x == 0) && (crop_size.y == 0) && (crop_size.z == 0))
            return arr;

        auto dims = arr.dims();
        dims_t new_dims = {dims.x - crop_size.x, dims.y - crop_size.y,
            dims.z - crop_size.z};
        Array<T> new_arr(new_dims);

        dims_t d = {crop_size.x / 2, crop_size.y / 2, crop_size.z / 2};
        if (pad_type == PadType::LEFT) d = {0, 0, 0};
        else if (pad_type == PadType::RIGHT)
            d = crop_size;

        for (int i = 0; i < new_dims.x; i++) {
            for (int j = 0; j < new_dims.y; j++)
                for (int k = 0; k < new_dims.z; k++)
                    new_arr[i, j, k] = arr[i + d.x, j + d.y, k + d.z];
        }
    }
} // namespace tomocam
#endif // PADDING__H
