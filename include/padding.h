/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals from the
 * U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 * IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 * the U.S. Government has been granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 * to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#ifndef PADDING__H
#define PADDING__H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <system_error>

#include "array.h"
#include "dtypes.h"

namespace tomocam {
    enum class PadType { LEFT, RIGHT, SYMMETRIC };

    template <typename T>
    size_t n_pad(size_t n, T factor) {
        auto n2 = static_cast<size_t>(n * factor);
        return 2 * ((n2 - n) / 2);
    }

    template <typename T>
    Array<T> pad2d(const Array<T> &arr, T factor, PadType pad_type) {

        // if factor is zero, return copy of input array
        if (factor - 1 < 1.e-06) { return arr.clone(); }

        size_t n2 = static_cast<size_t>(arr.nrows() * factor);
        size_t n3 = static_cast<size_t>(arr.ncols() * factor);

        // create and initialize return array
        dims_t dims{arr.nslices(), n2, n3};
        Array<T> arr2(dims);

        // assume symmetric padding by default
        size_t d1, d2;
        if (pad_type == PadType::RIGHT) {
            d1 = d2 = 0;
        } else if (pad_type == PadType::LEFT) {
            d1 = n2 - arr.nrows();
            d2 = n3 - arr.ncols();
        } else {
            d1 = (n2 - arr.nrows()) / 2;
            d2 = (n3 - arr.ncols()) / 2;
        }

        for (size_t i = 0; i < arr.nslices(); i++) {
            Slice<T> in = arr.slice(i, i + 1);
            Slice<T> out = arr2.slice(i, i + 1);
            for (size_t j = 0; j < arr.nrows(); j++) {
                std::copy(in.begin() + j * arr.ncols(),
                          in.begin() + j * arr.ncols() + arr.ncols(),
                          out.begin() + (j + d1) * arr2.ncols() + d2);
            }
        }
        return arr2;
    }

    template <typename T>
    Array<T> crop2d(const Array<T> &arr, dims_t new_dims, PadType pad_type) {

        auto crop_size = arr.dims() - new_dims;
        if ((crop_size.n2 == 0) && (crop_size.n3 == 0)) { return arr.clone(); }

        Array<T> arr2(new_dims);

        // assume symmetric padding by default
        size_t d1, d2;
        if (pad_type == PadType::RIGHT) {
            d1 = d2 = 0;
        } else if (pad_type == PadType::LEFT) {
            d1 = crop_size.n2;
            d2 = crop_size.n3;
        } else {
            d1 = crop_size.n2 / 2;
            d2 = crop_size.n3 / 2;
        }

        for (size_t i = 0; i < new_dims.n1; i++) {
            Slice<T> in = arr.slice(i, i + 1);
            Slice<T> out = arr2.slice(i, i + 1);
            for (size_t j = 0; j < new_dims.n2; j++) {
                std::copy(in.begin() + (j + d1) * arr.ncols() + d2,
                          in.begin() + (j + d1) * arr.ncols() + d2 + new_dims.n3,
                          out.begin() + j * arr2.ncols());
            }
        }
        return arr2;
    }

    template <typename T>
    auto pad3d(const Array<T> &arr, T factor, PadType pad_type) -> Array<T> {

        if (factor - 1 < 1.e-06) { return arr.clone(); }

        // calculate new dims
        auto n1 = static_cast<size_t>(arr.nslices() * factor);
        auto n2 = static_cast<size_t>(arr.nrows() * factor);
        auto n3 = static_cast<size_t>(arr.ncols() * factor);

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
        for (size_t i = 0; i < arr.nslices(); i++) {
            Slice<T> in = arr.slice(i, i + 1);
            Slice<T> out = arr2.slice(d.n1 + i, d.n1 + i + 1);
            for (size_t j = 0; j < arr.nrows(); j++) {
                std::copy(in.begin() + j * arr.ncols(),
                          in.begin() + j * arr.ncols() + arr.ncols(),
                          out.begin() + (j + d.n2) * arr2.ncols() + d.n3);
            }
        }
        return arr2;
    }

    template <typename T>
    Array<T> crop3d(const Array<T> &arr, dims_t new_dims, PadType pad_type) {

        auto crop_size = arr.dims() - new_dims;
        if ((crop_size.n1 == 0) && (crop_size.n2 == 0) && (crop_size.n3 == 0))
            return arr.clone();
        Array<T> arr2(new_dims);

        dims_t d = crop_size / 2;
        if (pad_type == PadType::RIGHT)
            d = dims_t(0, 0, 0);
        else if (pad_type == PadType::LEFT)
            d = crop_size;

        for (size_t i = 0; i < new_dims.n1; i++) {
            Slice<T> in = arr.slice(i + d.n1, i + d.n1 + 1);
            Slice<T> out = arr2.slice(i, i + 1);
            for (size_t j = 0; j < new_dims.n2; j++) {
                std::copy(in.begin() + (j + d.n2) * arr.ncols() + d.n3,
                          in.begin() + (j + d.n2) * arr.ncols() + d.n3 + new_dims.n3,
                          out.begin() + j * arr2.ncols());
            }
        }
        return arr2;
    }
} // namespace tomocam
#endif // PADDING__H
