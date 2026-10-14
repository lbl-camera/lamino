
#ifndef TOMOCAM_MASK_H
#define TOMOCAM_MASK_H

#include <cmath>
#include <execution>
#include <format>
#include <stdexcept>

#include "array.h"
#include <algorithm>

namespace tomocam {
    /**
     * @brief Creates a mask identifying finite values in an array.
     *
     * Returns a masked array where all the NaNs and Infs are set to 0.
     *
     * @tparam T Numeric type of the projection array elements
     * @param projs Input array to check for finite values
     * @return Array<T> Masked array with NaNs and Infs replaced by 0
     */
    template <typename T>
    Array<T> mask_infs_nans(const Array<T> &projs) {
        Array<T> masked = projs.clone();
        std::for_each(std::execution::par_unseq, masked.begin(), masked.end(),
                      [](T &val) {
                          if (!std::isfinite(val)) {
                              val = 0; // Set NaNs and Infs to 0
                          }
                      });
        return masked;
    }
    /**
     * @brief Creates a mask for laminography/tomography support region.
     *
     * Ptychographic projections are "infinite" in the sense that they
     * add successively add more date to the projection view as the object is
     * rotated. This extra data is missing from projections which are more normal to
     * the beam. This function assumes that the object is contained with rectangular
     * support.
     *
     * @tparam T Numeric type for coordinates and angles
     * @param dims Dimensions of the object being reconstructed (n1, n2, n3)
     * @param sup Dimensions of the support region (n1, n2, n3)
     * @return Array<T> with values outside the support replaced by 0
     */
    template <typename T>
    Array<T> mask_support(const dims_t &dims, const dims_t &sup) {

        const T xcen = static_cast<T>(dims.n3) / 2;
        const T ycen = static_cast<T>(dims.n2) / 2;
        const T zcen = static_cast<T>(dims.n1) / 2;

        Array<T> mask = Array<T>::zeros(dims);
        for (size_t z = 0; z < dims.n1; ++z) {
            for (size_t y = 0; y < dims.n2; ++y) {
                for (size_t x = 0; x < dims.n3; ++x) {
                    T dx = static_cast<T>(x) - xcen;
                    T dy = static_cast<T>(y) - ycen;
                    T dz = static_cast<T>(z) - zcen;

                    if (std::abs(dx) <= static_cast<T>(sup.n3) / 2 &&
                        std::abs(dy) <= static_cast<T>(sup.n2) / 2 &&
                        std::abs(dz) <= static_cast<T>(sup.n1) / 2) {
                        mask[{z, y, x}] = (T)1; // Inside the support region
                    }
                }
            }
        }
        return mask;
    }

    /**
     * @brief Binary support of the reconstructed block, laid out at recon_dims.
     *
     * A file mask of shape (1, n2, n3) (a 2D image) is extruded along n1; a
     * mask of shape (n1, n2, n3) is used as is. Values are binarised with
     * > 0.5. With an empty mask (no file given) every voxel is inside.
     *
     * @param file_mask Mask read from file, or an empty Array
     * @param recon_dims Dimensions of the reconstructed (output) volume
     */
    template <typename T>
    Array<T> support_mask(const Array<T> &file_mask, const dims_t &recon_dims) {
        Array<T> mask = Array<T>::ones(recon_dims);
        if (file_mask.size() == 0) { return mask; }

        const dims_t fd = file_mask.dims();
        if (fd.n2 != recon_dims.n2 || fd.n3 != recon_dims.n3 ||
            (fd.n1 != 1 && fd.n1 != recon_dims.n1)) {
            throw std::invalid_argument(std::format(
                "mask has dimensions [{}, {}, {}], expected [1 or {}, {}, {}]",
                fd.n1, fd.n2, fd.n3, recon_dims.n1, recon_dims.n2, recon_dims.n3));
        }
        const size_t slice = recon_dims.n2 * recon_dims.n3;
        for (size_t k = 0; k < recon_dims.n1; ++k) {
            const T *src = file_mask.begin() + (fd.n1 == 1 ? 0 : k * slice);
            T *dst = mask.begin() + k * slice;
            for (size_t i = 0; i < slice; ++i) {
                dst[i] = src[i] > static_cast<T>(0.5) ? T(1) : T(0);
            }
        }
        return mask;
    }

    /**
     * @brief Place a mask in the centre of a larger volume (the placement of
     * crop3d with PadType::SYMMETRIC); voxels outside it are 0, i.e. free.
     */
    template <typename T>
    Array<T> embed_mask(const Array<T> &mask, const dims_t &full) {
        const dims_t md = mask.dims();
        if (md.n1 > full.n1 || md.n2 > full.n2 || md.n3 > full.n3) {
            throw std::invalid_argument("mask is larger than the volume");
        }
        Array<T> out = Array<T>::zeros(full);
        const size_t o1 = (full.n1 - md.n1) / 2, o2 = (full.n2 - md.n2) / 2,
                     o3 = (full.n3 - md.n3) / 2;
        for (size_t k = 0; k < md.n1; ++k) {
            for (size_t j = 0; j < md.n2; ++j) {
                const T *src = mask.begin() + (k * md.n2 + j) * md.n3;
                T *dst = out.begin() + ((k + o1) * full.n2 + j + o2) * full.n3 + o3;
                std::copy(src, src + md.n3, dst);
            }
        }
        return out;
    }

} // namespace tomocam

#endif // TOMOCAM_MASK_H
