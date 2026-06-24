
#ifndef TOMOCAM_MASK_H
#define TOMOCAM_MASK_H

#include <cmath>
#include <execution>
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

                    if (std::abs(dx) <= sup.n3 / 2 && std::abs(dy) <= sup.n2 / 2 &&
                        std::abs(dz) <= sup.n1 / 2) {
                        mask[{z, y, x}] = (T)1; // Inside the support region
                    }
                }
            }
        }
        return mask;
    }

} // namespace tomocam

#endif // TOMOCAM_MASK_H
