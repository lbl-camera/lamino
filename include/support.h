#ifndef SUPPORT_H
#define SUPPORT_H

#include <algorithm>
#include <execution>

#include "array.h"

namespace tomocam {
    template <typename T>
    class Support {
      private:
        Array<T> support_;

      public:
        Support(dims_t dims, dims_t supp_dims) : support_(dims) {

            int half_n1 = (int)(supp_dims.n1 / 2);
            int half_n2 = (int)(supp_dims.n2 / 2);
            int half_n3 = (int)(supp_dims.n3 / 2);

            int center_i = (int)dims.n1 / 2;
            int center_j = (int)dims.n2 / 2;
            int center_k = (int)dims.n3 / 2;

            for (int i = 0; i < (int)dims.n1; i++) {
                for (int j = 0; j < (int)dims.n2; j++) {
                    for (int k = 0; k < (int)dims.n3; k++) {
                        int x = i - center_i;
                        int y = j - center_j;
                        int z = k - center_k;
                        if ((x + half_n1) >= 0 &&
                            (x + half_n1) < (int)supp_dims.n1 &&
                            (y + half_n2) >= 0 &&
                            (y + half_n2) < (int)supp_dims.n2 &&
                            (z + half_n3) >= 0 &&
                            (z + half_n3) < (int)supp_dims.n3) {
                            support_[{i, j, k}] = 1;
                        } else {
                            support_[{i, j, k}] = 0;
                        }
                    }
                }
            }
        }

        void apply(Array<T> &data) const {
            std::transform(std::execution::par_unseq, support_.begin(),
                           support_.end(), data.begin(), data.begin(),
                           std::multiplies<T>());
        }
    };
} // namespace tomocam
#endif // SUPPORT_H
