#include "array.h"
#include "array_ops.h"
#include "optimize.h"

// Power iteration to estimate the Lipschitz constant of A
namespace tomocam::opt {
    template <typename T>
    T lipschitz(const Function<T> &A, const Array<T> &ref, size_t max_iter, T tol) {
        Array<T> x = Array<T>::random_like(ref);
        x = x / array::norm2(x);
        T L = 0;
        for (int i = 0; i < max_iter; ++i) {
            auto Ax = A(x);
            T norm_Ax = array::norm2(Ax);
            // avoid division by almost zero
            if (std::abs(norm_Ax) < 1e-12) {
                L = 0;
                break;
            }
            x = Ax / norm_Ax;
            T L_new = norm_Ax;
            if (std::abs(L_new - L) < tol) { break; }
            L = L_new;
        }
        return L;
    }
    // explicit template instantiation
    template float lipschitz(const Function<float> &, const Array<float> &, size_t,
                             float);
    template double lipschitz(const Function<double> &, const Array<double> &,
                              size_t, double);
} // namespace tomocam::opt
