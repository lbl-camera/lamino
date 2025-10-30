#include <execution>
#include <format>
#include <functional>
#include <iostream>

#include "array.h"
#include "array_ops.h"

namespace tomocam::opt {

    template <typename T>
    Array<T> cgsolver(std::function<Array<T>(const Array<T> &)> A,
        const Array<T> &yT, size_t max_iter, T tol) {

        // initialize
        auto x = yT.clone();
        auto r = yT - A(x);
        auto p = r.clone();
        auto rs_old = array::dot(r, r);

        for (size_t iter = 0; iter < max_iter; iter++) {

            auto Ap = A(p);
            auto pAp = array::dot(p, Ap);
            if (std::abs(pAp) < 1.e-15) {
                std::cerr << "pAp is close to zero\n";
                break;
            }
            auto alpha = rs_old / pAp;
            x += p * alpha;
            r -= Ap * alpha;

            auto rs_new = array::norm2(r);
            std::cout << std::format("iter: {}, residual: {}\n", iter, rs_new);
            if (rs_new < tol) { break; }

            p += r + (p * (rs_new / rs_old));
            rs_old = rs_new;
        }
        return x;
    }

    // template instantiations
    template Array<float> cgsolver<float>(
        std::function<Array<float>(const Array<float> &)>, const Array<float> &,
        size_t, float);
    template Array<double> cgsolver<double>(
        std::function<Array<double>(const Array<double> &)>,
        const Array<double> &, size_t, double);

} // namespace tomocam::opt
