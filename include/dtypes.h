#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <tuple>

#ifndef DTYPES__H
    #define DTYPES__H

namespace tomocam {

    struct dims_t {
        size_t n1;
        size_t n2;
        size_t n3;

      public:
        dims_t() : n1(0), n2(0), n3() {}
        dims_t(size_t a, size_t b, size_t c) : n1(a), n2(b), n3(c) {}

        [[nodiscard]] std::tuple<size_t, size_t, size_t> unravel_idx(
            size_t idx) const {
            size_t i = idx / n2 / n3;
            size_t j = (idx / n3) % n2;
            size_t k = idx % n3;
            return std::make_tuple(i, j, k);
        }

        // change individual dims (needed in padding)
        void set_x(size_t v) { n1 = v; }
        void set_y(size_t v) { n2 = v; }
        void set_z(size_t v) { n3 = v; }

        [[nodiscard]] size_t size() const { return (n1 * n2 * n3); }

        [[nodiscard]] size_t flat_idx(size_t i, size_t j, size_t k) const {
            return (((i * n2 + j) * n3) + k);
        }
        [[nodiscard]] size_t x() const { return n1; }
        [[nodiscard]] size_t y() const { return n2; }
        [[nodiscard]] size_t z() const { return n3; }

        dims_t operator+(const dims_t &v) const {
            return {n1 + v.n1, n2 + v.n2, n3 + v.n3};
        }
        dims_t operator-(const dims_t &v) const {
            return {n1 - v.n1, n2 - v.n2, n3 - v.n3};
        }
        dims_t operator*(int v) const { return {n1 * v, n2 * v, n3 * v}; }
        dims_t operator/(int v) const {
            if (v == 0) { throw std::runtime_error("divide by zeros"); }
            return {n1 / v, n2 / v, n3 / v};
        }
    };

    inline std::ostream &operator<<(std::ostream &outs, dims_t d) {
        outs << "dims_t(" << d.x() << ", " << d.y() << ", " << d.z() << ")";
        return outs;
    }
} // namespace tomocam
#endif // DTYPES__H
