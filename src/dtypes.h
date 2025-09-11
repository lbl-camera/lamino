#include <cstdint>
#include <stdexcept>
#include <tuple>

#ifndef DTYPES__H
#define DTYPES__H

namespace tomocam {

    class dims_t {
        uint64_t n1;
        uint64_t n2;
        uint64_t n3;

      public:
        explicit dims_t(uint64_t a, uint64_t b, uint64_t c) : n1(a), n2(b), n3(c) {}

        std::tuple<uint64_t, uint64_t, uint64_t> unravel_idx(uint64_t idx) const {
            uint64_t i = idx / n2 / n3;
            uint64_t j = (idx / n3) % n2;
            uint64_t k = idx % n3;
            return std::make_tuple(i, j, k);
        }

        uint64_t size() const { return (n1 * n2 * n3); }

        uint64_t flat_idx(uint64_t i, uint64_t j, uint64_t k) const {
            return ((i * n2 + j) * n3 + k);
        }
        uint64_t x() const { return n1; }
        uint64_t y() const { return n2; }
        uint64_t z() const { return n3; }

        dims_t operator+(const dims_t &v) { return dims_t(n1 + v.n1, n2 + v.n2, n3 + v.n3); }
        dims_t operator-(const dims_t &v) { return dims_t(n1 - v.n1, n2 - v.n2, n3 - v.n3); }
        dims_t operator*(int v) { return dims_t(n1 * v, n2 * v, n3 * v); }
        dims_t operator/(int v) {
            if (v == 0)
                throw std::runtime_error("divide by zeros");
            return dims_t(n1 / v, n2 / v, n3 / v);
        }
    };
} // namespace tomocam
#endif // DTYPES__H
