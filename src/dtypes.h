#include <memory>
#include <stdexcept>

#ifndef DTYPES__H
#define DTYPES__H

namespace tomocam {

    struct dims_t {
        size_t x;
        size_t y;
        size_t z;
        size_t size() const { return x * y * z; }

        dims_t operator+(const dims_t &v) {
            return {this->x + v.x, this->y + v.y, this->z + v.z};
        }

        dims_t operator-(const dims_t &v) {
            return {this->x - v.x, this->y - v.y, this->z - v.z};
        }

        dims_t operator*(int v) {
            return {this->x * v, this->y * v, this->z * v};
        }

        dims_t operator/(int v) {
            if (v == 0) std::runtime_error("divide by zeros");
            return {this->x / v, this->y / v, this->z / v};
        }
    };

} // namespace tomocam
#endif // DTYPES__H
