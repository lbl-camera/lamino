
#include <chrono>
#include <iostream>
#include <random>

#ifndef UTILS__H
#define UTILS__H

namespace cam {
    class NumPyRandom {
      private:
        std::mt19937 gen_;

      public:
        NumPyRandom() { gen_ = std::mt19937(5489u); }

        template <typename T>
        T rand() {
            int a = gen_() >> 5;
            int b = gen_() >> 6;
            return static_cast<T>((a * 67108864.0 + b) / 9007199254740992.0);
        }
    };

} // namespace cam
#endif // UTILS__H
