
#include <cmath>
#include <tr1/cmath>

#ifndef KAISER_WINDOW__H
#define KAISER_WINDOW__H

class Kernel {
 private:
    float radius_;
    float beta_;

 public:
     Kernel(float r, float b):radius_(r), beta_(b) {
    } float radius() {
        return radius_;
    }
    float beta() {
        return beta_;
    }
    float window_size() {
        return 2 * radius_ + 1;
    }

    float weight(float r) {
        float weight = 0.f;
        float W = window_size();
        if (std::abs(r) < radius_) {
            float t0 = 2 * r / W;
            float t1 =
                std::tr1::cyl_bessel_i(0, beta_ * std::sqrt(1 - t0 * t0));
            float t2 = std::tr1::cyl_bessel_i(0, beta_);
            weight = t1 / t2;
        }
        return weight;
    }

    float kaiserFT(float k) {
        float W = window_size();
		float w0 = std::sinh(beta_) / beta_;
        float t1 = std::sqrt(std::pow(beta_, 2) - std::pow(k * M_PI * W, 2));
		return  std::sinh(t1) / t1 / w0;
    }
};

#endif // KAISER_WINDOW__H
