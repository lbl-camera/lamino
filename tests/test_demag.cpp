/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals from the
 * U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 * IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 * the U.S. Government has been granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 * to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#include <array>
#include <cmath>
#include <iostream>
#include <random>

#include "array.h"
#include "demag.h"
#include "fft.h"

using namespace tomocam;

template <typename T>
bool check_for_invalid_values(const std::array<Array<T>, 3> &data, const std::string &name) {
    bool has_invalid = false;
    
    for (int comp = 0; comp < 3; ++comp) {
        for (size_t i = 0; i < data[comp].size(); ++i) {
            if (std::isnan(data[comp][i]) || std::isinf(data[comp][i])) {
                std::cout << "  ERROR: " << name << "[" << comp << "][" << i << "] = " 
                          << data[comp][i] << " (invalid)" << std::endl;
                has_invalid = true;
                break;
            }
        }
        if (has_invalid) break;
    }
    
    return has_invalid;
}

template <typename T>
bool check_for_invalid_complex(const Array<std::complex<T>> &data, const std::string &name) {
    for (size_t i = 0; i < data.size(); ++i) {
        T real_part = data[i].real();
        T imag_part = data[i].imag();
        if (std::isnan(real_part) || std::isinf(real_part) || 
            std::isnan(imag_part) || std::isinf(imag_part)) {
            std::cout << "  ERROR: " << name << "[" << i << "] = " 
                      << data[i] << " (invalid)" << std::endl;
            return true;
        }
    }
    return false;
}

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist_float(0.0f, 1.0f);
    std::uniform_real_distribution<double> dist_double(0.0, 1.0);

    // Test FFT round-trip with float
    {
        std::cout << "test: FFT round-trip <float>: ";
        dims_t dims = {16, 16, 16};
        Array<float> m(dims);
        
        for (size_t i = 0; i < m.size(); ++i) {
            m[i] = dist_float(gen);
        }
        
        auto m_hat = fft::fft3_r2c(m);
        bool has_invalid = check_for_invalid_complex(m_hat, "m_hat_float");
        
        if (has_invalid) {
            std::cout << "failed in forward FFT" << std::endl;
        } else {
            auto m_back = fft::fft3_c2r(m_hat, dims);
            has_invalid = false;
            for (size_t i = 0; i < m_back.size(); ++i) {
                if (std::isnan(m_back[i]) || std::isinf(m_back[i])) {
                    std::cout << "  ERROR: m_back[" << i << "] = " << m_back[i] << std::endl;
                    has_invalid = true;
                    break;
                }
            }
            if (has_invalid) {
                std::cout << "failed in inverse FFT" << std::endl;
            } else {
                std::cout << "passed" << std::endl;
            }
        }
    }

    // Test with float
    {
        std::cout << "test: demag<float> with random input [0,1]: ";
        
        dims_t dims = {16, 16, 16};
        std::array<Array<float>, 3> m;
        
        for (int i = 0; i < 3; ++i) {
            m[i] = Array<float>(dims);
            for (size_t j = 0; j < m[i].size(); ++j) {
                m[i][j] = dist_float(gen);
            }
        }
        
        auto H = opt::demag(m);
        
        std::cout << std::endl;
        std::cout << "  H[0][0] = " << H[0][0] << std::endl;
        std::cout << "  H[1][0] = " << H[1][0] << std::endl;
        std::cout << "  H[2][0] = " << H[2][0] << std::endl;
        
        bool has_invalid = check_for_invalid_values(H, "H_float");
        
        if (has_invalid) {
            std::cout << "failed (NaN or Inf detected)" << std::endl;
        } else {
            std::cout << "passed (no NaN or Inf)" << std::endl;
        }
    }

    // Test with double
    {
        std::cout << "test: demag<double> with random input [0,1]: ";
        
        dims_t dims = {16, 16, 16};
        std::array<Array<double>, 3> m;
        
        for (int i = 0; i < 3; ++i) {
            m[i] = Array<double>(dims);
            for (size_t j = 0; j < m[i].size(); ++j) {
                m[i][j] = dist_double(gen);
            }
        }
        
        auto H = opt::demag(m);
        
        bool has_invalid = check_for_invalid_values(H, "H_double");
        
        if (has_invalid) {
            std::cout << "failed (NaN or Inf detected)" << std::endl;
        } else {
            std::cout << "passed (no NaN or Inf)" << std::endl;
        }
    }

    return 0;
}
