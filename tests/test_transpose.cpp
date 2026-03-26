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

#include <cstdint>
#include <iostream>

#include "array.h"
#include "array_ops.h"
#include "tomocam.h"

using namespace tomocam;

int main() {

    // Test 1: Identity permutation {0, 1, 2} - should return same array
    Array<float> a(2, 3, 4);
    for (size_t i = 0; i < a.size(); i++) { a[i] = static_cast<float>(i); }

    auto b = array::transpose(a, {0, 1, 2});

    bool test1_passed = true;
    if (b.dims() != a.dims()) {
        test1_passed = false;
    } else {
        for (size_t i = 0; i < a.size(); i++) {
            if (a[i] != b[i]) {
                test1_passed = false;
                break;
            }
        }
    }

    std::cout << "test: transpose identity {0,1,2}: .... ";
    if (test1_passed) {
        std::cout << "passed" << std::endl;
    } else {
        std::cout << "failed" << std::endl;
        exit(1);
    }

    // Test 2: Permutation {2, 1, 0} - swap first and last dimensions
    auto d = array::transpose(a, {2, 1, 0});

    bool test2_passed = true;
    dims_t expected_dims = {4, 3, 2};
    if (d.dims() != expected_dims) {
        test2_passed = false;
        std::cout << std::format("transformed dims: ({}, {}, {})\n", d.dims().n1,
                                 d.dims().n2, d.dims().n3);
        std::exit(1);
    } else {
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 3; j++) {
                for (size_t k = 0; k < 2; k++) {
                    if (d[{i, j, k}] != a[{k, j, i}]) {
                        std::cout << std::format(
                            "mismatch at ({}, {}, {}): expected {}, got {}\n", i, j,
                            k, a[{k, j, i}], d[{i, j, k}]);
                        test2_passed = false;
                        break;
                    }
                }
            }
        }
    }

    std::cout << "test: transpose {2,1,0}: .... ";
    if (test2_passed) {
        std::cout << "passed" << std::endl;
    } else {
        std::cout << "failed" << std::endl;
        std::exit(1);
    }

    // Test 3: Permutation {1, 0, 2} - swap first two dimensions
    auto f = array::transpose(a, {1, 0, 2});

    bool test3_passed = true;
    dims_t expected_dims3 = {3, 2, 4};
    if (f.dims() != expected_dims3) {
        test3_passed = false;
    } else {
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 2; j++) {
                for (size_t k = 0; k < 4; k++) {
                    if (f[{i, j, k}] != a[{j, i, k}]) {
                        test3_passed = false;
                        break;
                    }
                }
            }
        }
    }

    std::cout << "test: transpose {1,0,2}: .... ";
    if (test3_passed) {
        std::cout << "passed" << std::endl;
    } else {
        std::cout << "failed" << std::endl;
    }

    // Test 4: Permutation {0, 2, 1} - swap last two dimensions
    auto h = array::transpose(a, {1, 2, 0});

    bool test4_passed = true;
    dims_t expected_dims4 = {3, 4, 2};
    if (h.dims() != expected_dims4) {
        test4_passed = false;
        std::cout << std::format("transformed dims: ({}, {}, {})\n", h.dims().n1,
                                 h.dims().n2, h.dims().n3);
        std::exit(1);
    } else {
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 4; j++) {
                for (size_t k = 0; k < 2; k++) {
                    if (h[{i, j, k}] != a[{k, i, j}]) {
                        std::cout << std::format(
                            "mismatch at ({}, {}, {}): expected {}, got {}\n", i, j,
                            k, a[{k, i, j}], h[{i, j, k}]);
                        test4_passed = false;
                        break;
                    }
                }
            }
        }
    }

    std::cout << "test: transpose {1,2,0}: .... ";
    if (test4_passed) {
        std::cout << "passed" << std::endl;
    } else {
        std::cout << "failed" << std::endl;
    }

    bool all_passed = test1_passed && test2_passed && test3_passed && test4_passed;
    return all_passed ? 0 : 1;
}
