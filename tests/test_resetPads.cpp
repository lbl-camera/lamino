
#include <format>
#include <iostream>
#include <random>

#include "array.h"
#include "dtypes.h"
#include "padding.h"

using namespace tomocam;

template <typename T>
void print_array(const Array<T> &arr) {
    for (size_t i = 0; i < arr.nslices(); ++i) {
        for (size_t j = 0; j < arr.nrows(); ++j) {
            for (size_t k = 0; k < arr.ncols(); ++k) {
                std::cout << std::format(" {}, ", arr[{i, j, k}]);
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }
}

int main() {
    int test_errors = 0;
    dims_t dims = {3, 5, 5}; // Original dimensions
    float pad_factor = 1.0f; // Pad by 100% (double the size)

    auto random = []() {
        static std::mt19937 gen(42); // Fixed seed for reproducibility
        static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        return dist(gen);
    };

    {
        // Test 1: 2D padding/resetPads/cropping
        std::cout << "Test 1: 2D pad, resetPads, crop (3x16x16)\n";
        auto arr = Array<float>::random(dims); // Fill with random values

        // Pad 2D by 0.25 factor (symmetric)
        auto padded = pad2d(arr, pad_factor, PadType::SYMMETRIC);

        std::cout << "  After pad2d: " << padded.dims().n1 << "x" << padded.dims().n2
                  << "x" << padded.dims().n3 << ", pads: " << padded.pads().n1 << ","
                  << padded.pads().n2 << "," << padded.pads().n3 << "\n";

        // fill padded array with random values to ensure resetPads is actually
        // zeroing out
        std::fill(padded.begin(), padded.end(), random());

        std::cout << "  Before resetPads (padded values should be random):\n";
        print_array(padded);
        // Reset pads (zero out padded regions)
        padded.resetPads();
        std::cout << "  After resetPads (padded values should be zero):\n";
        print_array(padded);

        // Crop back to original size
        auto cropped = crop2d(padded, dims, PadType::SYMMETRIC);
        std::cout << "  After crop2d: " << cropped.dims().n1 << "x"
                  << cropped.dims().n2 << "x" << cropped.dims().n3
                  << ", pads: " << cropped.pads().n1 << "," << cropped.pads().n2
                  << "," << cropped.pads().n3 << "\n";

        // Verify dimensions match
        if (cropped.dims().n1 != dims.n1 || cropped.dims().n2 != dims.n2 ||
            cropped.dims().n3 != dims.n3) {
            std::cerr << "  ERROR: Dimensions don't match after crop\n";
            test_errors++;
        } else {
            std::cout << "  ✓ Dimensions preserved\n";
        }
    }
    std::cout << "\n==============\n";

    {
        std::cout << "Test 2:  3D pad, resetPads, crop (3x5x5)\n";
        auto arr = Array<float>::random(dims);

        // Set up pads: pad all dimensions by 2
        auto padded = pad3d(arr, pad_factor, PadType::SYMMETRIC);
        // Fill padded array with random values to ensure resetPads is actually
        // zeroing out
        std::fill(padded.begin(), padded.end(), random());
        std::cout << "  Before resetPads (padded values should be random):\n";
        print_array(padded);
        std::cout << "\n-----------------\n";
        padded.resetPads();
        std::cout << "  After resetPads (padded values should be 0, unpadded should "
                     "be 1):\n";
        print_array(padded);
        auto cropped = crop3d(padded, dims, PadType::SYMMETRIC);
        // print dims and pads of cropped
        std::cout << std::format("  Cropped dims: {}x{}x{}, pads: {},{},{}\n",
                                 cropped.dims().n1, cropped.dims().n2,
                                 cropped.dims().n3, cropped.pads().n1,
                                 cropped.pads().n2, cropped.pads().n3);
    }
}
