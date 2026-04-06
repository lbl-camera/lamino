
#include <algorithm>
#include <cmath>
#include <execution>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "array.h"
#include "mask.h"
#include "tiff.h"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.tif> <angles.txt>\n", argv[0]);
        return 1;
    }

    // Load the input image
    auto projs = tomocam::tiff::read(argv[1]);

    // --- Test mask_infs_nans ---
    std::cout << "--- Testing mask_infs_nans ---\n";

    // check if image has nans
    size_t num_nans = 0;
    std::for_each(std::execution::par_unseq, projs.begin(), projs.end(),
                  [&](float val) {
                      if (std::isnan(val)) { num_nans++; }
                  });
    std::cout << std::format("Image has {} NaNs\n", num_nans);

    // create a mask of the nans and infs
    auto m0 = tomocam::mask_infs_nans(projs);

    // apply the mask to the projections
    size_t masked_nans = 0;
    auto maked_projs = projs.clone();
    for (size_t i = 0; i < projs.size(); i++) {
        if (m0[i] == 1) {
            if (std::isnan(maked_projs[i])) { masked_nans++; }
        }
    }
    std::cout << std::format("Masked projections have {} NaNs\n", masked_nans);

    // --- Test mask_support ---
    std::cout << "\n--- Testing mask_support ---\n";
    size_t nprojs = projs.nslices();
    size_t nrows = projs.nrows();
    size_t ncols = projs.ncols();

    // adjust dims so that (digonal of [rows, cols] fits within the camera FOV at all
    // angles)
    double diag = std::sqrt(static_cast<double>(nrows * nrows + ncols * ncols));
    double scale = static_cast<double>(nrows) / diag;
    size_t new_n1 = static_cast<size_t>(std::ceil(nrows * scale));
    size_t new_n2 = static_cast<size_t>(std::ceil(ncols * scale));
    tomocam::dims_t dims{nprojs, new_n1, new_n2};
    std::cout << std::format("Masking dims: nprojs={}, rows={}, cols={}\n", nprojs,
                             new_n1, new_n2);

    // read angles from the text file
    std::string angle_file = argv[2];
    std::vector<float> theta;
    std::ifstream infile(angle_file);
    if (!infile) {
        std::cerr << "Error opening angle file: " << angle_file << std::endl;
        return 1;
    }
    float angle;
    while (infile >> angle) { theta.push_back(angle); }
    if (theta.size() != nprojs) {
        std::cerr << "Error: number of angles (" << theta.size()
                  << ") does not match number of projections (" << nprojs << ")\n";
        return 1;
    }
    // convert angles from degrees to radians
    std::transform(theta.begin(), theta.end(), theta.begin(),
                   [](float deg) { return deg * M_PI / 180.0f; });

    // Laminography tilt angle: -45 degrees
    const float gamma = -45.0f * M_PI / 180.0f;

    auto support_mask = tomocam::mask_support(projs, dims, theta, gamma);
    // convert mask to float for visualization
    tomocam::Array<float> support(support_mask.dims());
    std::transform(support_mask.begin(), support_mask.end(), support.begin(),
                   [](int m) { return m == 1 ? 1.0f : 0.0f; });

    // save mask as tiff for visualization
    std::string mask_output = "support_mask.tif";
    tomocam::tiff::write(mask_output, support);
    std::cout << std::format("Support mask saved to {}\n", mask_output);

    // --- Test combine_masks ---
    return 0;
}
