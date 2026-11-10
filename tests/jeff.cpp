#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "array_ops.h"
#include "tiff.h"
#include "tomocam.h"

using namespace tomocam;
int main(int argc, char *argv[]) {

    // load phantom data (121, 121, 17) from text file
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <phantom_file.txt>" << std::endl;
        return 1;
    }
    std::string phantom_path = argv[1];
    std::filesystem::path phantom_file(phantom_path);
    if (!std::filesystem::exists(phantom_file)) {
        std::cerr << "Phantom file does not exist: " << phantom_file << std::endl;
        return 1;
    }
    // read text file into a vector of floats
    std::vector<float> xval, yval, zval;
    std::ifstream infile(phantom_file);
    if (!infile.is_open()) {
        std::cerr << "Failed to open phantom file: " << phantom_file << std::endl;
        return 1;
    }
    // each lines has 3 float values, read them into separate vectors
    float v1, v2, v3;
    while (infile >> v1 >> v2 >> v3) {
        xval.push_back(v1);
        yval.push_back(v2);
        zval.push_back(v3);
    }

    // [x, y, z]
    tomocam::dims_t dims = {121, 121, 17};
    std::array<Array<float>, 3> phantom_data;
    for (int i = 0; i < 3; ++i) { phantom_data[i] = std::move(Array<float>(dims)); }
    std::copy(xval.begin(), xval.end(), phantom_data[0].begin());
    std::copy(yval.begin(), yval.end(), phantom_data[1].begin());
    std::copy(zval.begin(), zval.end(), phantom_data[2].begin());

    // write vti file for visualization in paraview
    tomocam::vti::write_vectors("phantom.vti", phantom_data);

    // transpose the data to (17, 121, 121)
    for (int i = 0; i < 3; ++i) {
        phantom_data[i] = array::transpose(phantom_data[i], {1, 2, 0});
    }

    // create a vector of projection angles -70 to 70 degrees in 1 degree increments
    std::vector<float> angles;
    for (int i = -70; i <= 70; ++i) { angles.push_back(i * M_PI / 180.0f); }

    size_t nrows = 2 * dims.n1 - 1;
    size_t ncols = 2 * dims.n2 - 1;

    // output directory
    std::filesystem::path home_dir = std::filesystem::path(getenv("HOME"));
    std::filesystem::path output_dir = home_dir / "data/jeff";
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }
    std::vector<float> orientations = {-30, 0, 30};
    // for (int i = -90; i <= 90; i += 15) { orientations.push_back(i); }
    std::cout << "Output directory: " << output_dir << std::endl;
    // do the forward projection
    for (int i = 0; i < orientations.size(); ++i) {
        float orientation = orientations[i];

        // pad the phantom data with zeros to  2x the size in each dimension
        for (int j = 0; j < 3; ++j) {
            phantom_data[j] =
                pad3d<float>(phantom_data[j], 1, tomocam::PadType::SYMMETRIC);
        }

        // create polar grid
        float angle_rad = orientation * M_PI / 180.0f;
        tomocam::PolarGrid polar_grid(angles, nrows, ncols, angle_rad);

        // create forward projection object
        auto proj = forward(phantom_data, polar_grid, orientation);

        // cropy projection to original size (121, 121)
        dims_t proj_dims = {angles.size(), dims.n1, dims.n2};
        proj = crop2d<float>(proj, proj_dims, tomocam::PadType::SYMMETRIC);

        // save the projection to tiff file
        std::string filename = "gamma" + std::to_string((int)orientation) + ".tiff";
        std::filesystem::path output_file = output_dir / filename;
        /*
        tomocam::tiff::write(output_file, proj);
        */
        std::cout << "Saved projection to: " << output_file << std::endl;
    }
    return 0;
}
