#include <cstdint>
#include <fstream>
#include <toml++/toml.h>
#include <vector>

#include "shepplogan.h"
#include "tiff.h"
#include "tomocam.h"

float bilinear(float x, float y, const tomocam::Slice<float> &img) {
    size_t i0 = (size_t)y;
    if ((i0 < 0) || (i0 >= img.dims[0])) return 0;
    size_t j0 = (size_t)x;
    if ((j0 < 0) || (j0 >= img.dims[1])) return 0;
    float dx = x - j0;
    float dy = y - i0;
    size_t i1 = i0 + 1;
    size_t j1 = j0 + 1;
    float val = 0;
    if ((i1 >= img.dims[0]) && (j1 >= img.dims[1]))
        val = img[{i0, j0}];
    else if ((i1 < img.dims[0]) && (j1 >= img.dims[1]))
        val = img[{i0, j0}] * dy + img[{i1, j0}] * (1 - dy);
    else if ((i1 >= img.dims[0]) && (j1 < img.dims[1]))
        val = img[{i0, j0}] * dx + img[{i0, j1}] * (1 - dx);
    else
        val = (img[{i0, j0}] * dy * dx + img[{i0, j1}] * dy * (1 - dx) +
               img[{i1, j0}] * (1 - dy) * dx + img[{i1, j1}] * (1 - dy) * (1 - dx));
    return val;
}

int main(int argc, char **argv) {

    if (argc < 1) {
        std::cerr << std::format("Usage: {} <config.toml>", argv[0]) << std::endl;
        return 1;
    }
    std::fstream toml_file(argv[1]);
    if (!toml_file.is_open()) {
        std::cerr << std::format("Error: Could not open TOML file '{}'", argv[1])
                  << std::endl;
        return 1;
    }
    auto params = toml::parse(toml_file);
    size_t nslcs = params["nslices"].value_or(21);
    size_t nrows = params["nrows"].value_or(511);
    size_t ncols = params["ncols"].value_or(511);
    float rot_ang = params["rot_angle"].value_or(0.0f);
    std::string outfile = params["output"].value_or("shepplogan.tiff");

    tomocam::Array<float> data(nslcs, nrows, ncols);

    auto sheppSlice = sheppLoganPhantom(nrows, ncols);
    for (size_t i = 0; i < nslcs; i++) {
        auto slc = data.slice(i);
        std::copy(sheppSlice.begin(), sheppSlice.end(), slc.ptr);
    }
    tomocam::tiff::write(outfile, data);
    return 0;
}
