#include <cstdint>
#include <fstream>
#include <vector>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "tiff.h"
#include "tomocam.h"

struct Image {
    int n1;
    int n2;
    std::vector<float> pxl;
    Image(uint64_t n1_, uint64_t n2_) :
        n1(n1_), n2(n2_), pxl(std::vector<float>(n1 * n2, 0)) {}
    float &operator[](int i, int j) { return pxl[i * n2 + j]; }
    const float &operator[](int i, int j) const { return pxl[i * n2 + j]; }
};

float bilinear(float x, float y, const Image &img) {
    int i0 = (int)y;
    if ((i0 < 0) || (i0 >= img.n1)) return 0;
    int j0 = (int)x;
    if ((j0 < 0) || (j0 >= img.n2)) return 0;
    float dx = x - j0;
    float dy = y - i0;
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    float val = 0;
    if ((i1 >= img.n1) && (j1 >= img.n2)) val = img[i0, j0];
    else if ((i1 < img.n1) && (j1 >= img.n2))
        val = img[i0, j0] * dy + img[i1, j0] * (1 - dy);
    else if ((i1 >= img.n1) && (j1 < img.n2))
        val = img[i0, j0] * dx + img[i0, j1] * (1 - dx);
    else
        val = (img[i0, j0] * dy * dx + img[i0, j1] * dy * (1 - dx) +
               img[i1, j0] * (1 - dy) * dx + img[i1, j1] * (1 - dy) * (1 - dx));
    return val;
}

Image rotate(const Image &img, float angle) {

    // copy for allocation
    auto newimg = Image(img.n1, img.n2);
    float ycen = img.n1 / 2.f;
    float xcen = img.n2 / 2.f;
    auto cost = std::cos(angle);
    auto sint = std::sin(angle);
    for (int y = 0; y < img.n1; y++) {
        for (int x = 0; x < img.n2; x++) {
            float dy = y - ycen;
            float dx = x - xcen;
            float dxr = cost * dx - sint * dy;
            float dyr = sint * dx + cost * dy;
            newimg[y, x] = bilinear(dxr + xcen, dyr + ycen, img);
        }
    }
    return newimg;
}

int main(int argc, char **argv) {

    int nslcs = 511;
    int nrows = 511;
    int ncols = 511;
    float rot_ang = 0.f;
    std::string outfile("chess.tif");

    if (argc > 1) {
        std::fstream json_file(argv[1]);
        if (json_file.is_open()) {
            auto params = json::parse(json_file);
            nslcs = params["nslices"];
            nrows = params["nrows"];
            ncols = params["ncols"];
            rot_ang = params["rot_angle"];
            outfile = params["output_filename"];
        }
    }

    tomocam::Array<float> data(nslcs, nrows, ncols);
    data.fill(0.f);

    int ibeg = nslcs / 2 - 10;
    int iend = nslcs / 2 + 11;
    int nbox = 32;
    for (int i = ibeg; i < iend; i++) {
        auto img = Image(nrows, ncols);
        for (int j = 2 * nbox; j < nrows - 2 * nbox; j++) {
            for (int k = 3 * nbox; k < ncols - 3 * nbox; k++) {
                int m = j / nbox;
                int n = k / nbox;
                if ((m + n) & 1) {
                    img[j, k] = 2.f;

                } else {
                    img[j, k] = 1.f;
                }
            }
        }
        img = rotate(img, rot_ang);
        std::copy(img.pxl.begin(), img.pxl.end(), data.slice(i));
    }
    tomocam::tiff::write(outfile, data);
    return 0;
}
