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
#ifndef TOMOCAM_TIFF__H
#define TOMOCAM_TIFF__H

#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <tiffio.h>
#include <vector>

#include "array.h"

namespace tomocam::tiff {

    namespace detail {
        struct TiffCloser {
            void operator()(TIFF *t) const { TIFFClose(t); }
        };
        using TiffPtr = std::unique_ptr<TIFF, TiffCloser>;

        inline TiffPtr open(const std::string &filename, const char *mode) {
            TiffPtr tif(TIFFOpen(filename.c_str(), mode));
            if (!tif) throw std::runtime_error("tiff: cannot open file " + filename);
            return tif;
        }
    } // namespace detail

    inline Array<float> read(std::string filename) {

        // count number of projections
        uint32_t npages = 0;
        {
            auto counter = detail::open(filename, "r");
            do { npages++; } while (TIFFReadDirectory(counter.get()));
        }

        // reopen file for reading
        auto tif = detail::open(filename, "r");
        TIFF *tif_ = tif.get();

        // get image size
        uint16_t bits, format;
        uint32_t w, h;
        TIFFGetField(tif_, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tif_, TIFFTAG_IMAGELENGTH, &h);
        TIFFGetField(tif_, TIFFTAG_BITSPERSAMPLE, &bits);
        TIFFGetFieldDefaulted(tif_, TIFFTAG_SAMPLEFORMAT, &format);

        // allocate memory
        size_t nscls = static_cast<size_t>(npages);
        size_t nrows = static_cast<size_t>(h);
        if (w % 2 == 0) {
            nrows -= 1; // ensure odd number of rows for symmetry around center
        }
        size_t ncols = static_cast<size_t>(w);
        if (h % 2 == 0) {
            ncols -= 1; // ensure odd number of cols for symmetry around center
        }
        Array<float> data(nscls, nrows, ncols);

        // supported pixel types: float32, float64, uint32
        enum class Kind { F32, F64, U32 };
        Kind kind;
        if (format == SAMPLEFORMAT_IEEEFP && bits == 32) {
            kind = Kind::F32;
        } else if (format == SAMPLEFORMAT_IEEEFP && bits == 64) {
            kind = Kind::F64;
        } else if (format == SAMPLEFORMAT_UINT && bits == 32) {
            kind = Kind::U32;
        } else {
            throw std::runtime_error(
                "tiff: unsupported format (bits = " + std::to_string(bits) +
                ", sample format = " + std::to_string(format) + ") in " + filename);
        }

        tsize_t line_size = TIFFScanlineSize(tif_);
        if (line_size != static_cast<tsize_t>(w * (bits / 8))) {
            throw std::runtime_error("tiff: line_size, width mismatch in " +
                                     filename);
        }
        std::vector<char> raw(line_size);
        void *buf = raw.data();

        for (size_t i = 0; i < nscls; i++) {
            TIFFSetDirectory(tif_, static_cast<tdir_t>(i));
            for (size_t j = 0; j < nrows; j++) {
                if (TIFFReadScanline(tif_, buf, static_cast<uint32_t>(j)) < 0) {
                    throw std::runtime_error("tiff: failed to read scanline " +
                                             std::to_string(j) + " of page " +
                                             std::to_string(i) + " in " + filename);
                }
                for (size_t k = 0; k < ncols; k++) {
                    switch (kind) {
                        case Kind::F32:
                            data[{i, j, k}] = static_cast<float *>(buf)[k];
                            break;
                        case Kind::F64:
                            data[{i, j, k}] =
                                static_cast<float>(static_cast<double *>(buf)[k]);
                            break;
                        case Kind::U32:
                            data[{i, j, k}] =
                                static_cast<float>(static_cast<uint32_t *>(buf)[k]);
                            break;
                    }
                }
            }
        }
        return data;
    }

    inline void write(std::string filename, const Array<float> &data) {

        // open file
        auto tif = detail::open(filename, "w");
        TIFF *tif_ = tif.get();

        uint32_t npages = static_cast<uint32_t>(data.nslices());
        uint32_t height = static_cast<uint32_t>(data.nrows());
        uint32_t width = static_cast<uint32_t>(data.ncols());

        std::vector<float> buf(width);

        for (uint32_t i = 0; i < npages; i++) {
            TIFFSetField(tif_, TIFFTAG_IMAGEWIDTH, width);
            TIFFSetField(tif_, TIFFTAG_IMAGELENGTH, height);
            TIFFSetField(tif_, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(tif_, TIFFTAG_BITSPERSAMPLE, 32);
            TIFFSetField(tif_, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
            TIFFSetField(tif_, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(tif_, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tif_, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            TIFFSetField(tif_, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(tif_, TIFFTAG_ROWSPERSTRIP, height); // one strip

            for (uint32_t j = 0; j < height; j++) {
                for (uint32_t k = 0; k < width; k++) buf[k] = data[{i, j, k}];

                if (TIFFWriteScanline(tif_, buf.data(), j) < 0) {
                    throw std::runtime_error("tiff: failed to write data to " +
                                             filename);
                }
            }
            TIFFWriteDirectory(tif_);
        }
    }

    template <typename T>
    inline void write_vectors(std::string filename,
                              const std::array<Array<T>, 3> &data) {

        // open file
        auto tif = detail::open(filename, "w");
        TIFF *tif_ = tif.get();

        uint32_t npages = static_cast<uint32_t>(data[0].nslices());
        uint32_t height = static_cast<uint32_t>(data[0].nrows());
        uint32_t width = static_cast<uint32_t>(data[0].ncols());

        for (uint32_t i = 0; i < npages; i++) {
            TIFFSetField(tif_, TIFFTAG_IMAGEWIDTH, width);
            TIFFSetField(tif_, TIFFTAG_IMAGELENGTH, height);
            TIFFSetField(tif_, TIFFTAG_SAMPLESPERPIXEL, 3);
            TIFFSetField(tif_, TIFFTAG_BITSPERSAMPLE, 32);
            TIFFSetField(tif_, TIFFTAG_SAMPLEFORMAT,
                         std::is_floating_point<T>::value ? SAMPLEFORMAT_IEEEFP
                                                          : SAMPLEFORMAT_UINT);
            TIFFSetField(tif_, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(tif_, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tif_, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
            TIFFSetField(tif_, TIFFTAG_COMPRESSION, COMPRESSION_LZW);

            // Metadata
            TIFFSetField(
                tif_, TIFFTAG_IMAGEDESCRIPTION,
                "3D vector field: x,y,z components; layout: [slice, row, col]");

            std::vector<T> buffer(width * 3);
            for (uint32_t j = 0; j < height; ++j) {
                for (uint32_t k = 0; k < width; ++k) {
                    buffer[3 * k + 0] = data[0][{i, j, k}];
                    buffer[3 * k + 1] = data[1][{i, j, k}];
                    buffer[3 * k + 2] = data[2][{i, j, k}];
                }

                if (TIFFWriteScanline(tif_, buffer.data(), j) < 0) {
                    throw std::runtime_error("tiff: failed to write data to " +
                                             filename);
                }
            }
            TIFFWriteDirectory(tif_);
        }
    }

    template <typename T>
    void write3(const std::string &name, const std::array<Array<T>, 3> &data) {
        // trucate extension if present
        std::string basename = name;
        // check for existing .tif or .tiff extension
        if (name.rfind(".tif") != std::string::npos) {
            basename = name.substr(0, name.rfind(".tif"));
        }
        for (size_t i = 0; i < 3; ++i) {
            std::string filename = basename + std::to_string(i) + ".tiff";
            tiff::write(filename, data[i]);
        }
    }
} // namespace tomocam::tiff
#endif // TOMOCAM_TIFF__H
