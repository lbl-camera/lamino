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
#include <format>
#include <iostream>
#include <string>
#include <tiffio.h>

#include "array.h"

namespace tomocam::tiff {

    inline Array<float> read(std::string filename) {

        TIFF *tif_ = TIFFOpen(filename.c_str(), "r");
        if (tif_ == nullptr) {
            throw std::runtime_error(
                std::format("Error: failed to open TIFF file: {}", filename));
        }

        // count number of projections
        uint32_t npages = 0;
        do { npages++; } while (TIFFReadDirectory(tif_));
        if (npages == 0) {
            TIFFClose(tif_);
            throw std::runtime_error("Error: no projections found in tiff file.");
        }
        TIFFSetDirectory(tif_, 0); // reset to first page

        // get image size
        uint16_t bits, format;
        uint32_t w, h;
        if (!TIFFGetField(tif_, TIFFTAG_IMAGEWIDTH, &w) ||
            !TIFFGetField(tif_, TIFFTAG_IMAGELENGTH, &h) ||
            !TIFFGetField(tif_, TIFFTAG_BITSPERSAMPLE, &bits)) {
            TIFFClose(tif_);
            throw std::runtime_error("Error: failed to read TIFF metadata");
        }
        try {
            TIFFGetField(tif_, TIFFTAG_SAMPLEFORMAT, &format);
        } catch (const std::exception &) {
            format = SAMPLEFORMAT_UINT; // default to uint if not specified
        }
        // force garbage format value to uint
        if (format != SAMPLEFORMAT_UINT && format != SAMPLEFORMAT_IEEEFP &&
            format != SAMPLEFORMAT_INT && format != SAMPLEFORMAT_VOID &&
            format != SAMPLEFORMAT_COMPLEXINT &&
            format != SAMPLEFORMAT_COMPLEXIEEEFP) {
            format = SAMPLEFORMAT_UINT;
        }

        // allocate memory
        size_t nscls = static_cast<size_t>(npages);
        size_t nrows = static_cast<size_t>(h);
        size_t ncols = static_cast<size_t>(w);
        // ensure odd dimensions for centered-fft
        if (nrows % 2 == 0) { nrows -= 1; }
        if (ncols % 2 == 0) { ncols -= 1; }
        size_t ptr_shift = bits / 8;
        Array<float> data(nscls, nrows, ncols);

        tsize_t line_size = TIFFScanlineSize(tif_);
        void *buf = _TIFFmalloc(line_size);
        if (buf == nullptr) {
            TIFFClose(tif_);
            throw std::runtime_error(
                "Error: failed to allocate buffer for scanline");
        }

        // Read scanlines and convert to float
        for (size_t i = 0; i < nscls; i++) {
            TIFFSetDirectory(tif_, static_cast<tdir_t>(i));
            for (size_t j = 0; j < nrows; j++) {
                if (TIFFReadScanline(tif_, buf, static_cast<uint32_t>(j)) < 0) {
                    _TIFFfree(buf);
                    TIFFClose(tif_);
                    throw std::runtime_error(
                        std::format("Error reading scanline {} of page {}.", j, i));
                }

                for (size_t k = 0; k < ncols; k++) {
                    if (format == SAMPLEFORMAT_UINT) {
                        if (bits == 8) {
                            data[{i, j, k}] =
                                static_cast<float>(((uint8_t *)buf)[k]);
                        } else if (bits == 16) {
                            data[{i, j, k}] =
                                static_cast<float>(((uint16_t *)buf)[k]);
                        } else if (bits == 32) {
                            data[{i, j, k}] =
                                static_cast<float>(((uint32_t *)buf)[k]);
                        }
                    } else if (format == SAMPLEFORMAT_IEEEFP) {
                        if (bits == 32) {
                            data[{i, j, k}] = ((float *)buf)[k];
                        } else if (bits == 64) {
                            data[{i, j, k}] = static_cast<float>(((double *)buf)[k]);
                        }
                    } else {
                        _TIFFfree(buf);
                        TIFFClose(tif_);
                        throw std::runtime_error(
                            std::format("Error: unsupported format/bits "
                                        "combination: format={}, bits={}",
                                        format, bits));
                    }
                }
            }
        }
        _TIFFfree(buf);
        TIFFClose(tif_);
        return data;
    }

    inline void write(std::string filename, const Array<float> &data) {

        // open file
        TIFF *tif_ = TIFFOpen(filename.c_str(), "w");

        uint32_t npages = static_cast<uint32_t>(data.nslices());
        uint32_t height = static_cast<uint32_t>(data.nrows());
        uint32_t width = static_cast<uint32_t>(data.ncols());

        float *buf = (float *)_TIFFmalloc(sizeof(float) * width);

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

                if (TIFFWriteScanline(tif_, buf, j) < 0) {
                    std::cerr << "Error wrtiting data to tif file." << std::endl;
                    std::exit(2);
                }
            }
            TIFFWriteDirectory(tif_);
        }
        _TIFFfree(buf);
        TIFFClose(tif_);
    }
} // namespace tomocam::tiff
#endif // TOMOCAM_TIFF__H
