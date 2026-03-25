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

#ifndef WRITE_VTI__H
#define WRITE_VTI__H

#include <array>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "array.h"

namespace tomocam::vti {
    template <typename T>
    void write_vectors(const std::string &filename,
                       const std::array<Array<T>, 3> &vectors,
                       const std::string &field_name = "magetization") {

        // replace file extension with .vti
        std::filesystem::path filepath(filename);
        filepath.replace_extension(".vti");

        std::ofstream file(filepath, std::ios::binary);
        if (!file) { throw std::runtime_error("Could not open file for writing"); }

        size_t nz = vectors[0].nslices();
        size_t ny = vectors[0].nrows();
        size_t nx = vectors[0].ncols();
        const size_t npoints = nx * ny * nz;

        // define origin at the center of the volume
        double ox = -0.5 * static_cast<double>(nx);
        double oy = -0.5 * static_cast<double>(ny);
        double oz = -0.5 * static_cast<double>(nz);

        // define spacing as 1.0
        double dx = 1.0;
        double dy = 1.0;
        double dz = 1.0;

        // Header
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"ImageData\" version=\"0.1\" "
                "byte_order=\"LittleEndian\">\n";

        file << std::format("  <ImageData WholeExtent=\"0 {} 0 {} 0 {}\" "
                            "Origin=\"{} {} {}\" Spacing=\"{} {} {}\">\n",
                            nx - 1, ny - 1, nz - 1, ox, oy, oz, dx, dy, dz);

        file << std::format("    <Piece Extent=\"0 {} 0 {} 0 {}\">\n", nx - 1,
                            ny - 1, nz - 1);

        // Point data with appended binary format
        const size_t data_size = npoints * 3 * sizeof(double);
        file << std::format("      <PointData Vectors=\"{}\">\n", field_name);
        file << std::format(
            "        <DataArray type=\"Float64\" Name=\"{}\" "
            "NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>\n",
            field_name);
        file << "      </PointData>\n";

        // Required empty cell data
        file << "      <CellData></CellData>\n";

        // Footer
        file << "    </Piece>\n";
        file << "  </ImageData>\n";

        // Appended data section
        file << "  <AppendedData encoding=\"raw\">\n";
        file << "   _";

        // Write size header (32-bit unsigned integer)
        uint32_t byte_count = static_cast<uint32_t>(data_size);
        file.write(reinterpret_cast<const char *>(&byte_count), sizeof(uint32_t));

        // Write binary vector data interleaved
        std::vector<double> buffer(3);
        for (size_t i = 0; i < npoints; ++i) {
            buffer[0] = static_cast<double>(vectors[0][i]);
            buffer[1] = static_cast<double>(vectors[1][i]);
            buffer[2] = static_cast<double>(vectors[2][i]);
            file.write(reinterpret_cast<const char *>(buffer.data()),
                       3 * sizeof(double));
        }

        file << "\n  </AppendedData>\n";
        file << "</VTKFile>\n";
    }
} // namespace tomocam::vti
#endif // WRITE_VTI__H
