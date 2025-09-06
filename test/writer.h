#include <H5Cpp.h>
#include <fstream>
#include <iostream>

#include "array.h"

#ifndef TOMOCAM_WRITER__H
#define TOMOCAM_WRITER__H

namespace tomocam {

    class H5Writer {

        private:
            H5::H5File fp_;

        public:
            H5Writer(const char *filename): fp_(H5::H5File(filename, H5F_ACC_TRUNC)) {}
            ~H5Writer() { fp_.close(); }

            void write(const char *dname, Array<float> &array) {

                // create dataspace
                auto d = array.dims();
                hsize_t dims[2];
                dims[0] = d.x;
                dims[1] = d.y;
                H5::DataSpace spc(2, dims);
                H5::DataSet dset = fp_.createDataSet(dname, H5::PredType::NATIVE_FLOAT, spc);
                dset.write(array.begin(), H5::PredType::NATIVE_FLOAT);
            }

            void write(const char *dname, Array<double> &array) {

                // create dataspace
                auto d = array.dims();
                hsize_t dims[2];
                dims[0] = d.x;
                dims[1] = d.y;
                H5::DataSpace spc(2, dims);
                H5::DataSet dset = fp_.createDataSet(dname, H5::PredType::NATIVE_DOUBLE, spc);
                dset.write(array.begin(), H5::PredType::NATIVE_DOUBLE);
            }

    };
} // namespace tomocam
#endif // TOMOCAM_WRITER__H
