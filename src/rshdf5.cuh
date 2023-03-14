/// rshdf5.h - Header file for HDF5 functions
/// Marc Brooker, 03 November 2006
/// Edited by Yaaseen Martin, 27 August 2019

#ifndef __RS_HDF5_H
#define __RS_HDF5_H

#include <complex>

/// For use of hid_t datatype
extern "C" {
    #include <hdf5.h>
    #include <H5LTpublic.h>
}

namespace rshdf5 {

    ///Open the HDF5 file for writing
    hid_t CreateFile(const std::string name);

    ///Add a dataset to the HDF5 file
    void AddChunkToFile(hid_t file, std::complex<rsFloat> *data, unsigned int size, rsFloat time, rsFloat rate, rsFloat fullscale, unsigned int count);

    ///Close the HDF5 file
    void CloseFile(hid_t file);

    /// Read an antenna gain pattern or RCS pattern from a file
    rsFloat **ReadPattern(const std::string &name, const std::string &dataset_name, unsigned int &azi_size, unsigned int &elev_size);

}

#endif
