/// rshdf5.cpp - Interact with the HDF5 file format
/// Marc Brooker, 03 November 2006
/// Edited by Yaaseen Martin, 27 August 2019

#include <config.h>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include "rshdf5.cuh"
#include "rsparameters.cuh"
#include "rsdebug.cuh"

using namespace rshdf5;

/// For use of hid_t datatype
extern "C" {
    #include <hdf5.h>
    #include <H5LTpublic.h>
}

/// Open the HDF5 file for reading
hid_t OpenFile(const std::string &name)
{
    // Read-only access; default I/O access parameters
    hid_t file = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Throw an error if file cannot be opened
    if (file < 0)
        throw std::runtime_error("ERROR: Could not open HDF5 file " + name + " to read pulse.");
    return file;
}

/// Open the HDF5 file for writing
hid_t rshdf5::CreateFile(const std::string name)
{
    // Read and write access; overwrite existing file; default I/O access parameters
    hid_t file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Throw an error if file cannot be created
    if (file < 0)
        throw std::runtime_error("ERROR: Could not create HDF5 file " + name + " for export");
    return file;
}

/// Add a dataset to the HDF5 file
void rshdf5::AddChunkToFile(hid_t file, std::complex<rsFloat> *data, unsigned int size, rsFloat time, rsFloat rate, rsFloat fullscale,
                            unsigned int count)
{
    // Create the name of the dataset
    std::ostringstream oss;
    oss << "chunk_" << std::setw(6) << std::setfill('0') << count;
    std::string I_chunk_name = oss.str()+"_I";
    std::string Q_chunk_name = oss.str()+"_Q";

    // Create the size variable needed by the HDF5 Lite API
    hsize_t datasize = size;

    // Separate I and Q data
    double *I = new double[size];
    double *Q = new double[size];
    for (unsigned int i = 0; i < size; i++) {
        I[i] = data[i].real();
        Q[i] = data[i].imag();
    }

    // Create the dataset using the HDF5 Lite API
    if (H5LTmake_dataset_double(file, I_chunk_name.c_str(), 1, &datasize, I) < 0)
        throw std::runtime_error("ERROR: Error while writing I data to HDF5 file");
    if (H5LTmake_dataset_double(file, Q_chunk_name.c_str(), 1, &datasize, Q) < 0)
        throw std::runtime_error("ERROR: Error while writing Q data to HDF5 file");

    // Add attributes to the dataset, with the attributes of the response
    if (H5LTset_attribute_double(file, I_chunk_name.c_str(), "time", &time, 1) < 0)
        throw std::runtime_error("ERROR: Error while setting attribute \"time\" on chunk " + I_chunk_name);
    if (H5LTset_attribute_double(file, I_chunk_name.c_str(), "rate", &rate, 1) < 0)
        throw std::runtime_error("ERROR: Error while setting attribute \"rate\" on chunk " + I_chunk_name);
    if (H5LTset_attribute_double(file, I_chunk_name.c_str(), "fullscale", &fullscale, 1) < 0)
        throw std::runtime_error("ERROR: Error while setting attribute \"fullscale\" on chunk " + I_chunk_name);
    if (H5LTset_attribute_double(file, Q_chunk_name.c_str(), "time", &time, 1) < 0)
        throw std::runtime_error("ERROR: Error while setting attribute \"time\" on chunk " + Q_chunk_name);
    if (H5LTset_attribute_double(file, Q_chunk_name.c_str(), "rate", &rate, 1) < 0)
        throw std::runtime_error("ERROR: Error while setting attribute \"rate\" on chunk " + Q_chunk_name);
    if (H5LTset_attribute_double(file, Q_chunk_name.c_str(), "fullscale", &fullscale, 1) < 0)
        throw std::runtime_error("ERROR: Error while setting attribute \"fullscale\" on chunk " + Q_chunk_name);

    // Free the buffers
    delete[] I;
    delete[] Q;
}

/// Close the HDF5 file
void rshdf5::CloseFile(hid_t file) {
    if (H5Fclose(file) < 0)
        throw std::runtime_error("ERROR: Error while closing HDF5 file");
}

/// Read an antenna gain pattern or RCS pattern from a file
rsFloat **rshdf5::ReadPattern(const std::string &name, const std::string &dataset_name, unsigned int &azi_size, unsigned int &elev_size)
{
    hid_t file_id;
    int rank;
    hsize_t dims[2];
    size_t type_size;
    H5T_class_t data_class;

    // Load the HDF5 file
    file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0)
        throw std::runtime_error("ERROR: Cannot open HDF5 file "+name + " to read antenna data");

    // Get the rank of the dataset
    herr_t err = H5LTget_dataset_ndims(file_id, dataset_name.c_str(), &rank);
    if (err < 0)
        throw std::runtime_error("ERROR: Could not get rank of dataset \""+dataset_name+"\" in file " + name);
    else if (rank != 2)
        throw std::runtime_error("ERROR: Dataset \""+dataset_name+"\" in file " + name + " does not have rank 2");

    // Get the dimensions of the file
    err = H5LTget_dataset_info(file_id, dataset_name.c_str(), &(dims[0]), &data_class, &type_size);
    if (err < 0)
        throw std::runtime_error("ERROR: Could not get dimensions of dataset \""+dataset_name+"\" in file " + name);
    if (type_size != sizeof(float))
        throw std::runtime_error("ERROR: Type size incorrect in dataset \""+dataset_name+"\" in file " + name);

    // Allocate memory for the pattern
    float *data = new float[dims[0]*dims[1]];

    // Load the pattern into memory
    err = H5LTread_dataset_float(file_id, dataset_name.c_str(), data);
    if (err < 0) {
        delete[] data;
        throw std::runtime_error("ERROR: Could not read float data from dataset \""+dataset_name+"\" in file" + name);
    }

    // Close the HDF5 file
    if (H5Fclose(file_id) < 0)
        throw std::runtime_error("ERROR: Error while closing HDF5 file "+name);

    // Break the data down into a 2D array
    azi_size = dims[0];
    elev_size = dims[1];
    rsFloat **ret = new rsFloat*[azi_size];
    for (unsigned int i = 0; i < azi_size; i++) {
        ret[i] = new rsFloat[elev_size];
        for (unsigned int j = 0; j < elev_size; j++)
            ret[i][j] = data[i*azi_size+j];
    }

    // Clean up
    delete[] data;
    return ret;
}
