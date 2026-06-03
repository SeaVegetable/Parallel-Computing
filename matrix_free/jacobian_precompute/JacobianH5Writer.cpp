#include "JacobianH5Writer.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <hdf5.h>

namespace
{
void ThrowIfError(const herr_t status, const std::string &message)
{
    if (status < 0)
    {
        throw std::runtime_error(message);
    }
}

void WriteIntDataset(const hid_t group_id,
    const std::string &name,
    const std::vector<int> &values)
{
    const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    if (dataspace_id < 0) throw std::runtime_error("Failed to create int dataspace for " + name);

    hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_INT,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        H5Sclose(dataspace_id);
        throw std::runtime_error("Failed to create int dataset " + name);
    }

    ThrowIfError(H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, values.data()), "Failed to write int dataset " + name);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
}

void WriteDoubleDataset1D(const hid_t group_id,
    const std::string &name,
    const std::vector<double> &values)
{
    const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    if (dataspace_id < 0) throw std::runtime_error("Failed to create double dataspace for " + name);

    hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_DOUBLE,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        H5Sclose(dataspace_id);
        throw std::runtime_error("Failed to create double dataset " + name);
    }

    ThrowIfError(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, values.data()), "Failed to write double dataset " + name);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
}

void WriteDoubleDataset2D(const hid_t group_id,
    const std::string &name,
    const std::vector<double> &values,
    const hsize_t dim0,
    const hsize_t dim1)
{
    const hsize_t dims[2] = {dim0, dim1};
    hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
    if (dataspace_id < 0) throw std::runtime_error("Failed to create 2D dataspace for " + name);

    hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_DOUBLE,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        H5Sclose(dataspace_id);
        throw std::runtime_error("Failed to create 2D dataset " + name);
    }

    ThrowIfError(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, values.data()), "Failed to write 2D dataset " + name);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
}

void WriteDoubleDataset4D(const hid_t group_id,
    const std::string &name,
    const std::vector<double> &values,
    const hsize_t dim0,
    const hsize_t dim1,
    const hsize_t dim2,
    const hsize_t dim3)
{
    const hsize_t dims[4] = {dim0, dim1, dim2, dim3};
    hid_t dataspace_id = H5Screate_simple(4, dims, nullptr);
    if (dataspace_id < 0) throw std::runtime_error("Failed to create 4D dataspace for " + name);

    hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_DOUBLE,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        H5Sclose(dataspace_id);
        throw std::runtime_error("Failed to create 4D dataset " + name);
    }

    ThrowIfError(H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, values.data()), "Failed to write 4D dataset " + name);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
}
}

void JacobianH5Writer::Write(const std::string &filename,
    const int &p,
    const int &q,
    const int &nlocalelemx,
    const int &nlocalelemy,
    const int &nqp1,
    const int &nqp2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    const std::vector<double> &qp1,
    const std::vector<double> &qp2,
    const std::vector<double> &w1,
    const std::vector<double> &w2,
    const std::vector<double> &inv_jacobian,
    const std::vector<double> &det_jacobian) const
{
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
    {
        throw std::runtime_error("Failed to create HDF5 file " + filename);
    }

    hid_t meta_group = H5Gcreate2(file_id, "/meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t quadrature_group = H5Gcreate2(file_id, "/quadrature", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t mesh_group = H5Gcreate2(file_id, "/mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t geometry_group = H5Gcreate2(file_id, "/geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (meta_group < 0 || quadrature_group < 0 || mesh_group < 0 || geometry_group < 0)
    {
        if (meta_group >= 0) H5Gclose(meta_group);
        if (quadrature_group >= 0) H5Gclose(quadrature_group);
        if (mesh_group >= 0) H5Gclose(mesh_group);
        if (geometry_group >= 0) H5Gclose(geometry_group);
        H5Fclose(file_id);
        throw std::runtime_error("Failed to create groups in HDF5 file " + filename);
    }

    WriteIntDataset(meta_group, "p", std::vector<int>(1, p));
    WriteIntDataset(meta_group, "q", std::vector<int>(1, q));
    WriteIntDataset(meta_group, "nlocalelemx", std::vector<int>(1, nlocalelemx));
    WriteIntDataset(meta_group, "nlocalelemy", std::vector<int>(1, nlocalelemy));
    WriteIntDataset(meta_group, "nqp1", std::vector<int>(1, nqp1));
    WriteIntDataset(meta_group, "nqp2", std::vector<int>(1, nqp2));

    WriteDoubleDataset1D(quadrature_group, "points_x", qp1);
    WriteDoubleDataset1D(quadrature_group, "points_y", qp2);
    WriteDoubleDataset1D(quadrature_group, "weights_x", w1);
    WriteDoubleDataset1D(quadrature_group, "weights_y", w2);

    WriteDoubleDataset1D(mesh_group, "elem_size1", elem_size1);
    WriteDoubleDataset1D(mesh_group, "elem_size2", elem_size2);

    const hsize_t nelem = static_cast<hsize_t>(nlocalelemx * nlocalelemy);
    const hsize_t nqp = static_cast<hsize_t>(nqp1 * nqp2);

    WriteDoubleDataset4D(geometry_group, "inv_jacobian", inv_jacobian, nelem, nqp, 2, 2);
    WriteDoubleDataset2D(geometry_group, "det_jacobian", det_jacobian, nelem, nqp);

    H5Gclose(geometry_group);
    H5Gclose(mesh_group);
    H5Gclose(quadrature_group);
    H5Gclose(meta_group);
    H5Fclose(file_id);
}
