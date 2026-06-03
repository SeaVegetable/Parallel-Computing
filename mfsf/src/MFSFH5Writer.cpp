#include "MFSFH5Writer.hpp"

#include <array>
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

hid_t CreateGroupOrThrow(const hid_t parent_id, const std::string &name)
{
    const hid_t group_id = H5Gcreate2(parent_id, name.c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0)
    {
        throw std::runtime_error("Failed to create HDF5 group " + name);
    }
    return group_id;
}

void WriteIntDataset(const hid_t group_id,
    const std::string &name,
    const std::vector<int> &values,
    const std::vector<hsize_t> &dims)
{
    const hid_t dataspace_id = H5Screate_simple(static_cast<int>(dims.size()),
        dims.data(), nullptr);
    if (dataspace_id < 0)
    {
        throw std::runtime_error("Failed to create int dataspace for " + name);
    }

    const hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_INT,
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

void WriteDoubleDataset(const hid_t group_id,
    const std::string &name,
    const std::vector<double> &values,
    const std::vector<hsize_t> &dims)
{
    const hid_t dataspace_id = H5Screate_simple(static_cast<int>(dims.size()),
        dims.data(), nullptr);
    if (dataspace_id < 0)
    {
        throw std::runtime_error("Failed to create double dataspace for " + name);
    }

    const hid_t dataset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_DOUBLE,
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

template<std::size_t N>
std::vector<int> ToVector(const std::array<int, N> &values)
{
    return std::vector<int>(values.begin(), values.end());
}

template<std::size_t N>
int Product(const std::array<int, N> &values)
{
    int product = 1;
    for (std::size_t axis = 0; axis < N; ++axis)
    {
        product *= values[axis];
    }
    return product;
}

template<int Dim>
void WritePrecomputeData(const std::string &filename,
    const MFSFPrecomputeData<Dim> &data)
{
    const hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
        H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
    {
        throw std::runtime_error("Failed to create HDF5 file " + filename);
    }

    try
    {
        const hid_t meta_group = CreateGroupOrThrow(file_id, "/meta");
        const hid_t topology_group = CreateGroupOrThrow(file_id, "/topology");
        const hid_t mesh_group = CreateGroupOrThrow(file_id, "/mesh");
        const hid_t quadrature_group = CreateGroupOrThrow(file_id, "/quadrature");
        const hid_t geometry_group = CreateGroupOrThrow(file_id, "/geometry");

        try
        {
            WriteIntDataset(meta_group, "dim", std::vector<int>(1, Dim),
                std::vector<hsize_t>(1, 1));
            WriteIntDataset(meta_group, "degrees", ToVector(data.degrees),
                std::vector<hsize_t>(1, static_cast<hsize_t>(Dim)));
            WriteIntDataset(meta_group, "num_funcs", ToVector(data.num_funcs),
                std::vector<hsize_t>(1, static_cast<hsize_t>(Dim)));
            WriteIntDataset(meta_group, "num_elements", ToVector(data.num_elements),
                std::vector<hsize_t>(1, static_cast<hsize_t>(Dim)));
            WriteIntDataset(meta_group, "nqp", ToVector(data.nqp),
                std::vector<hsize_t>(1, static_cast<hsize_t>(Dim)));

            const int total_elements = Product(data.num_elements);
            const int total_qp = Product(data.nqp);
            const int total_control_points = Product(data.num_funcs);
            int num_local_basis = 1;
            for (int axis = 0; axis < Dim; ++axis)
            {
                num_local_basis *= (data.degrees[axis] + 1);
            }

            WriteIntDataset(topology_group, "id", data.id,
                std::vector<hsize_t>(1, static_cast<hsize_t>(data.id.size())));
            WriteIntDataset(topology_group, "ien", data.ien,
                std::vector<hsize_t>{
                    static_cast<hsize_t>(total_elements),
                    static_cast<hsize_t>(num_local_basis)
                });
            WriteIntDataset(topology_group, "num_local_basis",
                std::vector<int>(1, num_local_basis),
                std::vector<hsize_t>(1, 1));

            WriteDoubleDataset(mesh_group, "control_points", data.control_points,
                std::vector<hsize_t>{
                    static_cast<hsize_t>(total_control_points),
                    static_cast<hsize_t>(Dim)
                });

            for (int axis = 0; axis < Dim; ++axis)
            {
                WriteDoubleDataset(mesh_group,
                    "element_size_axis_" + std::to_string(axis),
                    data.element_sizes[axis],
                    std::vector<hsize_t>{
                        static_cast<hsize_t>(data.element_sizes[axis].size())
                    });
                WriteDoubleDataset(quadrature_group,
                    "points_axis_" + std::to_string(axis),
                    data.quadrature_points[axis],
                    std::vector<hsize_t>{
                        static_cast<hsize_t>(data.quadrature_points[axis].size())
                    });
                WriteDoubleDataset(quadrature_group,
                    "weights_axis_" + std::to_string(axis),
                    data.quadrature_weights[axis],
                    std::vector<hsize_t>{
                        static_cast<hsize_t>(data.quadrature_weights[axis].size())
                    });
            }

            WriteDoubleDataset(geometry_group, "inv_jacobian", data.inv_jacobian,
                std::vector<hsize_t>{
                    static_cast<hsize_t>(total_elements),
                    static_cast<hsize_t>(total_qp),
                    static_cast<hsize_t>(Dim),
                    static_cast<hsize_t>(Dim)
                });
            WriteDoubleDataset(geometry_group, "det_jacobian", data.det_jacobian,
                std::vector<hsize_t>{
                    static_cast<hsize_t>(total_elements),
                    static_cast<hsize_t>(total_qp)
                });
            H5Gclose(geometry_group);
            H5Gclose(quadrature_group);
            H5Gclose(mesh_group);
            H5Gclose(topology_group);
            H5Gclose(meta_group);
        }
        catch (...)
        {
            H5Gclose(geometry_group);
            H5Gclose(quadrature_group);
            H5Gclose(mesh_group);
            H5Gclose(topology_group);
            H5Gclose(meta_group);
            throw;
        }

        H5Fclose(file_id);
    }
    catch (...)
    {
        H5Fclose(file_id);
        throw;
    }
}
}

void MFSFH5Writer::Write(const std::string &filename,
    const MFSFPrecomputeData<2> &data) const
{
    WritePrecomputeData(filename, data);
}

void MFSFH5Writer::Write(const std::string &filename,
    const MFSFPrecomputeData<3> &data) const
{
    WritePrecomputeData(filename, data);
}
