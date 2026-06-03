#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "BernsteinBasis.hpp"
#include "ElementMF.hpp"
#include "FileManager.hpp"
#include "QuadraturePoint.hpp"
#include "jacobian_precompute/JacobianH5Writer.hpp"

namespace
{
std::string GetDirectoryName(const std::string &path)
{
    const std::string::size_type pos = path.find_last_of("/\\");
    if (pos == std::string::npos)
    {
        return ".";
    }
    if (pos == 0)
    {
        return path.substr(0, 1);
    }
    return path.substr(0, pos);
}

std::string JoinPath(const std::string &lhs, const std::string &rhs)
{
    if (lhs.empty() || lhs == ".")
    {
        return rhs;
    }
    if (lhs[lhs.size() - 1] == '/')
    {
        return lhs + rhs;
    }
    return lhs + "/" + rhs;
}

bool FileExists(const std::string &path)
{
    std::ifstream file(path.c_str());
    return file.good();
}

void EnsureDirectory(const std::string &path)
{
    struct stat info;
    if (stat(path.c_str(), &info) == 0)
    {
        if ((info.st_mode & S_IFDIR) == 0)
        {
            throw std::runtime_error(path + " exists but is not a directory");
        }
        return;
    }

    if (mkdir(path.c_str(), 0755) != 0)
    {
        throw std::runtime_error("Failed to create directory " + path);
    }
}

std::string BuildOutputFilename(const std::string &partition_filename)
{
    const std::string suffix = ".txt";
    const std::string output_suffix = "_jacobian.h5";
    if (partition_filename.size() >= suffix.size()
        && partition_filename.substr(partition_filename.size() - suffix.size()) == suffix)
    {
        return partition_filename.substr(0, partition_filename.size() - suffix.size()) + output_suffix;
    }
    return partition_filename + output_suffix;
}
}

int main(int argc, char *argv[])
{
    try
    {
        const std::string info_file = (argc > 1) ? argv[1] : "info.txt";
        const std::string input_dir = (argc > 2) ? argv[2] : GetDirectoryName(info_file);
        const std::string output_dir = (argc > 3) ? argv[3] : JoinPath(input_dir, "jacobian_precompute_h5");

        EnsureDirectory(output_dir);

        int p, q, nElemX, nElemY, part_num_1d, dim;
        double Lx, Ly;
        std::string base_name;

        FileManager fm;
        fm.ReadPreprocessInfo(info_file, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

        QuadraturePoint quad1(p + 1, 0, 1);
        QuadraturePoint quad2(q + 1, 0, 1);
        const std::vector<double> qp1 = quad1.GetQuadraturePoint();
        const std::vector<double> qp2 = quad2.GetQuadraturePoint();
        const std::vector<double> w1 = quad1.GetWeight();
        const std::vector<double> w2 = quad2.GetWeight();

        BernsteinBasis bern1(p);
        BernsteinBasis bern2(q);
        bern1.GenerateBernsteinBasis(&quad1);
        bern2.GenerateBernsteinBasis(&quad2);
        const std::vector<double> &B1 = bern1.GetBernsteinBasis();
        const std::vector<double> &B2 = bern2.GetBernsteinBasis();
        const std::vector<double> &dB1 = bern1.GetBernsteinBasisDerivative();
        const std::vector<double> &dB2 = bern2.GetBernsteinBasisDerivative();

        ElementMF elemmf(p, q);
        JacobianH5Writer writer;

        const int nqp1 = quad1.GetNumQuadraturePoint();
        const int nqp2 = quad2.GetNumQuadraturePoint();
        const int nqp = nqp1 * nqp2;
        const int npartition = part_num_1d * part_num_1d;
        int processed_count = 0;

        for (int rank = 0; rank < npartition; ++rank)
        {
            const std::string partition_filename = fm.GetPartitionFilename(base_name, rank);
            const std::string partition_path = JoinPath(input_dir, partition_filename);

            if (!FileExists(partition_path))
            {
                std::cout << "Skipping missing partition file: " << partition_path << std::endl;
                continue;
            }

            int nlocalfunc;
            int nlocalelemx;
            int nlocalelemy;
            std::vector<int> ghostID;
            std::vector<double> CP;
            std::vector<int> ID;
            std::vector<int> IEN;
            std::vector<int> Dir;
            std::vector<double> elem_size1;
            std::vector<double> elem_size2;
            std::vector<double> NURBSExtraction1;
            std::vector<double> NURBSExtraction2;

            fm.ReadPartition(partition_path, nlocalfunc,
                nlocalelemx, nlocalelemy,
                elem_size1, elem_size2,
                CP, ID, ghostID, Dir, IEN,
                NURBSExtraction1, NURBSExtraction2);

            const int nLocBas = elemmf.GetNumLocalBasis();
            const int pp = elemmf.GetNumLocalBasis1D(0);
            const int qq = elemmf.GetNumLocalBasis1D(1);
            const int nelem = nlocalelemx * nlocalelemy;

            std::vector<double> inv_jacobian(nelem * nqp * 4, 0.0);
            std::vector<double> det_jacobian(nelem * nqp, 0.0);
            std::vector<double> eCP(2 * nLocBas, 0.0);
            std::vector<double> eNURBSExtraction1(pp * pp, 0.0);
            std::vector<double> eNURBSExtraction2(qq * qq, 0.0);

            for (int jy = 0; jy < nlocalelemy; ++jy)
            {
                for (int ix = 0; ix < nlocalelemx; ++ix)
                {
                    const int elem_index = jy * nlocalelemx + ix;
                    for (int local_basis = 0; local_basis < nLocBas; ++local_basis)
                    {
                        const int ien_index = IEN[elem_index * nLocBas + local_basis];
                        eCP[2 * local_basis] = CP[2 * ien_index];
                        eCP[2 * local_basis + 1] = CP[2 * ien_index + 1];
                    }

                    std::copy(NURBSExtraction1.begin() + ix * pp * pp,
                        NURBSExtraction1.begin() + (ix + 1) * pp * pp,
                        eNURBSExtraction1.begin());
                    std::copy(NURBSExtraction2.begin() + jy * qq * qq,
                        NURBSExtraction2.begin() + (jy + 1) * qq * qq,
                        eNURBSExtraction2.begin());

                    elemmf.SetElement(eNURBSExtraction1, eNURBSExtraction2, elem_size1[ix], elem_size2[jy]);

                    for (int qy = 0; qy < nqp2; ++qy)
                    {
                        const std::vector<double> b2(B2.begin() + qy * (q + 1), B2.begin() + (qy + 1) * (q + 1));
                        const std::vector<double> db2(dB2.begin() + qy * (q + 1), dB2.begin() + (qy + 1) * (q + 1));

                        for (int qx = 0; qx < nqp1; ++qx)
                        {
                            const std::vector<double> b1(B1.begin() + qx * (p + 1), B1.begin() + (qx + 1) * (p + 1));
                            const std::vector<double> db1(dB1.begin() + qx * (p + 1), dB1.begin() + (qx + 1) * (p + 1));

                            std::array<double, 4> jacobian_matrix{};
                            std::array<double, 4> inv_jacobian_matrix{};
                            double detJ = 0.0;

                            elemmf.ComputeJacobianDataSingleQP(b1, b2, db1, db2, eCP,
                                jacobian_matrix, inv_jacobian_matrix, detJ);

                            const int qp_index = qy * nqp1 + qx;
                            const int det_index = elem_index * nqp + qp_index;
                            const int inv_index = det_index * 4;

                            det_jacobian[det_index] = detJ;
                            for (int entry = 0; entry < 4; ++entry)
                            {
                                inv_jacobian[inv_index + entry] = inv_jacobian_matrix[entry];
                            }
                        }
                    }
                }
            }

            const std::string output_path = JoinPath(output_dir, BuildOutputFilename(partition_filename));
            writer.Write(output_path, p, q, nlocalelemx, nlocalelemy, nqp1, nqp2,
                elem_size1, elem_size2, qp1, qp2, w1, w2, inv_jacobian, det_jacobian);

            ++processed_count;
            std::cout << "Wrote " << output_path << std::endl;
        }

        if (processed_count == 0)
        {
            std::cerr << "No partition files were processed in " << input_dir << std::endl;
            return 1;
        }

        std::cout << "Processed " << processed_count << " partition files." << std::endl;
        return 0;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
}
