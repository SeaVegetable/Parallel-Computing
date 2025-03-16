#include "Partition.hpp"

void GeneratePartition(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2,
    const vector<double> &CP, const vector<int> &IEN, const vector<int> &ID,
    const vector<double> &NURBSExtraction1, const vector<double> &NURBSExtraction2)
{
    const int p = basis1->GetDegree();
    const int q = basis2->GetDegree();
    const std::vector<double> S = basis1->GetKnotVector();
    const std::vector<double> T = basis2->GetKnotVector();
    const int m = static_cast<int>(S.size()) - p - 1;
    const int n = static_cast<int>(T.size()) - q - 1;
    const int nElemX = m - p;
    const int nElemY = n - q;
    const int nLocBas = (p + 1) * (q + 1);

    const int localnElemX = nElemX / part_num_1d;
    const int localnElemY = nElemY / part_num_1d;

    if (nElemX % part_num_1d != 0 || nElemY % part_num_1d != 0)
    {
        std::cerr << "Error: Number of partitions must divide number of elements in each direction." << std::endl;
        exit(1);
    }

    int count = 0;

    for (int j = 0; j < nElemY; j += localnElemY)
    {
        for (int i = 0; i < nElemX; i += localnElemX)
        {
            const int localnElem = localnElemX * localnElemY;
            const int extSizeX = (p + 1) * (p + 1);
            const int extSizeY = (q + 1) * (q + 1);

            std::vector<double> localCP{};
            std::vector<int> local_to_global{};
            std::vector<int> localIEN{};
            std::vector<double> localNURBSExtraction1{};
            std::vector<double> localNURBSExtraction2{};

            for (int localj = 0; localj < localnElemY; ++localj)
            {
                for (int locali = 0; locali < localnElemX; ++locali)
                {
                    for (int ii = 0; i < nLocBas; ++ii)
                    {
                        int index = i + locali + (j + localj) * nElemX;
                        local_to_global.push_back(IEN[index * nLocBas + ii]);
                    }
                }
            }

            std::sort(local_to_global.begin(), local_to_global.end());
            local_to_global.erase(std::unique(local_to_global.begin(), local_to_global.end()), local_to_global.end());

            for (int localj = 0; localj < localnElemY; ++localj)
            {
                for (int locali = 0; locali < localnElemX; ++locali)
                {
                    for (int ii = 0; i < nLocBas; ++ii)
                    {
                        int index = i + locali + (j + localj) * nElemX;
                        localIEN.push_back(std::distance(local_to_global.begin(),
                            std::find(local_to_global.begin(), local_to_global.end(), IEN[index * nLocBas + ii])));
                    }

                    size_t start = (i + locali) * extSizeX;
                    size_t end = (i + locali + 1) * extSizeX;
                    
                    std::copy(NURBSExtraction1.begin() + start, NURBSExtraction1.begin() + end,
                        std::back_inserter(localNURBSExtraction1));
                }

                size_t start = (j + localj) * extSizeY;
                size_t end = (j + localj + 1) * extSizeY;

                std::copy(NURBSExtraction2.begin() + start, NURBSExtraction2.begin() + end,
                    std::back_inserter(localNURBSExtraction2));
            }

            for (int ii = 0; ii < local_to_global.size(); ++ii)
            {
                for (int jj = 0; jj < dim; ++jj)
                {
                    localCP.push_back(CP[local_to_global[ii] * dim + jj]);
                }
            }

            std::string filename = GetPartitionFilename("Partition", count);
            WritePartition(filename, local_to_global, localCP,
                localIEN, localNURBSExtraction1, localNURBSExtraction2);
            ++count;
        }
    }
}

void WritePartition(const std::string &filename, const std::vector<std::int> &local_to_global,
    const std::vector<double> &CP, const std::vector<int> &IEN,
    const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2) const
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    file << "LocalToGlobal" << std::endl;
    for (int i = 0; i < local_to_global.size(); ++i)
    {
        file << local_to_global[i] << std::endl;
    }

    file << "CP" << std::endl;
    for (int i = 0; i < CP.size(); ++i)
    {
        file << CP[i] << std::endl;
    }

    file << "IEN" << std::endl;
    for (int i = 0; i < IEN.size(); ++i)
    {
        file << IEN[i] << std::endl;
    }

    file << "NURBSExtraction1" << std::endl;
    for (int i = 0; i < NURBSExtraction1.size(); ++i)
    {
        file << NURBSExtraction1[i] << std::endl;
    }

    file << "NURBSExtraction2" << std::endl;
    for (int i = 0; i < NURBSExtraction2.size(); ++i)
    {
        file << NURBSExtraction2[i] << std::endl;
    }

    file.close();
}

std::string GetPartitionFilename(const std::string &base_name, const int &rank) const
{
    std::string filename = base_name + "_" + std::to_string(rank) + ".txt";
    return filename;
}