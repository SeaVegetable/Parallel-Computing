#include "Partition.hpp"

void Partition::GeneratePartition(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2,
    const std::vector<double> &CP, const std::vector<int> &IEN, const std::vector<int> &ID,
    const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2)
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

    FileManager * fm = new FileManager();

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
            std::cout << "Generating partition " << count << "..." << std::endl;
            const int localnElem = localnElemX * localnElemY;
            const int extSizeX = (p + 1) * (p + 1);
            const int extSizeY = (q + 1) * (q + 1);

            std::vector<double> localCP{};
            std::vector<int> local_to_global{};
            std::vector<int> localIEN{};
            std::vector<int> localID{};
            std::vector<double> localNURBSExtraction1{};
            std::vector<double> localNURBSExtraction2{};

            for (int localj = 0; localj < localnElemY; ++localj)
            {
                for (int locali = 0; locali < localnElemX; ++locali)
                {
                    for (int ii = 0; ii < nLocBas; ++ii)
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
                    for (int ii = 0; ii < nLocBas; ++ii)
                    {
                        int index = i + locali + (j + localj) * nElemX;
                        localIEN.push_back(std::distance(local_to_global.begin(),
                            std::find(local_to_global.begin(), local_to_global.end(), IEN[index * nLocBas + ii])));
                    }
                }
            }

            for (int locali = 0; locali < localnElemX; ++locali)
            {
                size_t start = (i + locali) * extSizeX;
                size_t end = (i + locali + 1) * extSizeX;
                    
                std::copy(NURBSExtraction1.begin() + start, NURBSExtraction1.begin() + end,
                    std::back_inserter(localNURBSExtraction1));
            }

            for (int localj = 0; localj < localnElemY; ++localj)
            {
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
                localID.push_back(ID[local_to_global[ii]]);
            }

            std::string filename = fm->GetPartitionFilename(base_name, count);
            fm->WritePartition(filename, local_to_global, localCP, localID,
                localIEN, localNURBSExtraction1, localNURBSExtraction2);

            std::cout << "Partition " << count << " generated." << std::endl;
            ++count;
        }
    }

    delete fm; fm = nullptr;
}