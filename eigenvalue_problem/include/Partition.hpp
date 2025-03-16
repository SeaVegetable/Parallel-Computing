#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include "BSplineBasis.hpp"

class Partition
{
    public:
        Partition(const int &part_num_1d, const int &dim, const std::string &base_name)
            : part_num_1d(part_num_1d), dim(dim), base_name(base_name) {};
        
        ~Partition();

        void GeneratePartition(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2,
            const vector<double> &CP, const vector<int> &IEN, const vector<int> &ID,
            const vector<int> &NURBSExtraction1, const vector<int> &NURBSExtraction2);
        
        void WritePartition(const std::string &filename, const std::vector<int> local_to_global,
            const std::vector<double> &CP, const std::vector<int> &ID, const std::vector<int> &IEN,
            const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2) const;

        std::string GetPartitionFilename(const std::string &base_name, const int &rank) const;

    private:
        const int part_num_1d;
        const int dim;
};

#endif