#ifndef FILEMANAGER_HPP
#define FILEMANAGER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class FileManager
{
public:
    FileManager(){}
    ~FileManager();

    void WritePartition(const std::string &filename, const std::vector<int> &local_to_global,
        const std::vector<double> &CP, const std::vector<int> &ID, const std::vector<int> &IEN,
        const std::vector<double> &NURBSExtraction1, const std::vector<double> &NURBSExtraction2) const;
    
    void ReadPartition(const std::string &filename, std::vector<int> &local_to_global,
        std::vector<double> &CP, std::vector<int> &ID, std::vector<int> &IEN,
        std::vector<double> &NURBSExtraction1, std::vector<double> &NURBSExtraction2) const;
    
    void WritePreprocessInfo(const std::string &filename, const int &p, const int &q, const double &Lx, const double &Ly,
        const int &nElemX, const int &nElemY, const int &part_num_1d, const int &dim, const std::string &base_name) const;

    void ReadPreprocessInfo(const std::string &filename, int &p, int &q, double &Lx, double &Ly,
        int &nElemX, int &nElemY, int &part_num_1d, int &dim, std::string &base_name) const;

    std::string GetPartitionFilename(const std::string &base_name, const int &rank) const;
};

#endif