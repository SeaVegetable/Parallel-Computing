#include "FileManager.hpp"

int main(int argc, char **argv)
{
    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    std::vector<int> ownlocalID;
    std::vector<int> ownghostID;
    std::vector<int> otherlocalID;
    std::vector<int> otherghostID;

    const int total_ranks = part_num_1d * part_num_1d;
    for (int rank = 0; rank < total_ranks; ++rank)
    {
        std::string filename = fm->GetPartitionFilename(base_name, rank);
        fm->ReadPartition(filename, ownlocalID, ownghostID);
        for (int other_rank = 0; other_rank < total_ranks; ++other_rank)
        {
            if (other_rank != rank)
            {
                std::string other_filename = fm->GetPartitionFilename(base_name, other_rank);
                fm->ReadPartition(other_filename, otherlocalID, otherghostID);
            }
        }
    }
}