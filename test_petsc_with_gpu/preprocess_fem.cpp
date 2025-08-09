// #if __has_include(<filesystem>)
//   #include <filesystem>
//   namespace fs = std::__fs::filesystem;
// #elif __has_include(<experimental/filesystem>)
//   #include <experimental/filesystem>
//   namespace fs = std::experimental::filesystem;
// #else
//   #error "No filesystem support"
// #endif
#include "AbscissaeGenerator.hpp"
#include "IENGenerator.hpp"
#include "IDGenerator.hpp"
#include "Partition.hpp"

int main(int argc, char *argv[])
{
    int p, q, nElemX, nElemY, part_num_1d, dim;
    double Lx, Ly;
    std::string base_name;

    std::string file_info = "info.txt";

    FileManager * fm = new FileManager();
    fm->ReadPreprocessInfo(file_info, p, q, Lx, Ly, nElemX, nElemY, part_num_1d, dim, base_name);

    std::string base_name_fem = "part_fem";
    
    // for (const auto& entry : fs::directory_iterator(fs::current_path())) {
    //     if (entry.is_regular_file()) {
    //         std::string filename = entry.path().filename().string();
    //         if (filename.find(base_name_fem) == 0) {
    //             std::cout << "Deleting file: " << filename << std::endl;
    //             fs::remove(entry.path());
    //         }
    //     }
    // }

    double hx = Lx / (nElemX + 2);
    double hy = Ly / (nElemY + 2);

    int nElem = nElemX * nElemY;
    int nLocBas = (p + 1) * (q + 1);
    
    int nFuncX = p + nElemX;
    int nFuncY = q + nElemY;
    int nFunc = nFuncX * nFuncY;

    std::vector<double> S;
    std::vector<double> T;

    for (int i = 0; i < p - 1; ++i) S.push_back(0.0);
    for (int i = 1; i < nElemX+2; ++i) S.push_back(i * hx);
    for (int i = 0; i < p - 1; ++i) S.push_back(Lx);

    for (int i = 0; i < q - 1; ++i) T.push_back(0.0);
    for (int i = 1; i < nElemY+2; ++i) T.push_back(i * hy);
    for (int i = 0; i < q - 1; ++i) T.push_back(Ly);

    AbscissaeGenerator * absgen = new AbscissaeGenerator();
    std::vector<double> CP = absgen->GenerateAbscissae2D(S, T, p-2, q-2);

    nElemX = nFuncX - 1;
    nElemY = nFuncY - 1;
    nElem = nElemX * nElemY;

    IENGenerator * ien = new IENGenerator();
    std::vector<int> IEN = ien->GenerateIEN2D(nElemX, nElemY);

    IDGenerator * idgen = new IDGenerator();
    std::vector<int> ID = idgen->GenerateID2D(nFuncX, nFuncY);

    Partition * part = new Partition(part_num_1d, part_num_1d, dim, base_name_fem);
    part->GeneratePartition(nElemX, nElemY, CP, IEN, ID);

    delete fm; fm = nullptr;
    delete ien; ien = nullptr;
    delete idgen; idgen = nullptr;
    delete absgen; absgen = nullptr;
    delete part; part = nullptr;
}