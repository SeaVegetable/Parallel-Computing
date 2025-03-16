#include <iostream>
#include "ControlPointGenerator.hpp"
#include "IENGenerator.hpp"
#include "IDGenerator.hpp"
#include "NURBSExtractionGenerator.hpp"
#include "Partition.hpp"

int main(int argc, char *argv[])
{
    int p = 2;
    int q = 2;

    double Lx = 1.0;
    double Ly = 1.0;

    int nElemX = 3;
    int nElemY = 3;

    std::string file_info = "info.txt";
    std::ofstream file(file_info.c_str());
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << file_info << std::endl;
        exit(1);
    }

    file << "p: " << p << std::endl;
    file << "q: " << q << std::endl;
    file << "Lx: " << Lx << std::endl;
    file << "Ly: " << Ly << std::endl;
    file << "nElemX: " << nElemX << std::endl;
    file << "nElemY: " << nElemY << std::endl;

    file.close();

    double hx = Lx / nElemX;
    double hy = Ly / nElemY;

    int nElem = nElemX * nElemY;
    int nLocBas = (p + 1) * (q + 1);
    
    int nFuncX = p + nElemX;
    int nFuncY = q + nElemY;
    int nFunc = nFuncX * nFuncY;

    std::vector<double> S;
    std::vector<double> T;

    for (int i = 0; i < p + 1; ++i) S.push_back(0.0);
    for (int i = 1; i < nElemX; ++i) S.push_back(i * hx);
    for (int i = 0; i < p + 1; ++i) S.push_back(Lx);

    for (int i = 0; i < q + 1; ++i) T.push_back(0.0);
    for (int i = 1; i < nElemY; ++i) T.push_back(i * hy);
    for (int i = 0; i < q + 1; ++i) T.push_back(Ly);

    ControlPointGenerator * cpg = new ControlPointGenerator();
    BSplineBasis * basis1 = new BSplineBasis(p, S);
    BSplineBasis * basis2 = new BSplineBasis(q, T);

    double P1min = 0.0;
    double P1max = Lx;
    double P2min = 0.0;
    double P2max = Ly;

    std::vector<double> CP = cpg->GenerateControlPoints2D(basis1, basis2, 0.0, Lx, 0.0, Ly);

    IENGenerator * igen = new IENGenerator();
    std::vector<int> IEN = igen->GenerateIEN2D(basis1, basis2);

    IDGenerator * idgen = new IDGenerator();
    std::vector<int> ID = idgen->GenerateID2D(basis1, basis2);

    NURBSExtractionGenerator * neg = new NURBSExtractionGenerator();
    std::vector<int> NURBSExtraction1 = neg->GenerateNURBSExtraction(basis1);
    std::vector<int> NURBSExtraction2 = neg->GenerateNURBSExtraction(basis2);

    std::string base_name = "part";
    Partition * part = new Partition(4, 2, base_name);
    part->GeneratePartition(basis1, basis2, CP, IEN, ID, NURBSExtraction1, NURBSExtraction2);

    return 0;
}