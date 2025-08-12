#include "Elem2COOGenerator.hpp"

Elem2COOGenerator::GenerateElem2COO(const std::vector<int> &IEN,
    const std::vector<int> &ID, const std::vector<int> &rows,
    const std::vector<int> &cols, std::vector<int> &elem2coo)
{
    elem2coo.clear();
    elem2coo.resize(IEN.size() * nLocBas);

    #pragma omp parallel for
    for (int elem = 0; elem < nlocalelemx * nlocalelemy; ++elem)
    {
        for (int i = 0; i < nLocBas; ++i)
        {
            int II = ID[IEN[elem * nLocBas + i]];
            int base_idx = elem * nLocBas * nLocBas + i * nLocBas;
            if(II >= 0)
            {
                auto it = std::find(rows.begin(), rows.end(), II);
                int Index = std::distance(rows.begin(), it);

                for (int j = 0; j < nLocBas; ++j)
                {
                    int JJ = ID[IEN[elem * nLocBas + j]];
                    if(JJ >= 0)
                    {
                        auto it2 = std::find(cols.begin() + Index, cols.end(), JJ);
                        int colIndex = std::distance(cols.begin() + Index, it2);
                        elem2coo[base_idx + j] = Index + colIndex;
                    }
                    else
                    {
                        elem2coo[base_idx + j] = -1; // If JJ is negative, we do not add it to the COO format
                    }
                }
            }
            else
            {
                for (int j = 0; j < nLocBas; ++j)
                {
                    elem2coo[base_idx + j] = -1; // If II is negative, we do not add it to the COO format
                }
            }
        }
    }
}

void Elem2COOGenerator::GenerateDir2COO(const std::vector<int> &Dir,
    const std::vector<int> &rows, std::vector<int> &dir2coo)
{
    dir2coo.clear();
    dir2coo.reserve(Dir.size());
    for (size_t i = 0; i < Dir.size(); ++i)
    {
        auto it = std::find(rows.begin(), rows.end(), Dir[i]);
        dir2coo.push_back(std::distance(rows.begin(), it));
    }
}