#include "InvLM.hpp"

InvLM::InvLM(const int &nLocBas,
    const int &nLocElem,
    const std::vector<int> local_to_global_wo_dir,
    const std::vector<int> ID,
    const std::vector<int> IEN)
{
    const int nlocfunc = local_to_global_wo_dir.size();

    elemNum.resize(nlocfunc, 0);
    std::vector<int> size_count(nlocfunc, 0);

    for (int i = 0; i < nLocElem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            const int global_func = ID[IEN[i*nLocBas + j]];
            for (int k = 0; k < nlocfunc; ++k)
            {
                if (local_to_global_wo_dir[k] == global_func)
                {
                    elemNum[k] += 1;
                }
            }
        }
    }

    int size = 0;
    for (int i = 0; i < nlocfunc; ++i)
         size += elemNum[i];
    
    elemIdx.resize(size, -1);
    baseIdx.resize(size, -1);

    for (int i = 0; i < nLocElem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            const int global_func = ID[IEN[i*nLocBas + j]];
            for (int k = 0; k < nlocfunc; ++k)
            {
                if (local_to_global_wo_dir[k] == global_func)
                {
                    elemIdx[size_count[k] + (k == 0 ? 0 : std::accumulate(elemNum.begin(), elemNum.begin() + k, 0))] = i;
                    baseIdx[size_count[k] + (k == 0 ? 0 : std::accumulate(elemNum.begin(), elemNum.begin() + k, 0))] = j;
                    size_count[k] += 1;
                }
            }
        }
    }
}