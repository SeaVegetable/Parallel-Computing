#include "InvLM.hpp"

InvLM::InvLM(const int &nLocBas,
    const int &nLocElem,
    const int &nlocfunc,
    const std::vector<int> IEN)
{
    elemNum.resize(nlocfunc, 0);
    offset.resize(nlocfunc, 0);
    std::vector<int> size_count(nlocfunc, 0);

    for (int i = 0; i < nLocElem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            const int global_func = IEN[i*nLocBas + j];
            for (int k = 0; k < nlocfunc; ++k)
            {
                if ( k == global_func )
                {
                    elemNum[k] += 1;
                }
            }
        }
    }

    int size = 0;
    for (int i = 0; i < nlocfunc; ++i)
    {
        offset[i] = size;
        size += elemNum[i];
    }

    elemIdx.resize(size, -1);
    baseIdx.resize(size, -1);

    for (int i = 0; i < nLocElem; ++i)
    {
        for (int j = 0; j < nLocBas; ++j)
        {
            const int global_func = IEN[i*nLocBas + j];
            for (int k = 0; k < nlocfunc; ++k)
            {
                if ( k == global_func )
                {
                    elemIdx[size_count[k] + offset[k]] = i;
                    baseIdx[size_count[k] + offset[k]] = j;
                    size_count[k] += 1;
                }

                if ( size_count[k] == elemNum[k] )
                {
                    break;
                }
            }
        }
    }
}