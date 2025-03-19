#ifndef GLOBALASSEMBLY_HPP
#define GLOBALASSEMBLY_HPP

#include <petscksp.h>
#include "LocalAssembly.hpp"

class GlobalAssembly
{
    public:
        Mat K;
        Vec F;

        GlobalAssembly(const int &nlocalfunc, const int &nLocBas);

        ~GlobalAssembly();

        void NonZeroCount(const Mat &K, std::vector<int> &dnz, std::vector<int> &onz);

        GetNonZeroEstimate();
    
    private:
        const int nlocalfunc;
};

#endif