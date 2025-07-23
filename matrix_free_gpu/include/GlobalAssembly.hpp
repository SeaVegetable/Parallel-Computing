#ifndef GLOBALASSEMBLY_HPP
#define GLOBALASSEMBLY_HPP

#include <petscmat.h>
#include "LocalAssembly.hpp"

class GlobalAssembly
{
    public:
        Mat K;
        Vec F;

        GlobalAssembly(const std::vector<int> &IEN, const std::vector<int> &ID,
            const std::vector<int> &Dir, LocalAssembly * const &locassem, const int &nLocBas,
            const int &nlocalfunc, const int &nlocalelemx, const int &nlocalelemy);

        ~GlobalAssembly();

        void NonZeroCount(const Mat &K, std::vector<int> &dnz, std::vector<int> &onz);

        void AssemNonZeroEstimate(LocalAssembly * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<int> &Dir);

        void AssemStiffnessLoad(LocalAssembly * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<int> &Dir,
            const std::vector<double> &CP,
            const std::vector<double> &NURBSExtraction1,
            const std::vector<double> &NURBSExtraction2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            Element * const &elem);
        
        void AssemStiffnessLoad(LocalAssembly * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<int> &Dir,
            const std::vector<double> &CP,
            ElementFEM * const &elem);
    
    private:
        const int nLocBas;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;

        void DirichletBCK(const std::vector<int> &Dir);

        void DirichletBC(const std::vector<int> &Dir);
};

#endif