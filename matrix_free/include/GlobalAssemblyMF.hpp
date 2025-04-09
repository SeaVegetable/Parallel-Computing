#ifndef GLOBALASSEMBLYMF_HPP
#define GLOBALASSEMBLYMF_HPP

#include <petscmat.h>
#include "LocalAssemblyMF.hpp"

class GlobalAssemblyMF
{
    public:
        Vec F;

        GlobalAssemblyMF(const int &nLocBas, const int &nlocalfunc,
            const int &nlocalelemx, const int &nlocalelemy,
            const std::vector<int> &ghostID);
            
        ~GlobalAssemblyMF()
        {
            VecDestroy(&F);
        }

        void AssemLoad(LocalAssemblyMF * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<double> &CP,
            const std::vector<double> &NURBSExtraction1,
            const std::vector<double> &NURBSExtraction2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            ElementMF * const &elemmf);
        
        void MatMulMF(LocalAssemblyMF * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<double> &CP,
            const std::vector<double> &NURBSExtraction1,
            const std::vector<double> &NURBSExtraction2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            ElementMF * const &elemmf,
            Vec x, Vec y);
    
    private:
        const int nLocBas;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;
};

#endif