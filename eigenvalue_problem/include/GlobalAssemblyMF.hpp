#ifndef GLOBALASSEMBLYMF_HPP
#define GLOBALASSEMBLYMF_HPP

#include <petscmat.h>
#include "LocalAssemblyMF.hpp"

class GlobalAssemblyMF
{
    public:
        Vec F;

        GlobalAssemblyMF(const std::vector<int> &IEN, const std::vector<int> &ID,
            LocalAssembly * const &locassem, const int &nLocBas,
            const int &nlocalfunc, const int &nlocalelemx, const int &nlocalelemy);

        ~GlobalAssemblyMF();

        void AssemLoad(LocalAssembly * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<double> &CP,
            const std::vector<double> &NURBSExtraction1,
            const std::vector<double> &NURBSExtraction2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            Element * const &elem);
        
        void MatMulMF(LocalAssembly * const &locassem,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<double> &CP,
            const std::vector<double> &NURBSExtraction1,
            const std::vector<double> &NURBSExtraction2,
            const std::vector<double> &elem_size1,
            const std::vector<double> &elem_size2,
            Element * const &elem,
            Vec x, Vec y);
    
    private:
        const int nLocBas;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;
};

#endif