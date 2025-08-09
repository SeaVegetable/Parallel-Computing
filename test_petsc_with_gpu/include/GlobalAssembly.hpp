#ifndef GLOBALASSEMBLY_HPP
#define GLOBALASSEMBLY_HPP

#include <petscmat.h>
#include "QuadraturePoint.hpp"

class GlobalAssembly
{
    public:
        Mat K;

        GlobalAssembly(const int &nLocBas,
            const int &nnz, const int &nlocalfunc,
            const int &nlocalelemx, const int &nlocalelemy,
            const std::vector<int> &rows,
            const std::vector<int> &cols);

        ~GlobalAssembly();

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
        
        void AssemStiffness(QuadraturePoint * const &quad1,
            QuadraturePoint * const &quad2,
            const std::vector<int> &IEN,
            const std::vector<int> &ID,
            const std::vector<double> &CP,
            
            ElementFEM * const &elem);
    
    private:
        const int nLocBas;
        const int nnz;
        const int nlocalfunc;
        const int nlocalelemx;
        const int nlocalelemy;

        void DirichletBCK(const std::vector<int> &Dir);

        void DirichletBC(const std::vector<int> &Dir);
};

#endif