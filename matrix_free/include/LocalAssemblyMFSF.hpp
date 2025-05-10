#ifndef LOCALASSEMBLYMFSF_HPP
#define LOCALASSEMBLYMFSF_HPP

#include <petscmat.h>
#include "ElementMFSF.hpp"

class LocalAssemblyMFSF
{
    public:
        PetscScalar *Floc_in;
        PetscScalar *Floc_out;
        PetscScalar *Floc;
        PetscScalar *Kloc;

        LocalAssemblyMFSF(const double &p, const double &q)
            : n((p+1)*(q+1))
        {
            quad1 = new QuadraturePoint(p+1, 0, 1);
            quad2 = new QuadraturePoint(q+1, 0, 1);
            Floc = new PetscScalar[n];
            Kloc = new PetscScalar[n*n];
            Floc_in = new PetscScalar[n];
            Floc_out = new PetscScalar[n];
        }

        ~LocalAssemblyMFSF()
        {
            delete quad1; quad1 = nullptr;
            delete quad2; quad2 = nullptr;
            delete[] Floc; Floc = nullptr;
            delete[] Kloc; Kloc = nullptr;
            delete[] Floc_in; Floc_in = nullptr;
            delete[] Floc_out; Floc_out = nullptr;
        }

        void AssemLocalLoad(ElementMFSF * const &elem,
            const std::vector<double> &eCP);
        
        void LocalMatMulMF(ElementMFSF * const &elem,
            const std::vector<double> &eCP);

    private:
        const int n;
        QuadraturePoint * quad1;
        QuadraturePoint * quad2;

        double Getf(const double &xi, const double &eta)
        {
            return xi*(1-xi)*eta*(1-eta);
        }

        void ResetLoad()
        {
            for (int i = 0; i < n; ++i)
                Floc[i] = 0.0;
        }

        void ResetStiffnessLoadOut()
        {
            for (int i = 0; i < n*n; ++i)
                Kloc[i] = 0.0;
            for (int i = 0; i < n; ++i)
                Floc_out[i] = 0.0;
        }
};

#endif