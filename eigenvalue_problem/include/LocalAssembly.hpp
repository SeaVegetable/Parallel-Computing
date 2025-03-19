#ifndef LOCALASSEMBLY_HPP
#define LOCALASSEMBLY_HPP

#include <petscmat.h>
#include "Element.hpp"
#include "QuadraturePoint.hpp"

class LocalAssembly
{
    public:
        LocalAssembly(const double &p, const double &q)
        {
            quad1 = new QuadraturePoint(p+1, 0, 1);
            quad2 = new QuadraturePoint(q+1, 0, 1);
            const int n = (p+1)*(q+1);
            Kloc = new PetscScalar[n*n];
            Floc = new PetscScalar[n];
        }

        ~LocalAssembly()
        {
            delete quad1; quad1 = nullptr;
            delete quad2; quad2 = nullptr;
            delete[] Kloc; Kloc = nullptr;
            delete[] Floc; Floc = nullptr;
        }

        void AssemLocalStiffness(const Element * const &elem,
            const std::vector<double> &eCP);

    private:
        QuadraturePoint * quad1;
        QuadraturePoint * quad2;
        PetscScalar *Kloc;
        PetscScalar *Floc;
};

#endif