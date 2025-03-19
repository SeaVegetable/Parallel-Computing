#ifndef LOCALASSEMBLY_HPP
#define LOCALASSEMBLY_HPP

class LocalAssembly
{
    public:
        LocalAssembly();
        ~LocalAssembly();

    private:
        PetscScalar *Kloc;
        PetscScalar *Floc;
};

#endif