#include "GlobalAssemblyMF.cuh"

GlobalAssemblyMF::GlobalAssemblyMF(const int &nLocBas, const int &nlocalfunc,
    const int &nlocalelemx, const int &nlocalelemy)
    : nLocBas(nLocBas), nlocalfunc(nlocalfunc),
      nlocalelemx(nlocalelemx), nlocalelemy(nlocalelemy)
{   
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, nlocalfunc, PETSC_DETERMINE);
    VecSetType(F, VECCUDA);
    VecSet(F, 0.0);
}

void GlobalAssemblyMF::AssemLoad(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf,
    BernsteinBasis * const &bernstein)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int nLocBas = elemmf->GetNumLocalBasis();
    const int pp = elemmf->get_p();
    const int qq = elemmf->get_q();

    std::vector<double> B1{};
    std::vector<double> B2{};
    std::vector<double> dB1{};
    std::vector<double> dB2{};

    for (int i = 0; i < nqp1; ++i)
    {
        B1.insert(B1.end(), bernstein->GetBernsteinBasisSingleQP(qp1[i]).begin(),
                          bernstein->GetBernsteinBasisSingleQP(qp1[i]).end());
        dB1.insert(dB1.end(), bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]).begin(),
                          bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]).end());
    }
    for (int j = 0; j < nqp2; ++j)
    {
        B2.insert(B2.end(), bernstein->GetBernsteinBasisSingleQP(qp2[j]).begin(),
                          bernstein->GetBernsteinBasisSingleQP(qp2[j]).end());
        dB2.insert(dB2.end(), bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]).begin(),
                          bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]).end());
    }

    double * d_B1, * d_B2, * d_dB1, * d_dB2;
    MallocDeviceMemory(&d_B1, B1.size());
    MallocDeviceMemory(&d_B2, B2.size());
    MallocDeviceMemory(&d_dB1, dB1.size());
    MallocDeviceMemory(&d_dB2, dB2.size());
    CopyToDeviceMemory(d_B1, B1.data(), B1.size());
    CopyToDeviceMemory(d_B2, B2.data(), B2.size());
    CopyToDeviceMemory(d_dB1, dB1.data(), dB1.size());
    CopyToDeviceMemory(d_dB2, dB2.data(), dB2.size());

    double * d_NURBSExtraction1, * d_NURBSExtraction2;
    MallocDeviceMemory(&d_NURBSExtraction1, NURBSExtraction1.size());
    MallocDeviceMemory(&d_NURBSExtraction2, NURBSExtraction2.size());
    CopyToDeviceMemory(d_NURBSExtraction1, NURBSExtraction1.data(), NURBSExtraction1.size());
    CopyToDeviceMemory(d_NURBSExtraction2, NURBSExtraction2.data(), NURBSExtraction2.size());

    int * d_IEN, * d_ID;
    double * d_CP;
    MallocDeviceMemory(&d_IEN, IEN.size());
    MallocDeviceMemory(&d_ID, ID.size());
    MallocDeviceMemory(&d_CP, CP.size());
    CopyToDeviceMemory(d_IEN, IEN.data(), IEN.size());
    CopyToDeviceMemory(d_ID, ID.data(), ID.size());
    CopyToDeviceMemory(d_CP, CP.data(), CP.size());

    double * d_elem_size1, * d_elem_size2;
    MallocDeviceMemory(&d_elem_size1, elem_size1.size());
    MallocDeviceMemory(&d_elem_size2, elem_size2.size());
    CopyToDeviceMemory(d_elem_size1, elem_size1.data(), elem_size1.size());
    CopyToDeviceMemory(d_elem_size2, elem_size2.data(), elem_size2.size());

    double * qw1, * qw2;
    MallocDeviceMemory(&qw1, nqp1);
    MallocDeviceMemory(&qw2, nqp2);
    CopyToDeviceMemory(qw1, w1.data(), nqp1);
    CopyToDeviceMemory(qw2, w2.data(), nqp2);

    VecSet(F, 0.0);
    double *d_F_array;
    VecCUDAGetArray(F, &d_F_array);

    AssembleLoadCUDA(pp, qq,
        nlocalelemx, nlocalelemy,
        d_B1, d_B2, d_dB1, d_dB2,
        d_NURBSExtraction1, d_NURBSExtraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2, d_F_array);
    
    DirichletBCCUDA(Dir.data(), static_cast<int>(Dir.size()), d_F_array, 0.0);

    VecCUDARestoreArray(F, &d_F_array);

    FreeDeviceMemory(d_B1);
    FreeDeviceMemory(d_B2);
    FreeDeviceMemory(d_dB1);
    FreeDeviceMemory(d_dB2);
    FreeDeviceMemory(d_NURBSExtraction1);
    FreeDeviceMemory(d_NURBSExtraction2);
    FreeDeviceMemory(d_IEN);
    FreeDeviceMemory(d_ID);
    FreeDeviceMemory(d_CP);
    FreeDeviceMemory(d_elem_size1);
    FreeDeviceMemory(d_elem_size2);
    FreeDeviceMemory(qw1);
    FreeDeviceMemory(qw2);
}

void GlobalAssemblyMF::MatMulMF(QuadraturePoint * const &quad1,
    QuadraturePoint * const &quad2,
    const std::vector<int> &IEN,
    const std::vector<int> &ID,
    const std::vector<int> &Dir,
    const std::vector<double> &CP,
    const std::vector<double> &NURBSExtraction1,
    const std::vector<double> &NURBSExtraction2,
    const std::vector<double> &elem_size1,
    const std::vector<double> &elem_size2,
    ElementMF * const &elemmf,
    BernsteinBasis * const &bernstein,
    Vec x, Vec y)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int nLocBas = elemmf->GetNumLocalBasis();
    const int pp = elemmf->get_p();
    const int qq = elemmf->get_q();

    std::vector<double> B1{};
    std::vector<double> B2{};
    std::vector<double> dB1{};
    std::vector<double> dB2{};

    for (int i = 0; i < nqp1; ++i)
    {
        B1.insert(B1.end(), bernstein->GetBernsteinBasisSingleQP(qp1[i]).begin(),
                          bernstein->GetBernsteinBasisSingleQP(qp1[i]).end());
        dB1.insert(dB1.end(), bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]).begin(),
                          bernstein->GetBernsteinBasisDerivativeSingleQP(qp1[i]).end());
    }
    for (int j = 0; j < nqp2; ++j)
    {
        B2.insert(B2.end(), bernstein->GetBernsteinBasisSingleQP(qp2[j]).begin(),
                          bernstein->GetBernsteinBasisSingleQP(qp2[j]).end());
        dB2.insert(dB2.end(), bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]).begin(),
                          bernstein->GetBernsteinBasisDerivativeSingleQP(qp2[j]).end());
    }

    double * d_B1, * d_B2, * d_dB1, * d_dB2;
    MallocDeviceMemory(&d_B1, B1.size());
    MallocDeviceMemory(&d_B2, B2.size());
    MallocDeviceMemory(&d_dB1, dB1.size());
    MallocDeviceMemory(&d_dB2, dB2.size());
    CopyToDeviceMemory(d_B1, B1.data(), B1.size());
    CopyToDeviceMemory(d_B2, B2.data(), B2.size());
    CopyToDeviceMemory(d_dB1, dB1.data(), dB1.size());
    CopyToDeviceMemory(d_dB2, dB2.data(), dB2.size());

    double * d_NURBSExtraction1, * d_NURBSExtraction2;
    MallocDeviceMemory(&d_NURBSExtraction1, NURBSExtraction1.size());
    MallocDeviceMemory(&d_NURBSExtraction2, NURBSExtraction2.size());
    CopyToDeviceMemory(d_NURBSExtraction1, NURBSExtraction1.data(),
        NURBSExtraction1.size());
    CopyToDeviceMemory(d_NURBSExtraction2, NURBSExtraction2.data(),
        NURBSExtraction2.size());
    
    int * d_IEN, * d_ID;
    double * d_CP;
    MallocDeviceMemory(&d_IEN, IEN.size());
    MallocDeviceMemory(&d_ID, ID.size());
    MallocDeviceMemory(&d_CP, CP.size());
    CopyToDeviceMemory(d_IEN, IEN.data(), IEN.size());
    CopyToDeviceMemory(d_ID, ID.data(), ID.size());
    CopyToDeviceMemory(d_CP, CP.data(), CP.size());

    double * d_elem_size1, * d_elem_size2;
    MallocDeviceMemory(&d_elem_size1, elem_size1.size());
    MallocDeviceMemory(&d_elem_size2, elem_size2.size());
    CopyToDeviceMemory(d_elem_size1, elem_size1.data(), elem_size1.size());
    CopyToDeviceMemory(d_elem_size2, elem_size2.data(), elem_size2.size());

    double * qw1, * qw2;
    MallocDeviceMemory(&qw1, nqp1.size());
    MallocDeviceMemory(&qw2, nqp2.size());
    CopyToDeviceMemory(qw1, w1.data(), nqp1.size());
    CopyToDeviceMemory(qw2, w2.data(), nqp2.size());

    VecSet(y, 0.0);
    double *d_x_array, *d_y_array;
    VecCUDAGetArray(x, &d_x_array);
    VecCUDAGetArray(y, &d_y_array);

    MatrixFreeMatMultCUDA(pp, qq,
        nlocalelemx, nlocalelemy,
        d_B1, d_B2, d_dB1, d_dB2,
        d_NURBSExtraction1, d_NURBSExtraction2,
        d_elem_size1, d_elem_size2,
        d_IEN, d_ID, d_CP,
        qw1, qw2, d_x_array, d_y_array);

    DirichletBCCUDA(Dir.data(), static_cast<int>(Dir.size()), d_y_array, 0.0);
    
    VecCUDARestoreArray(x, &d_x_array);
    VecCUDARestoreArray(y, &d_y_array);

    FreeDeviceMemory(d_B1);
    FreeDeviceMemory(d_B2);
    FreeDeviceMemory(d_dB1);
    FreeDeviceMemory(d_dB2);
    FreeDeviceMemory(d_NURBSExtraction1);
    FreeDeviceMemory(d_NURBSExtraction2);
    FreeDeviceMemory(d_IEN);
    FreeDeviceMemory(d_ID);
    FreeDeviceMemory(d_CP);
    FreeDeviceMemory(d_elem_size1);
    FreeDeviceMemory(d_elem_size2);
    FreeDeviceMemory(qw1);
    FreeDeviceMemory(qw2);
}