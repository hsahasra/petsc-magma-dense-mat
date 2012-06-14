

#ifndef __SEQSTRUCTGRIDGPU_H__
#define __SEQSTRUCTGRIDGPU_H__

//#include "../src/mat/impls/structgrid/matstructgrid.h"

#include "petsc-private/matimpl.h"


EXTERN_C_BEGIN
extern PetscErrorCode MatMult_SeqSGGPU(Mat, Vec, Vec);
extern PetscErrorCode MatCreate_SeqSGGPU(Mat);
extern PetscErrorCode MatCreate_SeqSG(Mat);
extern PetscErrorCode MatDestroy_SeqSG(Mat);
extern PetscErrorCode MatDestroy_SeqSGGPU(Mat);
EXTERN_C_END
extern PetscErrorCode SGCUDA_MatMult( PetscScalar* coeff, PetscScalar* x, PetscScalar* y,
PetscInt* idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,
PetscInt p, PetscInt nos, PetscCUSPFlag* fp,PetscInt DOF, PetscScalar ** gpuMat);


extern PetscErrorCode SGCUDA_MatMult_v2(PetscScalar* A, PetscScalar* B,
PetscScalar* X, struct Stencilparams P, PetscCUSPFlag* fp);

//EXTERN_C_END
#endif
