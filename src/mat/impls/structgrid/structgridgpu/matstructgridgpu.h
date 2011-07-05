

#ifndef __SEQSTRUCTGRIDGPU_H__
#define __SEQSTRUCTGRIDGPU_H__

//#include "../src/mat/impls/structgrid/matstructgrid.h"

#include "private/matimpl.h"





EXTERN_C_BEGIN
extern PetscErrorCode MatMult_SeqSGGPU(Mat, Vec, Vec);
extern PetscErrorCode MatCreate_SeqSGGPU(Mat);
extern PetscErrorCode MatCreate_SeqSG(Mat);
EXTERN_C_END
extern PetscErrorCode SGCUDA_MatMult( PetscScalar* coeff, PetscScalar* x, PetscScalar* y, PetscInt* idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos);


extern PetscErrorCode SGCUDA_MatMult_v2(PetscScalar* A, PetscScalar* B, PetscScalar* X, struct Stencilparams P);
//EXTERN_C_END
#endif
