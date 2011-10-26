

#ifndef __SEQSTRUCTGRIDGPU_H__
#define __SEQSTRUCTGRIDGPU_H__

//#include "../src/mat/impls/structgrid/matstructgrid.h"

#include "private/matimpl.h"
#include "private/vecimpl.h"
#include "petscconf.h"
#include "petscsys.h"
#include "petscvec.h"
#include "../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"
#include "../include/private/petscimpl.h"





EXTERN_C_BEGIN
extern PetscErrorCode MatMult_SeqSGGPU(Mat, Vec, Vec);
extern PetscErrorCode MatCreate_SeqSGGPU(Mat);
extern PetscErrorCode MatCreate_SeqSG(Mat);
extern PetscErrorCode MatDestroy_SeqSG(Mat);
extern PetscErrorCode MatDestroy_SeqSGGPU(Mat);
//extern PetscErrorCode SGCUDA_MatMult( PetscScalar* coeff, PetscScalar* x, PetscScalar* y, 
//PetscInt* idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,
//PetscInt p, PetscInt nos, PetscCUSPFlag* fp);

extern PetscErrorCode MatCheckCUDAError(const char *);
extern PetscErrorCode MatCheckCUDAStatus(cudaError_t,const char *);


PetscErrorCode SGCUDA_MatMult_v2(Mat mat, Vec x,Vec y);
EXTERN_C_END

#endif
