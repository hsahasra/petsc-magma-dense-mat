
#ifndef __SEQ_SGGPU_H__
#define __SEQ_SGGPU_H__

#include "petsc-private/matimpl.h"

// External decls
EXTERN_C_BEGIN
extern PetscErrorCode MatCreate_SeqSGGPU(Mat);
EXTERN_C_END


// Matrix type container
typedef struct {
  PetscScalar * hostData;   //< Host data
  PetscInt      stpoints;   //< Number of stencil points
  PetscInt      dof;        //< Degrees of freedom
  PetscInt      m;          //< Grid size (x)
  PetscInt      n;          //< Grid size (y)
  PetscInt      p;          //< Grid size (z)

  PetscScalar * deviceData; //< Device data
} Mat_SeqSGGPU;

#endif
