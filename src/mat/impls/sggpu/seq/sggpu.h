
#ifndef __SEQ_SGGPU_H__
#define __SEQ_SGGPU_H__

#include "petsc-private/matimpl.h"

#include <map>
#include <vector>
#include <cuda.h>

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

  PetscInt      non_zeros;  //< Count of non-zero entries

  PetscScalar * deviceData; //< Device data

  std::map<int, int> * diag_starts;
  std::vector<int> * diagonals;

  cudaStream_t stream;
} Mat_SeqSGGPU;

#endif
