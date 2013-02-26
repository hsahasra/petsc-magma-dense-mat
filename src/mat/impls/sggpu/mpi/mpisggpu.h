#ifndef __MPI_SGGPU_H__
#define __MPI_SGGPU_H__

#include "petsc-private/matimpl.h"

#include <../src/mat/impls/sggpu/seq/sggpu.h>


// External decls
EXTERN_C_BEGIN
extern PetscErrorCode MatCreate_MPISGGPU(Mat);
extern PetscErrorCode checkCudaError(cudaError_t err);
EXTERN_C_END



typedef struct {

  PetscInt nvec; 
  //Mat		A;			// Column major storage of diagonals */
  Mat_SeqSGGPU  *mat_seq;

  PetscMPIInt   size;                   /* size of communicator */
  PetscMPIInt   rank; 			/* rank of proc in communicator */

  /* The following variables are used for matrix-vector products */
  Vec           lvec;              	/* local vector */
  VecScatter    Mvctx;             	/* scatter context for vector */
  PetscBool     roworiented;       	/* if true, row-oriented input, default true */


} Mat_MPISGGPU;


// Prototypes
PetscErrorCode MatDestroy_MPISGGPU(Mat A);
PetscErrorCode MatSetGrid_MPISGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p);
PetscErrorCode MatMult_MPISGGPU(Mat mat, Vec x, Vec y);
PetscErrorCode MatSetValuesBlocked_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetValues_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetStencil_MPISGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof);
PetscErrorCode MatSetUp_MPISGGPU(Mat mat);
PetscErrorCode MatZeroEntries_MPISGGPU(Mat A);
PetscErrorCode MatGetDiagonal_MPISGGPU(Mat A, Vec v);
PetscErrorCode MatDiagonalScale_MPISGGPU(Mat A, Vec ll, Vec rr);
PetscErrorCode MatGetRow_MPISGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v);
PetscErrorCode MatRestoreRow_MPISGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v);
PetscErrorCode MatGetRowMaxAbs_MPISGGPU(Mat A, Vec v, PetscInt idx[]);
PetscErrorCode MatView_MPISGGPU(Mat A, PetscViewer viewer);
PetscErrorCode MatAssemblyBegin_MPISGGPU(Mat A, MatAssemblyType type);
PetscErrorCode MatAssemblyEnd_MPISGGPU(Mat A, MatAssemblyType type);
PetscErrorCode MatView_MPISGGPU_ASCII(Mat,PetscViewer);
extern PetscErrorCode MatFDColoringApply_MPISGGPU(Mat,MatFDColoring,Vec,MatStructure*,void*);
extern PetscErrorCode MatFDColoringCreate_MPISGGPU(Mat,ISColoring,MatFDColoring);
extern PetscErrorCode MatGetColumnIJ_MPISGGPU(Mat,PetscInt,PetscBool,PetscBool,PetscInt *,const PetscInt *[],const PetscInt *[],PetscBool*);
extern PetscErrorCode MatSetUpMultiply_MPISGGPU(Mat);
EXTERN_C_BEGIN
extern PetscErrorCode MatMPISGGPUSetPreallocation_MPISGGPU(Mat A,PetscInt nz, const PetscInt nnz[]);
EXTERN_C_END

#endif
