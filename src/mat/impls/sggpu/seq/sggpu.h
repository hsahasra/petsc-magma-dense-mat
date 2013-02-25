
#ifndef __SEQ_SGGPU_H__
#define __SEQ_SGGPU_H__

#include "petsc-private/matimpl.h"

#include <map>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>

// External decls
EXTERN_C_BEGIN
extern PetscErrorCode MatCreate_SeqSGGPU(Mat);
EXTERN_C_END


#if _TRACE
#define SGTrace printf("[SeqSGGPU] %s\n",__FUNCT__);
#else
#define SGTrace
#endif


// Matrix type container
typedef struct {
  PetscScalar * hostData;   //< Host data
  PetscInt      stpoints;   //< Number of stencil points
  PetscInt      dof;        //< Degrees of freedom
  PetscInt      m;          //< Grid size (x)
  PetscInt      n;          //< Grid size (y)
  PetscInt      p;          //< Grid size (z)
  PetscInt      dim;        //< Dimensionality
  PetscInt      stencil_type; /*0 for star, 1 for box */
  PetscInt      *ja, *ia;

  PetscInt      non_zeros;  //< Count of non-zero entries

  PetscScalar * deviceData; //< Device data

  PetscInt *diag_offsets;
  std::map<int, int> * diag_starts;
  std::vector<int> * diagonals;

  cudaStream_t stream;

  PetscScalar * deviceX;      //< Device pointer to X
  PetscScalar * deviceY;      //< Device pointer to y
  int         * deviceDiags;  //< Device pointer to diagonals specifier
} Mat_SeqSGGPU;

// Prototypes
PetscErrorCode MatDestroy_SeqSGGPU(Mat A);
PetscErrorCode MatSetGrid_SeqSGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p);
PetscErrorCode MatMult_SeqSGGPU(Mat mat, Vec x, Vec y);
PetscErrorCode MatSetValuesBlocked_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetValues_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetStencil_SeqSGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof);
PetscErrorCode MatSetUp_SeqSGGPU(Mat mat);
PetscErrorCode MatZeroEntries_SeqSGGPU(Mat A);
PetscErrorCode MatGetDiagonal_SeqSGGPU(Mat A, Vec v);
PetscErrorCode MatDiagonalScale_SeqSGGPU(Mat A, Vec ll, Vec rr);
PetscErrorCode MatGetRow_SeqSGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v);
PetscErrorCode MatRestoreRow_SeqSGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v);
PetscErrorCode MatGetRowMaxAbs_SeqSGGPU(Mat A, Vec v, PetscInt idx[]);
PetscErrorCode MatView_SeqSGGPU(Mat A, PetscViewer viewer);
PetscErrorCode MatAssemblyBegin_SeqSGGPU(Mat A, MatAssemblyType type);
PetscErrorCode MatAssemblyEnd_SeqSGGPU(Mat A, MatAssemblyType type);
PetscErrorCode MatLoad_SeqSGGPU(Mat newmat,PetscViewer viewer);
PetscErrorCode MatILUFactor_SeqSGGPU(Mat inA,IS row,IS col,const MatFactorInfo *info);
PetscErrorCode MatGetInfo_SeqSGGPU(Mat A,MatInfoType flag,MatInfo *info);
PetscErrorCode MatSolve_SeqSGGPU(Mat A,Vec bb,Vec xx);
PetscErrorCode MatSolve_SeqSGGPU_cpu(Mat A,Vec bb,Vec xx);
PetscErrorCode MatSolve_SeqSGGPU_gpu(Mat A,Vec bb,Vec xx);
PetscErrorCode MatSolve_SeqSGGPU_compare(Mat A,Vec bb,Vec xx);
PetscErrorCode MatView_SeqSGGPU_ASCII(Mat,PetscViewer);
extern PetscErrorCode MatFDColoringApply_SeqSGGPU(Mat,MatFDColoring,Vec,MatStructure*,void*);
extern PetscErrorCode MatFDColoringCreate_SeqSGGPU(Mat,ISColoring,MatFDColoring);
extern PetscErrorCode MatGetColumnIJ_SeqSGGPU(Mat,PetscInt,PetscBool,PetscBool,PetscInt *,const PetscInt *[],const PetscInt *[],PetscBool*);
EXTERN_C_BEGIN
extern PetscErrorCode MatSeqSGGPUSetPreallocation_SeqSGGPU(Mat A,PetscInt nz, const PetscInt nnz[]);
extern PetscErrorCode  MatConvert_SeqSGGPU_SeqAIJ(Mat, const MatType, MatReuse, Mat*);
extern PetscErrorCode  MatGetFactor_seqsggpu_petsc(Mat,MatFactorType,Mat*);
extern PetscErrorCode  MatConvertLU_SeqSGGPU_SeqAIJ(Mat A, Mat *AIJ_L, Mat *AIJ_U );
EXTERN_C_END
#endif
