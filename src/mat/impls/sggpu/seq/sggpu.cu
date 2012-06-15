/// SGGPU Matrix Type

#define PETSCMAT_DLL

#include "petsc-private/matimpl.h"
#include "sggpu.h"

#include <stdio.h>
#include <cuda.h>


// Debugging flags
#define _TRACE 1

// Prototypes
PetscErrorCode MatDestroy_SeqSGGPU(Mat A);
PetscErrorCode MatSetGrid_SeqSGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p);
PetscErrorCode MatMult_SeqSGGPU(Mat mat, Vec x, Vec y);
PetscErrorCode MatSetValuesBlocked_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetValues_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetStencil_SeqSGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof);
PetscErrorCode MatSetUpPreallocation_SeqSGGPU(Mat mat);
PetscErrorCode MatZeroEntries_SeqSGGPU(Mat A);
PetscErrorCode MatGetDiagonal_SeqSGGPU(Mat A, Vec v);
PetscErrorCode MatDiagonalScale_SeqSGGPU(Mat A, Vec ll, Vec rr);
PetscErrorCode MatGetRow_SeqSGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v);
PetscErrorCode MatRestoreRow_SeqSGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v);
PetscErrorCode MatGetRowMaxAbs_SeqSGGPU(Mat A, Vec v, PetscInt idx[]);
PetscErrorCode MatView_SeqSGGPU(Mat A, PetscViewer viewer);


// Matrix function table
static struct _MatOps MatOps_Values = {
/*0*/ MatSetValues_SeqSGGPU,MatGetRow_SeqSGGPU,MatRestoreRow_SeqSGGPU,MatMult_SeqSGGPU,0,
/*5*/0,0,0,0,0,
/*10*/0,0,0,0,0,
/*15*/0,0,MatGetDiagonal_SeqSGGPU,MatDiagonalScale_SeqSGGPU,0,
/*20*/0,0,0,MatZeroEntries_SeqSGGPU,0,
/*25*/0,0,0,0,MatSetUpPreallocation_SeqSGGPU,
/*30*/0,0,0,0,0,
/*35*/0,0,0,0,0,
/*40*/0,0,0,0,0,
/*45*/0,0,0,0,0,
/*50*/0,0,0,0,0,
/*55*/0,0,0,MatSetValuesBlocked_SeqSGGPU,0,
/*60*/MatDestroy_SeqSGGPU,MatView_SeqSGGPU,0,0,0,
/*65*/0,0,MatSetValues_SeqSGGPU,0,MatGetRowMaxAbs_SeqSGGPU,
/*70*/0,0,0,0,0,
/*75*/0,0,0,0,0,
/*80*/0,0,0,0,0,
/*85*/0,0,MatSetValuesBlocked_SeqSGGPU,0,0,
/*90*/0,0,0,0,0,
/*95*/0,0,0,0,0,
/*100*/0,0,0,0,0,
/*105*/0,0,0,0,0,
/*110*/0,0,0,0,0,
/*115*/MatCreate_SeqSGGPU,0,0,0,0,
/*120*/0,0,0,0,0,
/*125*/0,0,0,0,0,
/*130*/0,0,0,0,0,
/*135*/0,0,0,0,MatSetStencil_SeqSGGPU,
/*140*/MatSetGrid_SeqSGGPU
};


EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatCreate_SeqSGGPU"
PetscErrorCode MatCreate_SeqSGGPU(Mat A)
{
  Mat_SeqSGGPU * mat;
  PetscErrorCode ierr;
  PetscMPIInt size;

  PetscFunctionBegin;

#if _TRACE
  printf("[SeqSGGPU] MatCreate_SeqSGGPU\n");
#endif

  ierr = MPI_Comm_size(((PetscObject)A)->comm, &size); CHKERRQ(ierr);
  if (size > 1)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Comm must be size 1");

  // Create internal matrix structure
  ierr = PetscMalloc(sizeof(Mat_SeqSGGPU), &mat); CHKERRQ(ierr);
  memset(mat, 0, sizeof(Mat_SeqSGGPU));

  // Fill out PETSc matrix structure
  A->data = mat;
  memcpy(A->ops, &MatOps_Values, sizeof(struct _MatOps));
  A->same_nonzero= PETSC_FALSE;
  A->spptr = 0;

  PetscFunctionReturn(0);
}
EXTERN_C_END


#undef __FUNCT__
#define __FUNCT__ "MatDestroy_SeqSGGPU"
PetscErrorCode MatDestroy_SeqSGGPU(Mat A)
{
  Mat_SeqSGGPU * mat;
  PetscErrorCode ierr;

  PetscFunctionBegin;

#if _TRACE
  printf("[SeqSGGPU] MatDestroy_SeqSGGPU\n");
#endif

  mat = (Mat_SeqSGGPU*)A->data;

  if (mat->hostData) {
    ierr = PetscFree(mat->hostData); CHKERRQ(ierr);
  }
  if (mat->deviceData) {
    cudaFree(mat->deviceData);
  }
  PetscFree(mat); CHKERRQ(ierr);

  ierr = PetscObjectChangeTypeName((PetscObject)A, 0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetGrid_SeqSGGPU"
PetscErrorCode MatSetGrid_SeqSGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetGrid_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqSGGPU"
PetscErrorCode MatMult_SeqSGGPU(Mat mat, Vec x, Vec y)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatMult_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_SeqSGGPU"
PetscErrorCode MatSetValuesBlocked_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetValuesBlocked_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValues_SeqSGGPU"
PetscErrorCode MatSetValues_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetValues_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetStencil_SeqSGGPU"
PetscErrorCode MatSetStencil_SeqSGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetStencil_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetUpPreallocation_SeqSGGPU"
PetscErrorCode MatSetUpPreallocation_SeqSGGPU(Mat mat)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetUpPreallocation_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_SeqSGGPU"
PetscErrorCode MatZeroEntries_SeqSGGPU(Mat A)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatZeroEntries_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_SeqSGGPU"
PetscErrorCode MatGetDiagonal_SeqSGGPU(Mat A, Vec v)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatGetDiagonal_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_SeqSGGPU"
PetscErrorCode MatDiagonalScale_SeqSGGPU(Mat A, Vec ll, Vec rr)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatDiagonalScale_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRow_SeqSGGPU"
PetscErrorCode MatGetRow_SeqSGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatGetRow_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_SeqSGGPU"
PetscErrorCode MatRestoreRow_SeqSGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatRestoreRow_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_SeqSGGPU"
PetscErrorCode MatGetRowMaxAbs_SeqSGGPU(Mat A, Vec v, PetscInt idx[])
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatGetRowMaxAbs_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatView_SeqSGGPU"
PetscErrorCode MatView_SeqSGGPU(Mat A, PetscViewer viewer)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatView_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}
