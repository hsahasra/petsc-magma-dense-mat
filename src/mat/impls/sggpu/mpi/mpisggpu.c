#include <../src/mat/impls/sggpu/mpi/mpisggpu.h>





#undef __FUNCT__
#define __FUNCT__ "MatDestroy_MPISGGPU"
PetscErrorCode MatDestroy_MPISGGPU(Mat A) {
  MatDestroy_MPISGGPU(A);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetGrid_MPISGGPU"
PetscErrorCode MatSetGrid_MPISGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p) {
  MatSetGrid_SeqSGGPU(B, m, n, p);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatMult_MPISGGPU"
PetscErrorCode MatMult_MPISGGPU(Mat mat, Vec x, Vec y) {
  MatMult_MPISGGPU(mat, x, y); 
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_MPISGGPU"
PetscErrorCode MatSetValuesBlocked_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is) {
  MatSetValuesBlocked_SeqSGGPU(A, nrow, irow, ncol, icol, y, is);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValues_MPISGGPU"
PetscErrorCode MatSetValues_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is) {
  MatSetValues_MPISGGPU(A, nrow, irow, ncol, icol, y, is);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetStencil_MPISGGPU"
PetscErrorCode MatSetStencil_MPISGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof) {
  MatSetStencil_MPISGGPU(A, dim, dims, starts, dof);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetUp_MPISGGPU"
PetscErrorCode MatSetUp_MPISGGPU(Mat mat) {
  MatSetUp_SeqSGGPU(mat);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_MPISGGPU"
PetscErrorCode MatZeroEntries_MPISGGPU(Mat A) {
  MatZeroEntries_SeqSGGPU(A);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_MPISGGPU"
PetscErrorCode MatGetDiagonal_MPISGGPU(Mat A, Vec v) {
  MatGetDiagonal_SeqSGGPU(A, v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_MPISGGPU"
PetscErrorCode MatDiagonalScale_MPISGGPU(Mat A, Vec ll, Vec rr) {
  MatDiagonalScale_SeqSGGPU(A, ll, rr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRow_MPISGGPU"
PetscErrorCode MatGetRow_MPISGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v) {
  MatGetRow_SeqSGGPU(A, row, nz, idx , v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_MPISGGPU"
PetscErrorCode MatRestoreRow_MPISGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v) {
  MatRestoreRow_SeqSGGPU(A, row, nz, idx, v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_MPISGGPU"
PetscErrorCode MatGetRowMaxAbs_MPISGGPU(Mat A, Vec v, PetscInt idx[]) {
  MatGetRowMaxAbs_SeqSGGPU(A, v, idx);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatView_MPISGGPU"
PetscErrorCode MatView_MPISGGPU(Mat A, PetscViewer viewer) {
  MatView_SeqSGGPU(A, viewer);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyBegin_MPISGGPU"
PetscErrorCode MatAssemblyBegin_MPISGGPU(Mat A, MatAssemblyType type) {
  MatAssemblyBegin_SeqSGGPU(A, type);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyEnd_MPISGGPU"
PetscErrorCode MatAssemblyEnd_MPISGGPU(Mat A, MatAssemblyType type) {
  MatAssemblyEnd_SeqSGGPU(A, type);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatView_MPISGGPU_ASCII"
PetscErrorCode MatView_MPISGGPU_ASCII(Mat A, PetscViewer viewer) {
  MatView_SeqSGGPU_ASCII(A, viewer);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatFDColoringApply_MPISGGPU"
PetscErrorCode  MatFDColoringApply_MPISGGPU(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx) {
  MatFDColoringApply_SeqSGGPU(J, coloring, x1, flag, sctx);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatFDColoringCreate_MPISGGPU"
PetscErrorCode MatFDColoringCreate_MPISGGPU(Mat mat,ISColoring iscoloring,MatFDColoring c) {
  MatFDColoringCreate_SeqSGGPU(mat, iscoloring, c);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetColumnIJ_MPISGGPU"
PetscErrorCode MatGetColumnIJ_MPISGGPU(Mat A,PetscInt oshift,PetscBool  symmetric,PetscBool  inodecompressed,PetscInt *nn, const PetscInt *ia[], const PetscInt *ja[],PetscBool  *done) {
  MatGetColumnIJ_SeqSGGPU(A, oshift, symmetric, inodecompressed, nn, ia, ja, done);
  PetscFunctionReturn(0);
}
