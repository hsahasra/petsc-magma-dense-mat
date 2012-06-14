
#ifndef __SEQSTRUCTGRID_H__
#define __SEQSTRUCTGRID_H__

#include "petsc-private/matimpl.h"


/*
The following structure defines the structgrid datatype. It is a subclass
of generic matrix datatype and makes an implementation inheritance.
*/
typedef struct
{
PetscScalar * a; 	//data
PetscInt * idx;		//x indices
PetscInt * idy;		//y indices
PetscInt * idz;		//z indices
PetscInt stpoints;	//numberof stencil points
PetscInt dis;		//Stencil Width
PetscInt dof;		//Degrees of Freedom
PetscInt m;		//x- size
PetscInt n;		//y - size
PetscInt p;		//z size
PetscInt nz;		//number of grid elements

PetscScalar * xt;

PetscScalar * gpuMat;
}Mat_SeqSG;

extern PetscErrorCode MatCreate_SeqSG(Mat);
extern PetscErrorCode MatDestroy_SeqSG(Mat);
extern PetscErrorCode MatMult_SeqSG(Mat,Vec,Vec);
extern PetscErrorCode MatSetValues_SeqSG(Mat, PetscInt,const PetscInt[] , PetscInt,const PetscInt[],const PetscScalar[], InsertMode);
extern PetscErrorCode MatSetValuesBlocked_SeqSG(Mat, PetscInt,const PetscInt[], PetscInt,const PetscInt[],const PetscScalar[], InsertMode);
extern PetscErrorCode MatSetStencil_SeqSG(Mat, PetscInt,const PetscInt[],const PetscInt[], PetscInt );
extern PetscErrorCode MatSetUpPreallocation_SeqSG(Mat);
extern PetscErrorCode MatSetGrid_SeqSG(Mat,PetscInt, PetscInt, PetscInt);
extern PetscErrorCode MatZeroEntries_SeqSG(Mat);
extern PetscErrorCode MatGetDiagonal_SeqSG(Mat, Vec);
extern PetscErrorCode MatDiagonalScale_SeqSG(Mat, Vec, Vec);
extern PetscErrorCode MatGetRow_SeqSG(Mat , PetscInt, PetscInt *, PetscInt *[] , PetscScalar *[]);
extern PetscErrorCode MatRestoreRow_SeqSG(Mat , PetscInt, PetscInt *, PetscInt *[] , PetscScalar *[]);
extern PetscErrorCode MatGetRowMaxAbs_SeqSG(Mat, Vec, PetscInt[]);
extern PetscErrorCode MatView_SeqSG(Mat,PetscViewer);
extern PetscErrorCode MatConvert_SGtoAIJ(Mat, Mat*);
#endif

