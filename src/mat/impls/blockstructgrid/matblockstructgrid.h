
#ifndef __SEQBSTRUCTGRID_H__
#define __SEQBSTRUCTGRID_H__

#include "private/matimpl.h"


/*
The following structure defines the block structgrid datatype. It is a subclass 
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
PetscInt nz;		//length of grid elements

// Specific to block 
PetscInt bs; // block size , usually dof*dof
PetscInt tnz; // number of non - zero grid elements
}Mat_SeqBSG;

extern PetscErrorCode MatCreate_SeqBSG(Mat);
extern PetscErrorCode MatDestroy_SeqBSG(Mat);
extern PetscErrorCode MatMult_SeqBSG(Mat,Vec,Vec);
extern PetscErrorCode MatSetValues_SeqBSG(Mat, PetscInt,const PetscInt[] , PetscInt,const PetscInt[],const PetscScalar[], InsertMode); 
extern PetscErrorCode MatSetValuesBlocked_SeqBSG(Mat, PetscInt,const PetscInt[], PetscInt,const PetscInt[],const PetscScalar[], InsertMode);
extern PetscErrorCode MatSetStencil_SeqBSG(Mat, PetscInt,const PetscInt[],const PetscInt[], PetscInt );
extern PetscErrorCode MatSetUpPreallocation_SeqBSG(Mat);
extern PetscErrorCode MatSetGrid_SeqBSG(Mat,PetscInt, PetscInt, PetscInt);
extern PetscErrorCode MatZeroEntries_SeqBSG(Mat);
extern PetscErrorCode MatView_SeqBSG(Mat, PetscViewer);
#endif

