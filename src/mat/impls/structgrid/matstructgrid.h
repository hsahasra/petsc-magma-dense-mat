
#ifndef __SEQSTRUCTGRID_H__
#define __SEQSTRUCTGRID_H__

#include "petsc-private/matimpl.h"
#include "petsc-private/vecimpl.h"



typedef enum {MAT_SINGLE,MAT_PERSIST,MAT_DEALLOC,MAT_COLLECT} MatUsageGPUFlag;
typedef enum {MAT_UNALLOC,MAT_ALLOC,MAT_GPU,MAT_CPU,MAT_SYNCHED} MatGPUFlag;



/*
The following structure defines the structgrid datatype. It is a subclass 
of generic matrix datatype and makes an implementation inheritance. 
*/
typedef struct{                             
    PetscScalar* a; 	        //data
    PetscInt* idx;		//x indices
    PetscInt* idy;		//y indices
    PetscInt* idz;		//z indices
    PetscInt stpoints;	        //num of stencil points
    PetscInt dis;		//stencil Width
    PetscInt dof;		//degrees of freedom
    PetscInt m;		        //x size
    PetscInt n;		        //y size
    PetscInt p;		        //z size
    PetscInt nz;	        //number of grid elements
    PetscInt nos;
    PetscInt lda1;
    PetscInt lda2;
    PetscInt lda3;
    PetscInt matsize;
    PetscScalar * xt;
    PetscScalar* devptr;
    MatUsageGPUFlag lifetime;
    MatGPUFlag      syncState;
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
extern PetscErrorCode MatGetRow_SeqSG(Mat , PetscInt, PetscInt *, PetscInt *[] , PetscScalar *[]);
extern PetscErrorCode MatRestoreRow_SeqSG(Mat , PetscInt, PetscInt *, PetscInt *[] , PetscScalar *[]);
extern PetscErrorCode MatGetRowMaxAbs_SeqSG(Mat, Vec, PetscInt[]);
extern PetscErrorCode MatView_SeqSG(Mat,PetscViewer);
extern PetscErrorCode MatConvert_SGtoAIJ(Mat, Mat*);


#endif


