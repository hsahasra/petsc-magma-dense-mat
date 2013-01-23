#ifndef __SEQBSGCOMMON
#define __SEQBSGCOMMON

#define WORKINGSETSIZE 10000

PetscInt BSG_MatMult_1(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_2(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_3(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_4(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_5(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_6(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_7(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_Neven(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );
PetscInt BSG_MatMult_Nodd(PetscScalar **, const PetscScalar *, PetscScalar *,PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *, PetscInt , const PetscInt * , const PetscInt * , const PetscInt * );

PetscInt BSG_MatMult_8(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_10(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_12(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_14(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_16(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_18(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_20(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_22(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_24(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_26(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_28(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
PetscInt BSG_MatMult_30(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );

PetscErrorCode GetValues_Matrix_SeqBSG(Mat_SeqBSG *, PetscInt , const PetscInt[], const PetscInt[], PetscScalar*);

PetscErrorCode SetValues_Matrix_SeqBSG(Mat_SeqBSG *, PetscInt , const PetscInt[], const PetscInt[], const PetscInt[], const  PetscScalar[], InsertMode );

PetscErrorCode SetValuesBlocked_Matrix_SeqBSG(Mat_SeqBSG *, PetscInt , const PetscInt[], const PetscInt[], const  PetscScalar[], InsertMode );

PetscErrorCode MatGetSetValues_SeqBSG(Mat , PetscInt ,const PetscInt[], PetscInt,const PetscInt[],const PetscScalar *, InsertMode , PetscBool, PetscBool );

PetscErrorCode MatSetUpPreallocation_SubMatrix_SeqBSG(Mat);

PetscErrorCode MatSetUpSubRegion_SeqBSG(Mat);

PetscErrorCode MatSetUpRegion_SeqBSG(Mat);

PetscErrorCode MatSetUpPreallocation_Matrix_SeqBSG(Mat);

PetscErrorCode MatGetSubMatrix_SeqBSG_Private(Mat ,IS ,IS ,MatReuse,Mat *);
#endif
