#ifndef __SEQBSGCOMMON
#define __SEQBSGCOMMON


extern PetscInt BSG_MatMult_1(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );
extern PetscInt BSG_MatMult_2(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );
extern PetscInt BSG_MatMult_3(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );
extern PetscInt BSG_MatMult_4(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );
extern PetscInt BSG_MatMult_5(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );
extern PetscInt BSG_MatMult_6(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );
extern PetscInt BSG_MatMult_Neven(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );
extern PetscInt BSG_MatMult_Nodd(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt );

extern PetscErrorCode SetValues_SeqBSG(Mat_SeqBSG *, PetscInt , const PetscInt[], const PetscInt[] ,const PetscInt[], const PetscInt[], const PetscScalar[], InsertMode);


#endif
