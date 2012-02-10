#ifndef __SEQBSGCOMMON
#define __SEQBSGCOMMON

#define WORKINGSETSIZE 100

extern PetscInt BSG_MatMult_1(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *);
extern PetscInt BSG_MatMult_2(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *);
extern PetscInt BSG_MatMult_3(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt , const PetscInt *);
extern PetscInt BSG_MatMult_4(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_5(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_6(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_8(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_10(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_12(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_14(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_16(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_18(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_20(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_22(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_24(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_26(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_28(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_30(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );
extern PetscInt BSG_MatMult_Neven(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt ,const PetscInt *);
extern PetscInt BSG_MatMult_Nodd(PetscScalar **, const PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscInt *, PetscInt, PetscInt, PetscInt,PetscInt, PetscInt, PetscInt , PetscInt, const PetscInt * );

extern PetscErrorCode SetValues_SeqBSG(Mat_SeqBSG *, PetscInt , const PetscInt[], const PetscInt[] ,const PetscInt[], const PetscInt[], const PetscScalar[], InsertMode);


#endif
