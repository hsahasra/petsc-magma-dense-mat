#include <string.h>
#include <immintrin.h>
#include "../src/mat/impls/blockstructgrid/matblockstructgrid.h"
#include "../src/mat/impls/blockstructgrid/commonfunctions.h"

#define min(a,b) (a) < (b) ? (a) : (b)

#define setupct(pos,iter) \
                          count = stpoffset[pos] + iter*nos;\
                          ct0 = ctl[count + 0];\
                          ct1 = ctl[count + 1];\
                          ct2 = ctl[count + 2];\
                          ct3 = ctl[count + 3];\
                          ct4 = ctl[count + 4];\
                          ct5 = ctl[count + 5];\
                          ct6 = ctl[count + 6]

#define inline_8(l) \
                    mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                    mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                    mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                    mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                    mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                    mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                    mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                    mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                    mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                    mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                    mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                    mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                    mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                    mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                    mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                    mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                    mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                    mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                    mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                    mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                    mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                    mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                    mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                    mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                    mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                    mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                    mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                    mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                    mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                    mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                    mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                    mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                    mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                    mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                    mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                    mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc4));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc5));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc6));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc7));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc8));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc9));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc10));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc11));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc12));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc13));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc14));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc15));\
                    msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc16));\
                    msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc17));\
                    msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc18));\
                    msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc19));\
                    msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc20));\
                    msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc21));\
                    msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc22));\
                    msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc23));\
                    msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc24));\
                    msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc25));\
                    msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc26));\
                    msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc27));\
                    msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc28));\
                    msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc29));\
                    msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc30));\
                    msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc31))

#define setup_8(offset) \
                        t1 = k*dof; t2 = (k-(offset))*bs;\
                        msum0 = _mm_set_pd(0,0);\
                        msum1 = _mm_set_pd(0,0);\
                        msum2 = _mm_set_pd(0,0);\
                        msum3 = _mm_set_pd(0,0);\
                        msum4 = _mm_set_pd(0,0);\
                        msum5 = _mm_set_pd(0,0);\
                        msum6 = _mm_set_pd(0,0);\
                        msum7 = _mm_set_pd(0,0)

#define save_8() \
                 msum0 = _mm_hadd_pd(msum0, msum1);\
                 msum2 = _mm_hadd_pd(msum2, msum3);\
                 msum4 = _mm_hadd_pd(msum4, msum5);\
                 msum6 = _mm_hadd_pd(msum6, msum7);\
                 _mm_storeu_pd(y + t1 + 0,msum0);\
                 _mm_storeu_pd(y + t1 + 2,msum2);\
                 _mm_storeu_pd(y + t1 + 4,msum4);\
                 _mm_storeu_pd(y + t1 + 6,msum6)

PetscErrorCode BSG_MatMult_8(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_8(0);
            inline_8(3);
            inline_8(4);
            inline_8(5);
            inline_8(6);
            save_8();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_8(1);
            inline_8(2);
            inline_8(3);
            inline_8(4);
            inline_8(5);
            inline_8(6);
            save_8();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_8(lda3);
            inline_8(1);
            inline_8(2);
            inline_8(3);
            inline_8(4);
            inline_8(5);
            inline_8(6);
            save_8();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_8(lda2);
            inline_8(0);
            inline_8(1);
            inline_8(2);
            inline_8(3);
            inline_8(4);
            inline_8(5);
            inline_8(6);
            save_8();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_8(lda1 - lda2);
            inline_8(0);
            inline_8(1);
            inline_8(2);
            inline_8(3);
            inline_8(4);
            inline_8(5);
            save_8();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_8(lda1 - lda3);
            inline_8(0);
            inline_8(1);
            inline_8(2);
            inline_8(3);
            inline_8(4);
            save_8();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_8(lda1 - 1);
            inline_8(0);
            inline_8(1);
            inline_8(2);
            inline_8(3);
            save_8();
        }
    }

PetscFunctionReturn(0);
}

#define inline_10(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc5));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc6));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc7));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc8));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc9));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc10));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc11));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc12));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc13));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc14));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc15));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc16));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc17));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc18));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc19));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc20));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc21));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc22));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc23));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc24));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc25));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc26));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc27));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc28));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc29));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc30));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc31));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc32));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc33));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc34));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc35));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc36));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc37));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc38));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc39));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc40));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc41));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc42));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc43));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc44));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc45));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc46));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc47));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc48));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc49))

#define setup_10(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0)

#define save_10() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8)

PetscErrorCode BSG_MatMult_10(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_10(0);
            inline_10(3);
            inline_10(4);
            inline_10(5);
            inline_10(6);
            save_10();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_10(1);
            inline_10(2);
            inline_10(3);
            inline_10(4);
            inline_10(5);
            inline_10(6);
            save_10();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_10(lda3);
            inline_10(1);
            inline_10(2);
            inline_10(3);
            inline_10(4);
            inline_10(5);
            inline_10(6);
            save_10();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_10(lda2);
            inline_10(0);
            inline_10(1);
            inline_10(2);
            inline_10(3);
            inline_10(4);
            inline_10(5);
            inline_10(6);
            save_10();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_10(lda1 - lda2);
            inline_10(0);
            inline_10(1);
            inline_10(2);
            inline_10(3);
            inline_10(4);
            inline_10(5);
            save_10();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_10(lda1 - lda3);
            inline_10(0);
            inline_10(1);
            inline_10(2);
            inline_10(3);
            inline_10(4);
            save_10();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_10(lda1 - 1);
            inline_10(0);
            inline_10(1);
            inline_10(2);
            inline_10(3);
            save_10();
        }
    }

PetscFunctionReturn(0);
}

#define inline_12(l,offset) \
                            t1 = k * dof + 0;\
                            mx0 = _mm_loadu_pd(xt##l + t1);\
                            mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                            t2 = (k-(offset))*144+0 + 0;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                            msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*144+48 + 0;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc0, mc2));\
                            msum3 = _mm_add_pd(msum3, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*144+96 + 0;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum4 = _mm_add_pd(msum4, _mm_hadd_pd(mc0, mc2));\
                            msum5 = _mm_add_pd(msum5, _mm_hadd_pd(mc4, mc6));\
                            t1 = k * dof + 4;\
                            mx0 = _mm_loadu_pd(xt##l + t1);\
                            mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                            t2 = (k-(offset))*144+0 + 4;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                            msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*144+48 + 4;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc0, mc2));\
                            msum3 = _mm_add_pd(msum3, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*144+96 + 4;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum4 = _mm_add_pd(msum4, _mm_hadd_pd(mc0, mc2));\
                            msum5 = _mm_add_pd(msum5, _mm_hadd_pd(mc4, mc6));\
                            t1 = k * dof + 8;\
                            mx0 = _mm_loadu_pd(xt##l + t1);\
                            mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                            t2 = (k-(offset))*144+0 + 8;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                            msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*144+48 + 8;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc0, mc2));\
                            msum3 = _mm_add_pd(msum3, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*144+96 + 8;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 12);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 12 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 24);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 24 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 36);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 36 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum4 = _mm_add_pd(msum4, _mm_hadd_pd(mc0, mc2));\
                            msum5 = _mm_add_pd(msum5, _mm_hadd_pd(mc4, mc6))

#define setup_12() \
                   msum0 = _mm_set_pd(0,0);\
                   msum1 = _mm_set_pd(0,0);\
                   msum2 = _mm_set_pd(0,0);\
                   msum3 = _mm_set_pd(0,0);\
                   msum4 = _mm_set_pd(0,0);\
                   msum5 = _mm_set_pd(0,0)

#define save_12() \
                  _mm_storeu_pd(y + k * dof + 0, msum0);\
                  _mm_storeu_pd(y + k * dof + 2, msum1);\
                  _mm_storeu_pd(y + k * dof + 4, msum2);\
                  _mm_storeu_pd(y + k * dof + 6, msum3);\
                  _mm_storeu_pd(y + k * dof + 8, msum4);\
                  _mm_storeu_pd(y + k * dof + 10, msum5)

PetscErrorCode BSG_MatMult_12(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m128d mx0, mx1, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, msum0, msum1, msum2, msum3, msum4, msum5;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_12();
            inline_12(3,0);
            inline_12(4,0);
            inline_12(5,0);
            inline_12(6,0);
            save_12();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_12();
            inline_12(2,1);
            inline_12(3,1);
            inline_12(4,1);
            inline_12(5,1);
            inline_12(6,1);
            save_12();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_12();
            inline_12(1,lda3);
            inline_12(2,lda3);
            inline_12(3,lda3);
            inline_12(4,lda3);
            inline_12(5,lda3);
            inline_12(6,lda3);
            save_12();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_12();
            inline_12(0,lda2);
            inline_12(1,lda2);
            inline_12(2,lda2);
            inline_12(3,lda2);
            inline_12(4,lda2);
            inline_12(5,lda2);
            inline_12(6,lda2);
            save_12();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_12();
            inline_12(0,lda1 - lda2);
            inline_12(1,lda1 - lda2);
            inline_12(2,lda1 - lda2);
            inline_12(3,lda1 - lda2);
            inline_12(4,lda1 - lda2);
            inline_12(5,lda1 - lda2);
            save_12();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_12();
            inline_12(0,lda1 - lda3);
            inline_12(1,lda1 - lda3);
            inline_12(2,lda1 - lda3);
            inline_12(3,lda1 - lda3);
            inline_12(4,lda1 - lda3);
            save_12();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_12();
            inline_12(0,lda1 - 1);
            inline_12(1,lda1 - 1);
            inline_12(2,lda1 - 1);
            inline_12(3,lda1 - 1);
            save_12();
        }
    }

PetscFunctionReturn(0);
}

#define inline_14(l,offset) \
                            t1 = k * dof + 0;\
                            mx0 = _mm_loadu_pd(xt##l + t1);\
                            mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                            t2 = (k-(offset))*196+0 + 0;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                            msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+56 + 0;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc0, mc2));\
                            msum3 = _mm_add_pd(msum3, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+112 + 0;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum4 = _mm_add_pd(msum4, _mm_hadd_pd(mc0, mc2));\
                            msum5 = _mm_add_pd(msum5, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+168 + 0;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            msum6 = _mm_add_pd(msum6, _mm_hadd_pd(mc0, mc2));\
                            t1 = k * dof + 4;\
                            mx0 = _mm_loadu_pd(xt##l + t1);\
                            mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                            t2 = (k-(offset))*196+0 + 4;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                            msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+56 + 4;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc0, mc2));\
                            msum3 = _mm_add_pd(msum3, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+112 + 4;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum4 = _mm_add_pd(msum4, _mm_hadd_pd(mc0, mc2));\
                            msum5 = _mm_add_pd(msum5, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+168 + 4;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            msum6 = _mm_add_pd(msum6, _mm_hadd_pd(mc0, mc2));\
                            t1 = k * dof + 8;\
                            mx0 = _mm_loadu_pd(xt##l + t1);\
                            mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                            t2 = (k-(offset))*196+0 + 8;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                            msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+56 + 8;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc0, mc2));\
                            msum3 = _mm_add_pd(msum3, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+112 + 8;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc5 = _mm_loadu_pd(ct##l + t2 + 28 + 2);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc7 = _mm_loadu_pd(ct##l + t2 + 42 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            mc4 = _mm_add_pd(_mm_mul_pd(mx0, mc4), _mm_mul_pd(mx1, mc5));\
                            mc6 = _mm_add_pd(_mm_mul_pd(mx0, mc6), _mm_mul_pd(mx1, mc7));\
                            msum4 = _mm_add_pd(msum4, _mm_hadd_pd(mc0, mc2));\
                            msum5 = _mm_add_pd(msum5, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+168 + 8;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc3 = _mm_loadu_pd(ct##l + t2 + 14 + 2);\
                            mc0 = _mm_add_pd(_mm_mul_pd(mx0, mc0), _mm_mul_pd(mx1, mc1));\
                            mc2 = _mm_add_pd(_mm_mul_pd(mx0, mc2), _mm_mul_pd(mx1, mc3));\
                            msum6 = _mm_add_pd(msum6, _mm_hadd_pd(mc0, mc2));\
                            t1 = k * dof + 12;\
                            mx0 = _mm_loadu_pd(xt##l + t1);\
                            t2 = (k-(offset))*196+0 + 12;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc0 = _mm_mul_pd(mx0, mc0);\
                            mc2 = _mm_mul_pd(mx0, mc2);\
                            mc4 = _mm_mul_pd(mx0, mc4);\
                            mc6 = _mm_mul_pd(mx0, mc6);\
                            msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                            msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+56 + 12;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc0 = _mm_mul_pd(mx0, mc0);\
                            mc2 = _mm_mul_pd(mx0, mc2);\
                            mc4 = _mm_mul_pd(mx0, mc4);\
                            mc6 = _mm_mul_pd(mx0, mc6);\
                            msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc0, mc2));\
                            msum3 = _mm_add_pd(msum3, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+112 + 12;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc4 = _mm_loadu_pd(ct##l + t2 + 28);\
                            mc6 = _mm_loadu_pd(ct##l + t2 + 42);\
                            mc0 = _mm_mul_pd(mx0, mc0);\
                            mc2 = _mm_mul_pd(mx0, mc2);\
                            mc4 = _mm_mul_pd(mx0, mc4);\
                            mc6 = _mm_mul_pd(mx0, mc6);\
                            msum4 = _mm_add_pd(msum4, _mm_hadd_pd(mc0, mc2));\
                            msum5 = _mm_add_pd(msum5, _mm_hadd_pd(mc4, mc6));\
                            t2 = (k-(offset))*196+168 + 12;\
                            mc0 = _mm_loadu_pd(ct##l + t2);\
                            mc2 = _mm_loadu_pd(ct##l + t2 + 14);\
                            mc0 = _mm_mul_pd(mx0, mc0);\
                            mc2 = _mm_mul_pd(mx0, mc2);\
                            msum6 = _mm_add_pd(msum6, _mm_hadd_pd(mc0, mc2))

#define setup_14() \
                   msum0 = _mm_set_pd(0,0);\
                   msum1 = _mm_set_pd(0,0);\
                   msum2 = _mm_set_pd(0,0);\
                   msum3 = _mm_set_pd(0,0);\
                   msum4 = _mm_set_pd(0,0);\
                   msum5 = _mm_set_pd(0,0);\
                   msum6 = _mm_set_pd(0,0)

#define save_14() \
                  _mm_storeu_pd(y + k * dof + 0, msum0);\
                  _mm_storeu_pd(y + k * dof + 2, msum1);\
                  _mm_storeu_pd(y + k * dof + 4, msum2);\
                  _mm_storeu_pd(y + k * dof + 6, msum3);\
                  _mm_storeu_pd(y + k * dof + 8, msum4);\
                  _mm_storeu_pd(y + k * dof + 10, msum5);\
                  _mm_storeu_pd(y + k * dof + 12, msum6)

PetscErrorCode BSG_MatMult_14(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m128d mx0, mx1, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, msum0, msum1, msum2, msum3, msum4, msum5, msum6;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_14();
            inline_14(3,0);
            inline_14(4,0);
            inline_14(5,0);
            inline_14(6,0);
            save_14();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_14();
            inline_14(2,1);
            inline_14(3,1);
            inline_14(4,1);
            inline_14(5,1);
            inline_14(6,1);
            save_14();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_14();
            inline_14(1,lda3);
            inline_14(2,lda3);
            inline_14(3,lda3);
            inline_14(4,lda3);
            inline_14(5,lda3);
            inline_14(6,lda3);
            save_14();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_14();
            inline_14(0,lda2);
            inline_14(1,lda2);
            inline_14(2,lda2);
            inline_14(3,lda2);
            inline_14(4,lda2);
            inline_14(5,lda2);
            inline_14(6,lda2);
            save_14();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_14();
            inline_14(0,lda1 - lda2);
            inline_14(1,lda1 - lda2);
            inline_14(2,lda1 - lda2);
            inline_14(3,lda1 - lda2);
            inline_14(4,lda1 - lda2);
            inline_14(5,lda1 - lda2);
            save_14();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_14();
            inline_14(0,lda1 - lda3);
            inline_14(1,lda1 - lda3);
            inline_14(2,lda1 - lda3);
            inline_14(3,lda1 - lda3);
            inline_14(4,lda1 - lda3);
            save_14();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_14();
            inline_14(0,lda1 - 1);
            inline_14(1,lda1 - 1);
            inline_14(2,lda1 - 1);
            inline_14(3,lda1 - 1);
            save_14();
        }
    }

PetscFunctionReturn(0);
}


#define inline_16(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc8));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc9));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc10));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc11));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc12));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc13));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc15));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc16));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc17));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc18));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc19));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc20));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc21));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc22));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc23));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc24));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc25));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc26));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc27));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc28));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc29));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc30));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc31));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc32));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc33));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc34));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc35));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc36));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc37));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc38));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc39));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc40));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc41));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc42));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc43));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc44));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc45));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc46));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc47));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc48));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc49));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc50));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc51));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc52));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc53));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc54));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc55));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc56));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc57));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc58));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc59));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc60));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc61));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc62));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc63));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc64));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc65));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc66));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc67));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc68));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc69));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc70));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc71));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc72));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc73));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc74));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc75));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc76));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc77));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc78));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc79));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc80));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc81));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc82));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc83));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc84));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc85));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc86));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc87));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc88));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc89));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc90));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc91));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc92));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc93));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc94));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc95));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc96));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc97));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc98));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc99));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc100));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc101));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc102));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc103));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc104));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc105));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc106));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc107));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc108));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc109));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc110));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc111));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc112));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc113));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc114));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc115));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc116));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc117));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc118));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc119));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc120));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc121));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc122));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc123));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc124));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc125));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc126));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc127))

#define setup_16(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0)

#define save_16() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14)

PetscErrorCode BSG_MatMult_16(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_16(0);
            inline_16(3);
            inline_16(4);
            inline_16(5);
            inline_16(6);
            save_16();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_16(1);
            inline_16(2);
            inline_16(3);
            inline_16(4);
            inline_16(5);
            inline_16(6);
            save_16();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_16(lda3);
            inline_16(1);
            inline_16(2);
            inline_16(3);
            inline_16(4);
            inline_16(5);
            inline_16(6);
            save_16();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_16(lda2);
            inline_16(0);
            inline_16(1);
            inline_16(2);
            inline_16(3);
            inline_16(4);
            inline_16(5);
            inline_16(6);
            save_16();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_16(lda1 - lda2);
            inline_16(0);
            inline_16(1);
            inline_16(2);
            inline_16(3);
            inline_16(4);
            inline_16(5);
            save_16();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_16(lda1 - lda3);
            inline_16(0);
            inline_16(1);
            inline_16(2);
            inline_16(3);
            inline_16(4);
            save_16();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_16(lda1 - 1);
            inline_16(0);
            inline_16(1);
            inline_16(2);
            inline_16(3);
            save_16();
        }
    }

PetscFunctionReturn(0);
}

#define inline_18(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mx8 = _mm_loadu_pd(xt##l + t1 + 16);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     mc128 = _mm_loadu_pd(ct##l + t2 + 256);\
                     mc129 = _mm_loadu_pd(ct##l + t2 + 258);\
                     mc130 = _mm_loadu_pd(ct##l + t2 + 260);\
                     mc131 = _mm_loadu_pd(ct##l + t2 + 262);\
                     mc132 = _mm_loadu_pd(ct##l + t2 + 264);\
                     mc133 = _mm_loadu_pd(ct##l + t2 + 266);\
                     mc134 = _mm_loadu_pd(ct##l + t2 + 268);\
                     mc135 = _mm_loadu_pd(ct##l + t2 + 270);\
                     mc136 = _mm_loadu_pd(ct##l + t2 + 272);\
                     mc137 = _mm_loadu_pd(ct##l + t2 + 274);\
                     mc138 = _mm_loadu_pd(ct##l + t2 + 276);\
                     mc139 = _mm_loadu_pd(ct##l + t2 + 278);\
                     mc140 = _mm_loadu_pd(ct##l + t2 + 280);\
                     mc141 = _mm_loadu_pd(ct##l + t2 + 282);\
                     mc142 = _mm_loadu_pd(ct##l + t2 + 284);\
                     mc143 = _mm_loadu_pd(ct##l + t2 + 286);\
                     mc144 = _mm_loadu_pd(ct##l + t2 + 288);\
                     mc145 = _mm_loadu_pd(ct##l + t2 + 290);\
                     mc146 = _mm_loadu_pd(ct##l + t2 + 292);\
                     mc147 = _mm_loadu_pd(ct##l + t2 + 294);\
                     mc148 = _mm_loadu_pd(ct##l + t2 + 296);\
                     mc149 = _mm_loadu_pd(ct##l + t2 + 298);\
                     mc150 = _mm_loadu_pd(ct##l + t2 + 300);\
                     mc151 = _mm_loadu_pd(ct##l + t2 + 302);\
                     mc152 = _mm_loadu_pd(ct##l + t2 + 304);\
                     mc153 = _mm_loadu_pd(ct##l + t2 + 306);\
                     mc154 = _mm_loadu_pd(ct##l + t2 + 308);\
                     mc155 = _mm_loadu_pd(ct##l + t2 + 310);\
                     mc156 = _mm_loadu_pd(ct##l + t2 + 312);\
                     mc157 = _mm_loadu_pd(ct##l + t2 + 314);\
                     mc158 = _mm_loadu_pd(ct##l + t2 + 316);\
                     mc159 = _mm_loadu_pd(ct##l + t2 + 318);\
                     mc160 = _mm_loadu_pd(ct##l + t2 + 320);\
                     mc161 = _mm_loadu_pd(ct##l + t2 + 322);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx8, mc8));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc9));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc10));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc11));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc12));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc13));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc15));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc16));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx8, mc17));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc18));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc19));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc20));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc21));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc22));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc23));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc24));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc25));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx8, mc26));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc27));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc28));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc29));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc30));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc31));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc32));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc33));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc34));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx8, mc35));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc36));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc37));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc38));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc39));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc40));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc41));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc42));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc43));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx8, mc44));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc45));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc46));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc47));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc48));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc49));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc50));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc51));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc52));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx8, mc53));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc54));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc55));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc56));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc57));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc58));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc59));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc60));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc61));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx8, mc62));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc63));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc64));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc65));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc66));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc67));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc68));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc69));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc70));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx8, mc71));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc72));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc73));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc74));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc75));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc76));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc77));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc78));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc79));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx8, mc80));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc81));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc82));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc83));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc84));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc85));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc86));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc87));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc88));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx8, mc89));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc90));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc91));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc92));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc93));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc94));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc95));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc96));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc97));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx8, mc98));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc99));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc100));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc101));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc102));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc103));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc104));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc105));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc106));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx8, mc107));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc108));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc109));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc110));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc111));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc112));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc113));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc114));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc115));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx8, mc116));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc117));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc118));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc119));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc120));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc121));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc122));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc123));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc124));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx8, mc125));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc126));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc127));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc128));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc129));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc130));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc131));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc132));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc133));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx8, mc134));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc135));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc136));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc137));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc138));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc139));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc140));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc141));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc142));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx8, mc143));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx0, mc144));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx1, mc145));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx2, mc146));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx3, mc147));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx4, mc148));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx5, mc149));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx6, mc150));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx7, mc151));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx8, mc152));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx0, mc153));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx1, mc154));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx2, mc155));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx3, mc156));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx4, mc157));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx5, mc158));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx6, mc159));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx7, mc160));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx8, mc161))

#define setup_18(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0);\
                         msum16 = _mm_set_pd(0,0);\
                         msum17 = _mm_set_pd(0,0)

#define save_18() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  msum16 = _mm_hadd_pd(msum16, msum17);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14);\
                  _mm_storeu_pd(y + t1 + 16,msum16)

PetscErrorCode BSG_MatMult_18(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, mc128, mc129, mc130, mc131, mc132, mc133, mc134, mc135, mc136, mc137, mc138, mc139, mc140, mc141, mc142, mc143, mc144, mc145, mc146, mc147, mc148, mc149, mc150, mc151, mc152, mc153, mc154, mc155, mc156, mc157, mc158, mc159, mc160, mc161, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15, msum16, msum17;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_18(0);
            inline_18(3);
            inline_18(4);
            inline_18(5);
            inline_18(6);
            save_18();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_18(1);
            inline_18(2);
            inline_18(3);
            inline_18(4);
            inline_18(5);
            inline_18(6);
            save_18();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_18(lda3);
            inline_18(1);
            inline_18(2);
            inline_18(3);
            inline_18(4);
            inline_18(5);
            inline_18(6);
            save_18();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_18(lda2);
            inline_18(0);
            inline_18(1);
            inline_18(2);
            inline_18(3);
            inline_18(4);
            inline_18(5);
            inline_18(6);
            save_18();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_18(lda1 - lda2);
            inline_18(0);
            inline_18(1);
            inline_18(2);
            inline_18(3);
            inline_18(4);
            inline_18(5);
            save_18();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_18(lda1 - lda3);
            inline_18(0);
            inline_18(1);
            inline_18(2);
            inline_18(3);
            inline_18(4);
            save_18();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_18(lda1 - 1);
            inline_18(0);
            inline_18(1);
            inline_18(2);
            inline_18(3);
            save_18();
        }
    }

PetscFunctionReturn(0);
}

#define inline_20(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mx8 = _mm_loadu_pd(xt##l + t1 + 16);\
                     mx9 = _mm_loadu_pd(xt##l + t1 + 18);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     mc128 = _mm_loadu_pd(ct##l + t2 + 256);\
                     mc129 = _mm_loadu_pd(ct##l + t2 + 258);\
                     mc130 = _mm_loadu_pd(ct##l + t2 + 260);\
                     mc131 = _mm_loadu_pd(ct##l + t2 + 262);\
                     mc132 = _mm_loadu_pd(ct##l + t2 + 264);\
                     mc133 = _mm_loadu_pd(ct##l + t2 + 266);\
                     mc134 = _mm_loadu_pd(ct##l + t2 + 268);\
                     mc135 = _mm_loadu_pd(ct##l + t2 + 270);\
                     mc136 = _mm_loadu_pd(ct##l + t2 + 272);\
                     mc137 = _mm_loadu_pd(ct##l + t2 + 274);\
                     mc138 = _mm_loadu_pd(ct##l + t2 + 276);\
                     mc139 = _mm_loadu_pd(ct##l + t2 + 278);\
                     mc140 = _mm_loadu_pd(ct##l + t2 + 280);\
                     mc141 = _mm_loadu_pd(ct##l + t2 + 282);\
                     mc142 = _mm_loadu_pd(ct##l + t2 + 284);\
                     mc143 = _mm_loadu_pd(ct##l + t2 + 286);\
                     mc144 = _mm_loadu_pd(ct##l + t2 + 288);\
                     mc145 = _mm_loadu_pd(ct##l + t2 + 290);\
                     mc146 = _mm_loadu_pd(ct##l + t2 + 292);\
                     mc147 = _mm_loadu_pd(ct##l + t2 + 294);\
                     mc148 = _mm_loadu_pd(ct##l + t2 + 296);\
                     mc149 = _mm_loadu_pd(ct##l + t2 + 298);\
                     mc150 = _mm_loadu_pd(ct##l + t2 + 300);\
                     mc151 = _mm_loadu_pd(ct##l + t2 + 302);\
                     mc152 = _mm_loadu_pd(ct##l + t2 + 304);\
                     mc153 = _mm_loadu_pd(ct##l + t2 + 306);\
                     mc154 = _mm_loadu_pd(ct##l + t2 + 308);\
                     mc155 = _mm_loadu_pd(ct##l + t2 + 310);\
                     mc156 = _mm_loadu_pd(ct##l + t2 + 312);\
                     mc157 = _mm_loadu_pd(ct##l + t2 + 314);\
                     mc158 = _mm_loadu_pd(ct##l + t2 + 316);\
                     mc159 = _mm_loadu_pd(ct##l + t2 + 318);\
                     mc160 = _mm_loadu_pd(ct##l + t2 + 320);\
                     mc161 = _mm_loadu_pd(ct##l + t2 + 322);\
                     mc162 = _mm_loadu_pd(ct##l + t2 + 324);\
                     mc163 = _mm_loadu_pd(ct##l + t2 + 326);\
                     mc164 = _mm_loadu_pd(ct##l + t2 + 328);\
                     mc165 = _mm_loadu_pd(ct##l + t2 + 330);\
                     mc166 = _mm_loadu_pd(ct##l + t2 + 332);\
                     mc167 = _mm_loadu_pd(ct##l + t2 + 334);\
                     mc168 = _mm_loadu_pd(ct##l + t2 + 336);\
                     mc169 = _mm_loadu_pd(ct##l + t2 + 338);\
                     mc170 = _mm_loadu_pd(ct##l + t2 + 340);\
                     mc171 = _mm_loadu_pd(ct##l + t2 + 342);\
                     mc172 = _mm_loadu_pd(ct##l + t2 + 344);\
                     mc173 = _mm_loadu_pd(ct##l + t2 + 346);\
                     mc174 = _mm_loadu_pd(ct##l + t2 + 348);\
                     mc175 = _mm_loadu_pd(ct##l + t2 + 350);\
                     mc176 = _mm_loadu_pd(ct##l + t2 + 352);\
                     mc177 = _mm_loadu_pd(ct##l + t2 + 354);\
                     mc178 = _mm_loadu_pd(ct##l + t2 + 356);\
                     mc179 = _mm_loadu_pd(ct##l + t2 + 358);\
                     mc180 = _mm_loadu_pd(ct##l + t2 + 360);\
                     mc181 = _mm_loadu_pd(ct##l + t2 + 362);\
                     mc182 = _mm_loadu_pd(ct##l + t2 + 364);\
                     mc183 = _mm_loadu_pd(ct##l + t2 + 366);\
                     mc184 = _mm_loadu_pd(ct##l + t2 + 368);\
                     mc185 = _mm_loadu_pd(ct##l + t2 + 370);\
                     mc186 = _mm_loadu_pd(ct##l + t2 + 372);\
                     mc187 = _mm_loadu_pd(ct##l + t2 + 374);\
                     mc188 = _mm_loadu_pd(ct##l + t2 + 376);\
                     mc189 = _mm_loadu_pd(ct##l + t2 + 378);\
                     mc190 = _mm_loadu_pd(ct##l + t2 + 380);\
                     mc191 = _mm_loadu_pd(ct##l + t2 + 382);\
                     mc192 = _mm_loadu_pd(ct##l + t2 + 384);\
                     mc193 = _mm_loadu_pd(ct##l + t2 + 386);\
                     mc194 = _mm_loadu_pd(ct##l + t2 + 388);\
                     mc195 = _mm_loadu_pd(ct##l + t2 + 390);\
                     mc196 = _mm_loadu_pd(ct##l + t2 + 392);\
                     mc197 = _mm_loadu_pd(ct##l + t2 + 394);\
                     mc198 = _mm_loadu_pd(ct##l + t2 + 396);\
                     mc199 = _mm_loadu_pd(ct##l + t2 + 398);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx8, mc8));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx9, mc9));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc10));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc11));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc12));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc13));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc15));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc16));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc17));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx8, mc18));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx9, mc19));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc20));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc21));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc22));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc23));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc24));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc25));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc26));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc27));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx8, mc28));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx9, mc29));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc30));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc31));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc32));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc33));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc34));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc35));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc36));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc37));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx8, mc38));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx9, mc39));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc40));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc41));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc42));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc43));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc44));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc45));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc46));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc47));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx8, mc48));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx9, mc49));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc50));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc51));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc52));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc53));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc54));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc55));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc56));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc57));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx8, mc58));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx9, mc59));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc60));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc61));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc62));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc63));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc64));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc65));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc66));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc67));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx8, mc68));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx9, mc69));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc70));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc71));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc72));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc73));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc74));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc75));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc76));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc77));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx8, mc78));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx9, mc79));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc80));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc81));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc82));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc83));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc84));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc85));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc86));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc87));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx8, mc88));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx9, mc89));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc90));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc91));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc92));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc93));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc94));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc95));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc96));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc97));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx8, mc98));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx9, mc99));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc100));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc101));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc102));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc103));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc104));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc105));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc106));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc107));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx8, mc108));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx9, mc109));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc110));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc111));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc112));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc113));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc114));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc115));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc116));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc117));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx8, mc118));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx9, mc119));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc120));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc121));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc122));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc123));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc124));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc125));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc126));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc127));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx8, mc128));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx9, mc129));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc130));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc131));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc132));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc133));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc134));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc135));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc136));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc137));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx8, mc138));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx9, mc139));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc140));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc141));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc142));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc143));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc144));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc145));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc146));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc147));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx8, mc148));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx9, mc149));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc150));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc151));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc152));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc153));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc154));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc155));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc156));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc157));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx8, mc158));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx9, mc159));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx0, mc160));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx1, mc161));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx2, mc162));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx3, mc163));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx4, mc164));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx5, mc165));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx6, mc166));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx7, mc167));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx8, mc168));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx9, mc169));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx0, mc170));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx1, mc171));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx2, mc172));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx3, mc173));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx4, mc174));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx5, mc175));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx6, mc176));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx7, mc177));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx8, mc178));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx9, mc179));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx0, mc180));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx1, mc181));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx2, mc182));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx3, mc183));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx4, mc184));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx5, mc185));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx6, mc186));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx7, mc187));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx8, mc188));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx9, mc189));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx0, mc190));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx1, mc191));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx2, mc192));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx3, mc193));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx4, mc194));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx5, mc195));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx6, mc196));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx7, mc197));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx8, mc198));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx9, mc199))

#define setup_20(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0);\
                         msum16 = _mm_set_pd(0,0);\
                         msum17 = _mm_set_pd(0,0);\
                         msum18 = _mm_set_pd(0,0);\
                         msum19 = _mm_set_pd(0,0)

#define save_20() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  msum16 = _mm_hadd_pd(msum16, msum17);\
                  msum18 = _mm_hadd_pd(msum18, msum19);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14);\
                  _mm_storeu_pd(y + t1 + 16,msum16);\
                  _mm_storeu_pd(y + t1 + 18,msum18)

PetscErrorCode BSG_MatMult_20(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8, mx9, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, mc128, mc129, mc130, mc131, mc132, mc133, mc134, mc135, mc136, mc137, mc138, mc139, mc140, mc141, mc142, mc143, mc144, mc145, mc146, mc147, mc148, mc149, mc150, mc151, mc152, mc153, mc154, mc155, mc156, mc157, mc158, mc159, mc160, mc161, mc162, mc163, mc164, mc165, mc166, mc167, mc168, mc169, mc170, mc171, mc172, mc173, mc174, mc175, mc176, mc177, mc178, mc179, mc180, mc181, mc182, mc183, mc184, mc185, mc186, mc187, mc188, mc189, mc190, mc191, mc192, mc193, mc194, mc195, mc196, mc197, mc198, mc199, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15, msum16, msum17, msum18, msum19;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_20(0);
            inline_20(3);
            inline_20(4);
            inline_20(5);
            inline_20(6);
            save_20();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_20(1);
            inline_20(2);
            inline_20(3);
            inline_20(4);
            inline_20(5);
            inline_20(6);
            save_20();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_20(lda3);
            inline_20(1);
            inline_20(2);
            inline_20(3);
            inline_20(4);
            inline_20(5);
            inline_20(6);
            save_20();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_20(lda2);
            inline_20(0);
            inline_20(1);
            inline_20(2);
            inline_20(3);
            inline_20(4);
            inline_20(5);
            inline_20(6);
            save_20();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_20(lda1 - lda2);
            inline_20(0);
            inline_20(1);
            inline_20(2);
            inline_20(3);
            inline_20(4);
            inline_20(5);
            save_20();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_20(lda1 - lda3);
            inline_20(0);
            inline_20(1);
            inline_20(2);
            inline_20(3);
            inline_20(4);
            save_20();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_20(lda1 - 1);
            inline_20(0);
            inline_20(1);
            inline_20(2);
            inline_20(3);
            save_20();
        }
    }

PetscFunctionReturn(0);
}

#define inline_22(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mx8 = _mm_loadu_pd(xt##l + t1 + 16);\
                     mx9 = _mm_loadu_pd(xt##l + t1 + 18);\
                     mx10 = _mm_loadu_pd(xt##l + t1 + 20);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     mc128 = _mm_loadu_pd(ct##l + t2 + 256);\
                     mc129 = _mm_loadu_pd(ct##l + t2 + 258);\
                     mc130 = _mm_loadu_pd(ct##l + t2 + 260);\
                     mc131 = _mm_loadu_pd(ct##l + t2 + 262);\
                     mc132 = _mm_loadu_pd(ct##l + t2 + 264);\
                     mc133 = _mm_loadu_pd(ct##l + t2 + 266);\
                     mc134 = _mm_loadu_pd(ct##l + t2 + 268);\
                     mc135 = _mm_loadu_pd(ct##l + t2 + 270);\
                     mc136 = _mm_loadu_pd(ct##l + t2 + 272);\
                     mc137 = _mm_loadu_pd(ct##l + t2 + 274);\
                     mc138 = _mm_loadu_pd(ct##l + t2 + 276);\
                     mc139 = _mm_loadu_pd(ct##l + t2 + 278);\
                     mc140 = _mm_loadu_pd(ct##l + t2 + 280);\
                     mc141 = _mm_loadu_pd(ct##l + t2 + 282);\
                     mc142 = _mm_loadu_pd(ct##l + t2 + 284);\
                     mc143 = _mm_loadu_pd(ct##l + t2 + 286);\
                     mc144 = _mm_loadu_pd(ct##l + t2 + 288);\
                     mc145 = _mm_loadu_pd(ct##l + t2 + 290);\
                     mc146 = _mm_loadu_pd(ct##l + t2 + 292);\
                     mc147 = _mm_loadu_pd(ct##l + t2 + 294);\
                     mc148 = _mm_loadu_pd(ct##l + t2 + 296);\
                     mc149 = _mm_loadu_pd(ct##l + t2 + 298);\
                     mc150 = _mm_loadu_pd(ct##l + t2 + 300);\
                     mc151 = _mm_loadu_pd(ct##l + t2 + 302);\
                     mc152 = _mm_loadu_pd(ct##l + t2 + 304);\
                     mc153 = _mm_loadu_pd(ct##l + t2 + 306);\
                     mc154 = _mm_loadu_pd(ct##l + t2 + 308);\
                     mc155 = _mm_loadu_pd(ct##l + t2 + 310);\
                     mc156 = _mm_loadu_pd(ct##l + t2 + 312);\
                     mc157 = _mm_loadu_pd(ct##l + t2 + 314);\
                     mc158 = _mm_loadu_pd(ct##l + t2 + 316);\
                     mc159 = _mm_loadu_pd(ct##l + t2 + 318);\
                     mc160 = _mm_loadu_pd(ct##l + t2 + 320);\
                     mc161 = _mm_loadu_pd(ct##l + t2 + 322);\
                     mc162 = _mm_loadu_pd(ct##l + t2 + 324);\
                     mc163 = _mm_loadu_pd(ct##l + t2 + 326);\
                     mc164 = _mm_loadu_pd(ct##l + t2 + 328);\
                     mc165 = _mm_loadu_pd(ct##l + t2 + 330);\
                     mc166 = _mm_loadu_pd(ct##l + t2 + 332);\
                     mc167 = _mm_loadu_pd(ct##l + t2 + 334);\
                     mc168 = _mm_loadu_pd(ct##l + t2 + 336);\
                     mc169 = _mm_loadu_pd(ct##l + t2 + 338);\
                     mc170 = _mm_loadu_pd(ct##l + t2 + 340);\
                     mc171 = _mm_loadu_pd(ct##l + t2 + 342);\
                     mc172 = _mm_loadu_pd(ct##l + t2 + 344);\
                     mc173 = _mm_loadu_pd(ct##l + t2 + 346);\
                     mc174 = _mm_loadu_pd(ct##l + t2 + 348);\
                     mc175 = _mm_loadu_pd(ct##l + t2 + 350);\
                     mc176 = _mm_loadu_pd(ct##l + t2 + 352);\
                     mc177 = _mm_loadu_pd(ct##l + t2 + 354);\
                     mc178 = _mm_loadu_pd(ct##l + t2 + 356);\
                     mc179 = _mm_loadu_pd(ct##l + t2 + 358);\
                     mc180 = _mm_loadu_pd(ct##l + t2 + 360);\
                     mc181 = _mm_loadu_pd(ct##l + t2 + 362);\
                     mc182 = _mm_loadu_pd(ct##l + t2 + 364);\
                     mc183 = _mm_loadu_pd(ct##l + t2 + 366);\
                     mc184 = _mm_loadu_pd(ct##l + t2 + 368);\
                     mc185 = _mm_loadu_pd(ct##l + t2 + 370);\
                     mc186 = _mm_loadu_pd(ct##l + t2 + 372);\
                     mc187 = _mm_loadu_pd(ct##l + t2 + 374);\
                     mc188 = _mm_loadu_pd(ct##l + t2 + 376);\
                     mc189 = _mm_loadu_pd(ct##l + t2 + 378);\
                     mc190 = _mm_loadu_pd(ct##l + t2 + 380);\
                     mc191 = _mm_loadu_pd(ct##l + t2 + 382);\
                     mc192 = _mm_loadu_pd(ct##l + t2 + 384);\
                     mc193 = _mm_loadu_pd(ct##l + t2 + 386);\
                     mc194 = _mm_loadu_pd(ct##l + t2 + 388);\
                     mc195 = _mm_loadu_pd(ct##l + t2 + 390);\
                     mc196 = _mm_loadu_pd(ct##l + t2 + 392);\
                     mc197 = _mm_loadu_pd(ct##l + t2 + 394);\
                     mc198 = _mm_loadu_pd(ct##l + t2 + 396);\
                     mc199 = _mm_loadu_pd(ct##l + t2 + 398);\
                     mc200 = _mm_loadu_pd(ct##l + t2 + 400);\
                     mc201 = _mm_loadu_pd(ct##l + t2 + 402);\
                     mc202 = _mm_loadu_pd(ct##l + t2 + 404);\
                     mc203 = _mm_loadu_pd(ct##l + t2 + 406);\
                     mc204 = _mm_loadu_pd(ct##l + t2 + 408);\
                     mc205 = _mm_loadu_pd(ct##l + t2 + 410);\
                     mc206 = _mm_loadu_pd(ct##l + t2 + 412);\
                     mc207 = _mm_loadu_pd(ct##l + t2 + 414);\
                     mc208 = _mm_loadu_pd(ct##l + t2 + 416);\
                     mc209 = _mm_loadu_pd(ct##l + t2 + 418);\
                     mc210 = _mm_loadu_pd(ct##l + t2 + 420);\
                     mc211 = _mm_loadu_pd(ct##l + t2 + 422);\
                     mc212 = _mm_loadu_pd(ct##l + t2 + 424);\
                     mc213 = _mm_loadu_pd(ct##l + t2 + 426);\
                     mc214 = _mm_loadu_pd(ct##l + t2 + 428);\
                     mc215 = _mm_loadu_pd(ct##l + t2 + 430);\
                     mc216 = _mm_loadu_pd(ct##l + t2 + 432);\
                     mc217 = _mm_loadu_pd(ct##l + t2 + 434);\
                     mc218 = _mm_loadu_pd(ct##l + t2 + 436);\
                     mc219 = _mm_loadu_pd(ct##l + t2 + 438);\
                     mc220 = _mm_loadu_pd(ct##l + t2 + 440);\
                     mc221 = _mm_loadu_pd(ct##l + t2 + 442);\
                     mc222 = _mm_loadu_pd(ct##l + t2 + 444);\
                     mc223 = _mm_loadu_pd(ct##l + t2 + 446);\
                     mc224 = _mm_loadu_pd(ct##l + t2 + 448);\
                     mc225 = _mm_loadu_pd(ct##l + t2 + 450);\
                     mc226 = _mm_loadu_pd(ct##l + t2 + 452);\
                     mc227 = _mm_loadu_pd(ct##l + t2 + 454);\
                     mc228 = _mm_loadu_pd(ct##l + t2 + 456);\
                     mc229 = _mm_loadu_pd(ct##l + t2 + 458);\
                     mc230 = _mm_loadu_pd(ct##l + t2 + 460);\
                     mc231 = _mm_loadu_pd(ct##l + t2 + 462);\
                     mc232 = _mm_loadu_pd(ct##l + t2 + 464);\
                     mc233 = _mm_loadu_pd(ct##l + t2 + 466);\
                     mc234 = _mm_loadu_pd(ct##l + t2 + 468);\
                     mc235 = _mm_loadu_pd(ct##l + t2 + 470);\
                     mc236 = _mm_loadu_pd(ct##l + t2 + 472);\
                     mc237 = _mm_loadu_pd(ct##l + t2 + 474);\
                     mc238 = _mm_loadu_pd(ct##l + t2 + 476);\
                     mc239 = _mm_loadu_pd(ct##l + t2 + 478);\
                     mc240 = _mm_loadu_pd(ct##l + t2 + 480);\
                     mc241 = _mm_loadu_pd(ct##l + t2 + 482);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx8, mc8));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx9, mc9));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx10, mc10));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc11));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc12));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc13));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc15));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc16));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc17));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc18));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx8, mc19));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx9, mc20));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx10, mc21));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc22));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc23));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc24));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc25));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc26));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc27));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc28));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc29));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx8, mc30));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx9, mc31));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx10, mc32));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc33));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc34));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc35));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc36));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc37));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc38));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc39));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc40));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx8, mc41));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx9, mc42));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx10, mc43));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc44));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc45));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc46));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc47));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc48));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc49));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc50));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc51));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx8, mc52));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx9, mc53));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx10, mc54));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc55));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc56));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc57));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc58));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc59));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc60));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc61));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc62));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx8, mc63));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx9, mc64));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx10, mc65));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc66));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc67));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc68));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc69));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc70));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc71));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc72));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc73));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx8, mc74));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx9, mc75));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx10, mc76));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc77));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc78));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc79));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc80));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc81));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc82));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc83));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc84));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx8, mc85));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx9, mc86));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx10, mc87));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc88));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc89));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc90));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc91));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc92));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc93));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc94));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc95));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx8, mc96));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx9, mc97));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx10, mc98));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc99));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc100));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc101));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc102));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc103));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc104));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc105));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc106));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx8, mc107));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx9, mc108));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx10, mc109));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc110));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc111));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc112));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc113));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc114));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc115));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc116));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc117));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx8, mc118));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx9, mc119));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx10, mc120));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc121));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc122));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc123));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc124));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc125));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc126));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc127));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc128));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx8, mc129));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx9, mc130));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx10, mc131));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc132));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc133));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc134));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc135));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc136));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc137));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc138));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc139));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx8, mc140));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx9, mc141));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx10, mc142));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc143));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc144));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc145));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc146));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc147));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc148));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc149));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc150));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx8, mc151));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx9, mc152));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx10, mc153));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc154));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc155));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc156));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc157));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc158));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc159));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc160));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc161));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx8, mc162));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx9, mc163));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx10, mc164));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc165));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc166));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc167));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc168));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc169));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc170));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc171));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc172));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx8, mc173));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx9, mc174));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx10, mc175));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx0, mc176));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx1, mc177));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx2, mc178));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx3, mc179));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx4, mc180));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx5, mc181));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx6, mc182));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx7, mc183));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx8, mc184));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx9, mc185));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx10, mc186));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx0, mc187));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx1, mc188));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx2, mc189));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx3, mc190));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx4, mc191));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx5, mc192));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx6, mc193));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx7, mc194));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx8, mc195));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx9, mc196));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx10, mc197));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx0, mc198));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx1, mc199));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx2, mc200));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx3, mc201));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx4, mc202));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx5, mc203));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx6, mc204));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx7, mc205));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx8, mc206));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx9, mc207));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx10, mc208));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx0, mc209));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx1, mc210));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx2, mc211));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx3, mc212));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx4, mc213));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx5, mc214));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx6, mc215));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx7, mc216));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx8, mc217));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx9, mc218));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx10, mc219));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx0, mc220));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx1, mc221));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx2, mc222));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx3, mc223));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx4, mc224));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx5, mc225));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx6, mc226));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx7, mc227));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx8, mc228));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx9, mc229));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx10, mc230));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx0, mc231));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx1, mc232));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx2, mc233));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx3, mc234));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx4, mc235));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx5, mc236));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx6, mc237));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx7, mc238));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx8, mc239));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx9, mc240));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx10, mc241))

#define setup_22(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0);\
                         msum16 = _mm_set_pd(0,0);\
                         msum17 = _mm_set_pd(0,0);\
                         msum18 = _mm_set_pd(0,0);\
                         msum19 = _mm_set_pd(0,0);\
                         msum20 = _mm_set_pd(0,0);\
                         msum21 = _mm_set_pd(0,0)

#define save_22() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  msum16 = _mm_hadd_pd(msum16, msum17);\
                  msum18 = _mm_hadd_pd(msum18, msum19);\
                  msum20 = _mm_hadd_pd(msum20, msum21);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14);\
                  _mm_storeu_pd(y + t1 + 16,msum16);\
                  _mm_storeu_pd(y + t1 + 18,msum18);\
                  _mm_storeu_pd(y + t1 + 20,msum20)

PetscErrorCode BSG_MatMult_22(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8, mx9, mx10, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, mc128, mc129, mc130, mc131, mc132, mc133, mc134, mc135, mc136, mc137, mc138, mc139, mc140, mc141, mc142, mc143, mc144, mc145, mc146, mc147, mc148, mc149, mc150, mc151, mc152, mc153, mc154, mc155, mc156, mc157, mc158, mc159, mc160, mc161, mc162, mc163, mc164, mc165, mc166, mc167, mc168, mc169, mc170, mc171, mc172, mc173, mc174, mc175, mc176, mc177, mc178, mc179, mc180, mc181, mc182, mc183, mc184, mc185, mc186, mc187, mc188, mc189, mc190, mc191, mc192, mc193, mc194, mc195, mc196, mc197, mc198, mc199, mc200, mc201, mc202, mc203, mc204, mc205, mc206, mc207, mc208, mc209, mc210, mc211, mc212, mc213, mc214, mc215, mc216, mc217, mc218, mc219, mc220, mc221, mc222, mc223, mc224, mc225, mc226, mc227, mc228, mc229, mc230, mc231, mc232, mc233, mc234, mc235, mc236, mc237, mc238, mc239, mc240, mc241, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15, msum16, msum17, msum18, msum19, msum20, msum21;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_22(0);
            inline_22(3);
            inline_22(4);
            inline_22(5);
            inline_22(6);
            save_22();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_22(1);
            inline_22(2);
            inline_22(3);
            inline_22(4);
            inline_22(5);
            inline_22(6);
            save_22();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_22(lda3);
            inline_22(1);
            inline_22(2);
            inline_22(3);
            inline_22(4);
            inline_22(5);
            inline_22(6);
            save_22();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_22(lda2);
            inline_22(0);
            inline_22(1);
            inline_22(2);
            inline_22(3);
            inline_22(4);
            inline_22(5);
            inline_22(6);
            save_22();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_22(lda1 - lda2);
            inline_22(0);
            inline_22(1);
            inline_22(2);
            inline_22(3);
            inline_22(4);
            inline_22(5);
            save_22();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_22(lda1 - lda3);
            inline_22(0);
            inline_22(1);
            inline_22(2);
            inline_22(3);
            inline_22(4);
            save_22();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_22(lda1 - 1);
            inline_22(0);
            inline_22(1);
            inline_22(2);
            inline_22(3);
            save_22();
        }
    }

PetscFunctionReturn(0);
}

#define inline_24(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mx8 = _mm_loadu_pd(xt##l + t1 + 16);\
                     mx9 = _mm_loadu_pd(xt##l + t1 + 18);\
                     mx10 = _mm_loadu_pd(xt##l + t1 + 20);\
                     mx11 = _mm_loadu_pd(xt##l + t1 + 22);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     mc128 = _mm_loadu_pd(ct##l + t2 + 256);\
                     mc129 = _mm_loadu_pd(ct##l + t2 + 258);\
                     mc130 = _mm_loadu_pd(ct##l + t2 + 260);\
                     mc131 = _mm_loadu_pd(ct##l + t2 + 262);\
                     mc132 = _mm_loadu_pd(ct##l + t2 + 264);\
                     mc133 = _mm_loadu_pd(ct##l + t2 + 266);\
                     mc134 = _mm_loadu_pd(ct##l + t2 + 268);\
                     mc135 = _mm_loadu_pd(ct##l + t2 + 270);\
                     mc136 = _mm_loadu_pd(ct##l + t2 + 272);\
                     mc137 = _mm_loadu_pd(ct##l + t2 + 274);\
                     mc138 = _mm_loadu_pd(ct##l + t2 + 276);\
                     mc139 = _mm_loadu_pd(ct##l + t2 + 278);\
                     mc140 = _mm_loadu_pd(ct##l + t2 + 280);\
                     mc141 = _mm_loadu_pd(ct##l + t2 + 282);\
                     mc142 = _mm_loadu_pd(ct##l + t2 + 284);\
                     mc143 = _mm_loadu_pd(ct##l + t2 + 286);\
                     mc144 = _mm_loadu_pd(ct##l + t2 + 288);\
                     mc145 = _mm_loadu_pd(ct##l + t2 + 290);\
                     mc146 = _mm_loadu_pd(ct##l + t2 + 292);\
                     mc147 = _mm_loadu_pd(ct##l + t2 + 294);\
                     mc148 = _mm_loadu_pd(ct##l + t2 + 296);\
                     mc149 = _mm_loadu_pd(ct##l + t2 + 298);\
                     mc150 = _mm_loadu_pd(ct##l + t2 + 300);\
                     mc151 = _mm_loadu_pd(ct##l + t2 + 302);\
                     mc152 = _mm_loadu_pd(ct##l + t2 + 304);\
                     mc153 = _mm_loadu_pd(ct##l + t2 + 306);\
                     mc154 = _mm_loadu_pd(ct##l + t2 + 308);\
                     mc155 = _mm_loadu_pd(ct##l + t2 + 310);\
                     mc156 = _mm_loadu_pd(ct##l + t2 + 312);\
                     mc157 = _mm_loadu_pd(ct##l + t2 + 314);\
                     mc158 = _mm_loadu_pd(ct##l + t2 + 316);\
                     mc159 = _mm_loadu_pd(ct##l + t2 + 318);\
                     mc160 = _mm_loadu_pd(ct##l + t2 + 320);\
                     mc161 = _mm_loadu_pd(ct##l + t2 + 322);\
                     mc162 = _mm_loadu_pd(ct##l + t2 + 324);\
                     mc163 = _mm_loadu_pd(ct##l + t2 + 326);\
                     mc164 = _mm_loadu_pd(ct##l + t2 + 328);\
                     mc165 = _mm_loadu_pd(ct##l + t2 + 330);\
                     mc166 = _mm_loadu_pd(ct##l + t2 + 332);\
                     mc167 = _mm_loadu_pd(ct##l + t2 + 334);\
                     mc168 = _mm_loadu_pd(ct##l + t2 + 336);\
                     mc169 = _mm_loadu_pd(ct##l + t2 + 338);\
                     mc170 = _mm_loadu_pd(ct##l + t2 + 340);\
                     mc171 = _mm_loadu_pd(ct##l + t2 + 342);\
                     mc172 = _mm_loadu_pd(ct##l + t2 + 344);\
                     mc173 = _mm_loadu_pd(ct##l + t2 + 346);\
                     mc174 = _mm_loadu_pd(ct##l + t2 + 348);\
                     mc175 = _mm_loadu_pd(ct##l + t2 + 350);\
                     mc176 = _mm_loadu_pd(ct##l + t2 + 352);\
                     mc177 = _mm_loadu_pd(ct##l + t2 + 354);\
                     mc178 = _mm_loadu_pd(ct##l + t2 + 356);\
                     mc179 = _mm_loadu_pd(ct##l + t2 + 358);\
                     mc180 = _mm_loadu_pd(ct##l + t2 + 360);\
                     mc181 = _mm_loadu_pd(ct##l + t2 + 362);\
                     mc182 = _mm_loadu_pd(ct##l + t2 + 364);\
                     mc183 = _mm_loadu_pd(ct##l + t2 + 366);\
                     mc184 = _mm_loadu_pd(ct##l + t2 + 368);\
                     mc185 = _mm_loadu_pd(ct##l + t2 + 370);\
                     mc186 = _mm_loadu_pd(ct##l + t2 + 372);\
                     mc187 = _mm_loadu_pd(ct##l + t2 + 374);\
                     mc188 = _mm_loadu_pd(ct##l + t2 + 376);\
                     mc189 = _mm_loadu_pd(ct##l + t2 + 378);\
                     mc190 = _mm_loadu_pd(ct##l + t2 + 380);\
                     mc191 = _mm_loadu_pd(ct##l + t2 + 382);\
                     mc192 = _mm_loadu_pd(ct##l + t2 + 384);\
                     mc193 = _mm_loadu_pd(ct##l + t2 + 386);\
                     mc194 = _mm_loadu_pd(ct##l + t2 + 388);\
                     mc195 = _mm_loadu_pd(ct##l + t2 + 390);\
                     mc196 = _mm_loadu_pd(ct##l + t2 + 392);\
                     mc197 = _mm_loadu_pd(ct##l + t2 + 394);\
                     mc198 = _mm_loadu_pd(ct##l + t2 + 396);\
                     mc199 = _mm_loadu_pd(ct##l + t2 + 398);\
                     mc200 = _mm_loadu_pd(ct##l + t2 + 400);\
                     mc201 = _mm_loadu_pd(ct##l + t2 + 402);\
                     mc202 = _mm_loadu_pd(ct##l + t2 + 404);\
                     mc203 = _mm_loadu_pd(ct##l + t2 + 406);\
                     mc204 = _mm_loadu_pd(ct##l + t2 + 408);\
                     mc205 = _mm_loadu_pd(ct##l + t2 + 410);\
                     mc206 = _mm_loadu_pd(ct##l + t2 + 412);\
                     mc207 = _mm_loadu_pd(ct##l + t2 + 414);\
                     mc208 = _mm_loadu_pd(ct##l + t2 + 416);\
                     mc209 = _mm_loadu_pd(ct##l + t2 + 418);\
                     mc210 = _mm_loadu_pd(ct##l + t2 + 420);\
                     mc211 = _mm_loadu_pd(ct##l + t2 + 422);\
                     mc212 = _mm_loadu_pd(ct##l + t2 + 424);\
                     mc213 = _mm_loadu_pd(ct##l + t2 + 426);\
                     mc214 = _mm_loadu_pd(ct##l + t2 + 428);\
                     mc215 = _mm_loadu_pd(ct##l + t2 + 430);\
                     mc216 = _mm_loadu_pd(ct##l + t2 + 432);\
                     mc217 = _mm_loadu_pd(ct##l + t2 + 434);\
                     mc218 = _mm_loadu_pd(ct##l + t2 + 436);\
                     mc219 = _mm_loadu_pd(ct##l + t2 + 438);\
                     mc220 = _mm_loadu_pd(ct##l + t2 + 440);\
                     mc221 = _mm_loadu_pd(ct##l + t2 + 442);\
                     mc222 = _mm_loadu_pd(ct##l + t2 + 444);\
                     mc223 = _mm_loadu_pd(ct##l + t2 + 446);\
                     mc224 = _mm_loadu_pd(ct##l + t2 + 448);\
                     mc225 = _mm_loadu_pd(ct##l + t2 + 450);\
                     mc226 = _mm_loadu_pd(ct##l + t2 + 452);\
                     mc227 = _mm_loadu_pd(ct##l + t2 + 454);\
                     mc228 = _mm_loadu_pd(ct##l + t2 + 456);\
                     mc229 = _mm_loadu_pd(ct##l + t2 + 458);\
                     mc230 = _mm_loadu_pd(ct##l + t2 + 460);\
                     mc231 = _mm_loadu_pd(ct##l + t2 + 462);\
                     mc232 = _mm_loadu_pd(ct##l + t2 + 464);\
                     mc233 = _mm_loadu_pd(ct##l + t2 + 466);\
                     mc234 = _mm_loadu_pd(ct##l + t2 + 468);\
                     mc235 = _mm_loadu_pd(ct##l + t2 + 470);\
                     mc236 = _mm_loadu_pd(ct##l + t2 + 472);\
                     mc237 = _mm_loadu_pd(ct##l + t2 + 474);\
                     mc238 = _mm_loadu_pd(ct##l + t2 + 476);\
                     mc239 = _mm_loadu_pd(ct##l + t2 + 478);\
                     mc240 = _mm_loadu_pd(ct##l + t2 + 480);\
                     mc241 = _mm_loadu_pd(ct##l + t2 + 482);\
                     mc242 = _mm_loadu_pd(ct##l + t2 + 484);\
                     mc243 = _mm_loadu_pd(ct##l + t2 + 486);\
                     mc244 = _mm_loadu_pd(ct##l + t2 + 488);\
                     mc245 = _mm_loadu_pd(ct##l + t2 + 490);\
                     mc246 = _mm_loadu_pd(ct##l + t2 + 492);\
                     mc247 = _mm_loadu_pd(ct##l + t2 + 494);\
                     mc248 = _mm_loadu_pd(ct##l + t2 + 496);\
                     mc249 = _mm_loadu_pd(ct##l + t2 + 498);\
                     mc250 = _mm_loadu_pd(ct##l + t2 + 500);\
                     mc251 = _mm_loadu_pd(ct##l + t2 + 502);\
                     mc252 = _mm_loadu_pd(ct##l + t2 + 504);\
                     mc253 = _mm_loadu_pd(ct##l + t2 + 506);\
                     mc254 = _mm_loadu_pd(ct##l + t2 + 508);\
                     mc255 = _mm_loadu_pd(ct##l + t2 + 510);\
                     mc256 = _mm_loadu_pd(ct##l + t2 + 512);\
                     mc257 = _mm_loadu_pd(ct##l + t2 + 514);\
                     mc258 = _mm_loadu_pd(ct##l + t2 + 516);\
                     mc259 = _mm_loadu_pd(ct##l + t2 + 518);\
                     mc260 = _mm_loadu_pd(ct##l + t2 + 520);\
                     mc261 = _mm_loadu_pd(ct##l + t2 + 522);\
                     mc262 = _mm_loadu_pd(ct##l + t2 + 524);\
                     mc263 = _mm_loadu_pd(ct##l + t2 + 526);\
                     mc264 = _mm_loadu_pd(ct##l + t2 + 528);\
                     mc265 = _mm_loadu_pd(ct##l + t2 + 530);\
                     mc266 = _mm_loadu_pd(ct##l + t2 + 532);\
                     mc267 = _mm_loadu_pd(ct##l + t2 + 534);\
                     mc268 = _mm_loadu_pd(ct##l + t2 + 536);\
                     mc269 = _mm_loadu_pd(ct##l + t2 + 538);\
                     mc270 = _mm_loadu_pd(ct##l + t2 + 540);\
                     mc271 = _mm_loadu_pd(ct##l + t2 + 542);\
                     mc272 = _mm_loadu_pd(ct##l + t2 + 544);\
                     mc273 = _mm_loadu_pd(ct##l + t2 + 546);\
                     mc274 = _mm_loadu_pd(ct##l + t2 + 548);\
                     mc275 = _mm_loadu_pd(ct##l + t2 + 550);\
                     mc276 = _mm_loadu_pd(ct##l + t2 + 552);\
                     mc277 = _mm_loadu_pd(ct##l + t2 + 554);\
                     mc278 = _mm_loadu_pd(ct##l + t2 + 556);\
                     mc279 = _mm_loadu_pd(ct##l + t2 + 558);\
                     mc280 = _mm_loadu_pd(ct##l + t2 + 560);\
                     mc281 = _mm_loadu_pd(ct##l + t2 + 562);\
                     mc282 = _mm_loadu_pd(ct##l + t2 + 564);\
                     mc283 = _mm_loadu_pd(ct##l + t2 + 566);\
                     mc284 = _mm_loadu_pd(ct##l + t2 + 568);\
                     mc285 = _mm_loadu_pd(ct##l + t2 + 570);\
                     mc286 = _mm_loadu_pd(ct##l + t2 + 572);\
                     mc287 = _mm_loadu_pd(ct##l + t2 + 574);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx8, mc8));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx9, mc9));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx10, mc10));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx11, mc11));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc12));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc13));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc15));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc16));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc17));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc18));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc19));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx8, mc20));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx9, mc21));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx10, mc22));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx11, mc23));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc24));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc25));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc26));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc27));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc28));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc29));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc30));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc31));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx8, mc32));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx9, mc33));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx10, mc34));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx11, mc35));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc36));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc37));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc38));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc39));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc40));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc41));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc42));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc43));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx8, mc44));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx9, mc45));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx10, mc46));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx11, mc47));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc48));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc49));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc50));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc51));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc52));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc53));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc54));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc55));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx8, mc56));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx9, mc57));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx10, mc58));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx11, mc59));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc60));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc61));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc62));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc63));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc64));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc65));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc66));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc67));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx8, mc68));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx9, mc69));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx10, mc70));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx11, mc71));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc72));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc73));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc74));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc75));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc76));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc77));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc78));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc79));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx8, mc80));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx9, mc81));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx10, mc82));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx11, mc83));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc84));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc85));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc86));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc87));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc88));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc89));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc90));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc91));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx8, mc92));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx9, mc93));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx10, mc94));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx11, mc95));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc96));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc97));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc98));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc99));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc100));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc101));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc102));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc103));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx8, mc104));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx9, mc105));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx10, mc106));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx11, mc107));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc108));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc109));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc110));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc111));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc112));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc113));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc114));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc115));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx8, mc116));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx9, mc117));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx10, mc118));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx11, mc119));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc120));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc121));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc122));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc123));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc124));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc125));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc126));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc127));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx8, mc128));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx9, mc129));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx10, mc130));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx11, mc131));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc132));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc133));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc134));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc135));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc136));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc137));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc138));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc139));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx8, mc140));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx9, mc141));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx10, mc142));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx11, mc143));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc144));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc145));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc146));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc147));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc148));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc149));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc150));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc151));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx8, mc152));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx9, mc153));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx10, mc154));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx11, mc155));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc156));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc157));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc158));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc159));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc160));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc161));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc162));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc163));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx8, mc164));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx9, mc165));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx10, mc166));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx11, mc167));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc168));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc169));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc170));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc171));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc172));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc173));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc174));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc175));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx8, mc176));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx9, mc177));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx10, mc178));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx11, mc179));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc180));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc181));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc182));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc183));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc184));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc185));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc186));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc187));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx8, mc188));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx9, mc189));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx10, mc190));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx11, mc191));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx0, mc192));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx1, mc193));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx2, mc194));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx3, mc195));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx4, mc196));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx5, mc197));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx6, mc198));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx7, mc199));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx8, mc200));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx9, mc201));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx10, mc202));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx11, mc203));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx0, mc204));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx1, mc205));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx2, mc206));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx3, mc207));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx4, mc208));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx5, mc209));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx6, mc210));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx7, mc211));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx8, mc212));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx9, mc213));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx10, mc214));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx11, mc215));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx0, mc216));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx1, mc217));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx2, mc218));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx3, mc219));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx4, mc220));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx5, mc221));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx6, mc222));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx7, mc223));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx8, mc224));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx9, mc225));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx10, mc226));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx11, mc227));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx0, mc228));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx1, mc229));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx2, mc230));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx3, mc231));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx4, mc232));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx5, mc233));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx6, mc234));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx7, mc235));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx8, mc236));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx9, mc237));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx10, mc238));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx11, mc239));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx0, mc240));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx1, mc241));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx2, mc242));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx3, mc243));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx4, mc244));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx5, mc245));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx6, mc246));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx7, mc247));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx8, mc248));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx9, mc249));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx10, mc250));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx11, mc251));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx0, mc252));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx1, mc253));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx2, mc254));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx3, mc255));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx4, mc256));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx5, mc257));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx6, mc258));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx7, mc259));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx8, mc260));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx9, mc261));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx10, mc262));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx11, mc263));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx0, mc264));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx1, mc265));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx2, mc266));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx3, mc267));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx4, mc268));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx5, mc269));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx6, mc270));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx7, mc271));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx8, mc272));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx9, mc273));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx10, mc274));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx11, mc275));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx0, mc276));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx1, mc277));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx2, mc278));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx3, mc279));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx4, mc280));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx5, mc281));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx6, mc282));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx7, mc283));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx8, mc284));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx9, mc285));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx10, mc286));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx11, mc287))

#define setup_24(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0);\
                         msum16 = _mm_set_pd(0,0);\
                         msum17 = _mm_set_pd(0,0);\
                         msum18 = _mm_set_pd(0,0);\
                         msum19 = _mm_set_pd(0,0);\
                         msum20 = _mm_set_pd(0,0);\
                         msum21 = _mm_set_pd(0,0);\
                         msum22 = _mm_set_pd(0,0);\
                         msum23 = _mm_set_pd(0,0)

#define save_24() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  msum16 = _mm_hadd_pd(msum16, msum17);\
                  msum18 = _mm_hadd_pd(msum18, msum19);\
                  msum20 = _mm_hadd_pd(msum20, msum21);\
                  msum22 = _mm_hadd_pd(msum22, msum23);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14);\
                  _mm_storeu_pd(y + t1 + 16,msum16);\
                  _mm_storeu_pd(y + t1 + 18,msum18);\
                  _mm_storeu_pd(y + t1 + 20,msum20);\
                  _mm_storeu_pd(y + t1 + 22,msum22)

PetscErrorCode BSG_MatMult_24(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8, mx9, mx10, mx11, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, mc128, mc129, mc130, mc131, mc132, mc133, mc134, mc135, mc136, mc137, mc138, mc139, mc140, mc141, mc142, mc143, mc144, mc145, mc146, mc147, mc148, mc149, mc150, mc151, mc152, mc153, mc154, mc155, mc156, mc157, mc158, mc159, mc160, mc161, mc162, mc163, mc164, mc165, mc166, mc167, mc168, mc169, mc170, mc171, mc172, mc173, mc174, mc175, mc176, mc177, mc178, mc179, mc180, mc181, mc182, mc183, mc184, mc185, mc186, mc187, mc188, mc189, mc190, mc191, mc192, mc193, mc194, mc195, mc196, mc197, mc198, mc199, mc200, mc201, mc202, mc203, mc204, mc205, mc206, mc207, mc208, mc209, mc210, mc211, mc212, mc213, mc214, mc215, mc216, mc217, mc218, mc219, mc220, mc221, mc222, mc223, mc224, mc225, mc226, mc227, mc228, mc229, mc230, mc231, mc232, mc233, mc234, mc235, mc236, mc237, mc238, mc239, mc240, mc241, mc242, mc243, mc244, mc245, mc246, mc247, mc248, mc249, mc250, mc251, mc252, mc253, mc254, mc255, mc256, mc257, mc258, mc259, mc260, mc261, mc262, mc263, mc264, mc265, mc266, mc267, mc268, mc269, mc270, mc271, mc272, mc273, mc274, mc275, mc276, mc277, mc278, mc279, mc280, mc281, mc282, mc283, mc284, mc285, mc286, mc287, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15, msum16, msum17, msum18, msum19, msum20, msum21, msum22, msum23;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_24(0);
            inline_24(3);
            inline_24(4);
            inline_24(5);
            inline_24(6);
            save_24();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_24(1);
            inline_24(2);
            inline_24(3);
            inline_24(4);
            inline_24(5);
            inline_24(6);
            save_24();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_24(lda3);
            inline_24(1);
            inline_24(2);
            inline_24(3);
            inline_24(4);
            inline_24(5);
            inline_24(6);
            save_24();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_24(lda2);
            inline_24(0);
            inline_24(1);
            inline_24(2);
            inline_24(3);
            inline_24(4);
            inline_24(5);
            inline_24(6);
            save_24();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_24(lda1 - lda2);
            inline_24(0);
            inline_24(1);
            inline_24(2);
            inline_24(3);
            inline_24(4);
            inline_24(5);
            save_24();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_24(lda1 - lda3);
            inline_24(0);
            inline_24(1);
            inline_24(2);
            inline_24(3);
            inline_24(4);
            save_24();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_24(lda1 - 1);
            inline_24(0);
            inline_24(1);
            inline_24(2);
            inline_24(3);
            save_24();
        }
    }

PetscFunctionReturn(0);
}

#define inline_26(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mx8 = _mm_loadu_pd(xt##l + t1 + 16);\
                     mx9 = _mm_loadu_pd(xt##l + t1 + 18);\
                     mx10 = _mm_loadu_pd(xt##l + t1 + 20);\
                     mx11 = _mm_loadu_pd(xt##l + t1 + 22);\
                     mx12 = _mm_loadu_pd(xt##l + t1 + 24);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     mc128 = _mm_loadu_pd(ct##l + t2 + 256);\
                     mc129 = _mm_loadu_pd(ct##l + t2 + 258);\
                     mc130 = _mm_loadu_pd(ct##l + t2 + 260);\
                     mc131 = _mm_loadu_pd(ct##l + t2 + 262);\
                     mc132 = _mm_loadu_pd(ct##l + t2 + 264);\
                     mc133 = _mm_loadu_pd(ct##l + t2 + 266);\
                     mc134 = _mm_loadu_pd(ct##l + t2 + 268);\
                     mc135 = _mm_loadu_pd(ct##l + t2 + 270);\
                     mc136 = _mm_loadu_pd(ct##l + t2 + 272);\
                     mc137 = _mm_loadu_pd(ct##l + t2 + 274);\
                     mc138 = _mm_loadu_pd(ct##l + t2 + 276);\
                     mc139 = _mm_loadu_pd(ct##l + t2 + 278);\
                     mc140 = _mm_loadu_pd(ct##l + t2 + 280);\
                     mc141 = _mm_loadu_pd(ct##l + t2 + 282);\
                     mc142 = _mm_loadu_pd(ct##l + t2 + 284);\
                     mc143 = _mm_loadu_pd(ct##l + t2 + 286);\
                     mc144 = _mm_loadu_pd(ct##l + t2 + 288);\
                     mc145 = _mm_loadu_pd(ct##l + t2 + 290);\
                     mc146 = _mm_loadu_pd(ct##l + t2 + 292);\
                     mc147 = _mm_loadu_pd(ct##l + t2 + 294);\
                     mc148 = _mm_loadu_pd(ct##l + t2 + 296);\
                     mc149 = _mm_loadu_pd(ct##l + t2 + 298);\
                     mc150 = _mm_loadu_pd(ct##l + t2 + 300);\
                     mc151 = _mm_loadu_pd(ct##l + t2 + 302);\
                     mc152 = _mm_loadu_pd(ct##l + t2 + 304);\
                     mc153 = _mm_loadu_pd(ct##l + t2 + 306);\
                     mc154 = _mm_loadu_pd(ct##l + t2 + 308);\
                     mc155 = _mm_loadu_pd(ct##l + t2 + 310);\
                     mc156 = _mm_loadu_pd(ct##l + t2 + 312);\
                     mc157 = _mm_loadu_pd(ct##l + t2 + 314);\
                     mc158 = _mm_loadu_pd(ct##l + t2 + 316);\
                     mc159 = _mm_loadu_pd(ct##l + t2 + 318);\
                     mc160 = _mm_loadu_pd(ct##l + t2 + 320);\
                     mc161 = _mm_loadu_pd(ct##l + t2 + 322);\
                     mc162 = _mm_loadu_pd(ct##l + t2 + 324);\
                     mc163 = _mm_loadu_pd(ct##l + t2 + 326);\
                     mc164 = _mm_loadu_pd(ct##l + t2 + 328);\
                     mc165 = _mm_loadu_pd(ct##l + t2 + 330);\
                     mc166 = _mm_loadu_pd(ct##l + t2 + 332);\
                     mc167 = _mm_loadu_pd(ct##l + t2 + 334);\
                     mc168 = _mm_loadu_pd(ct##l + t2 + 336);\
                     mc169 = _mm_loadu_pd(ct##l + t2 + 338);\
                     mc170 = _mm_loadu_pd(ct##l + t2 + 340);\
                     mc171 = _mm_loadu_pd(ct##l + t2 + 342);\
                     mc172 = _mm_loadu_pd(ct##l + t2 + 344);\
                     mc173 = _mm_loadu_pd(ct##l + t2 + 346);\
                     mc174 = _mm_loadu_pd(ct##l + t2 + 348);\
                     mc175 = _mm_loadu_pd(ct##l + t2 + 350);\
                     mc176 = _mm_loadu_pd(ct##l + t2 + 352);\
                     mc177 = _mm_loadu_pd(ct##l + t2 + 354);\
                     mc178 = _mm_loadu_pd(ct##l + t2 + 356);\
                     mc179 = _mm_loadu_pd(ct##l + t2 + 358);\
                     mc180 = _mm_loadu_pd(ct##l + t2 + 360);\
                     mc181 = _mm_loadu_pd(ct##l + t2 + 362);\
                     mc182 = _mm_loadu_pd(ct##l + t2 + 364);\
                     mc183 = _mm_loadu_pd(ct##l + t2 + 366);\
                     mc184 = _mm_loadu_pd(ct##l + t2 + 368);\
                     mc185 = _mm_loadu_pd(ct##l + t2 + 370);\
                     mc186 = _mm_loadu_pd(ct##l + t2 + 372);\
                     mc187 = _mm_loadu_pd(ct##l + t2 + 374);\
                     mc188 = _mm_loadu_pd(ct##l + t2 + 376);\
                     mc189 = _mm_loadu_pd(ct##l + t2 + 378);\
                     mc190 = _mm_loadu_pd(ct##l + t2 + 380);\
                     mc191 = _mm_loadu_pd(ct##l + t2 + 382);\
                     mc192 = _mm_loadu_pd(ct##l + t2 + 384);\
                     mc193 = _mm_loadu_pd(ct##l + t2 + 386);\
                     mc194 = _mm_loadu_pd(ct##l + t2 + 388);\
                     mc195 = _mm_loadu_pd(ct##l + t2 + 390);\
                     mc196 = _mm_loadu_pd(ct##l + t2 + 392);\
                     mc197 = _mm_loadu_pd(ct##l + t2 + 394);\
                     mc198 = _mm_loadu_pd(ct##l + t2 + 396);\
                     mc199 = _mm_loadu_pd(ct##l + t2 + 398);\
                     mc200 = _mm_loadu_pd(ct##l + t2 + 400);\
                     mc201 = _mm_loadu_pd(ct##l + t2 + 402);\
                     mc202 = _mm_loadu_pd(ct##l + t2 + 404);\
                     mc203 = _mm_loadu_pd(ct##l + t2 + 406);\
                     mc204 = _mm_loadu_pd(ct##l + t2 + 408);\
                     mc205 = _mm_loadu_pd(ct##l + t2 + 410);\
                     mc206 = _mm_loadu_pd(ct##l + t2 + 412);\
                     mc207 = _mm_loadu_pd(ct##l + t2 + 414);\
                     mc208 = _mm_loadu_pd(ct##l + t2 + 416);\
                     mc209 = _mm_loadu_pd(ct##l + t2 + 418);\
                     mc210 = _mm_loadu_pd(ct##l + t2 + 420);\
                     mc211 = _mm_loadu_pd(ct##l + t2 + 422);\
                     mc212 = _mm_loadu_pd(ct##l + t2 + 424);\
                     mc213 = _mm_loadu_pd(ct##l + t2 + 426);\
                     mc214 = _mm_loadu_pd(ct##l + t2 + 428);\
                     mc215 = _mm_loadu_pd(ct##l + t2 + 430);\
                     mc216 = _mm_loadu_pd(ct##l + t2 + 432);\
                     mc217 = _mm_loadu_pd(ct##l + t2 + 434);\
                     mc218 = _mm_loadu_pd(ct##l + t2 + 436);\
                     mc219 = _mm_loadu_pd(ct##l + t2 + 438);\
                     mc220 = _mm_loadu_pd(ct##l + t2 + 440);\
                     mc221 = _mm_loadu_pd(ct##l + t2 + 442);\
                     mc222 = _mm_loadu_pd(ct##l + t2 + 444);\
                     mc223 = _mm_loadu_pd(ct##l + t2 + 446);\
                     mc224 = _mm_loadu_pd(ct##l + t2 + 448);\
                     mc225 = _mm_loadu_pd(ct##l + t2 + 450);\
                     mc226 = _mm_loadu_pd(ct##l + t2 + 452);\
                     mc227 = _mm_loadu_pd(ct##l + t2 + 454);\
                     mc228 = _mm_loadu_pd(ct##l + t2 + 456);\
                     mc229 = _mm_loadu_pd(ct##l + t2 + 458);\
                     mc230 = _mm_loadu_pd(ct##l + t2 + 460);\
                     mc231 = _mm_loadu_pd(ct##l + t2 + 462);\
                     mc232 = _mm_loadu_pd(ct##l + t2 + 464);\
                     mc233 = _mm_loadu_pd(ct##l + t2 + 466);\
                     mc234 = _mm_loadu_pd(ct##l + t2 + 468);\
                     mc235 = _mm_loadu_pd(ct##l + t2 + 470);\
                     mc236 = _mm_loadu_pd(ct##l + t2 + 472);\
                     mc237 = _mm_loadu_pd(ct##l + t2 + 474);\
                     mc238 = _mm_loadu_pd(ct##l + t2 + 476);\
                     mc239 = _mm_loadu_pd(ct##l + t2 + 478);\
                     mc240 = _mm_loadu_pd(ct##l + t2 + 480);\
                     mc241 = _mm_loadu_pd(ct##l + t2 + 482);\
                     mc242 = _mm_loadu_pd(ct##l + t2 + 484);\
                     mc243 = _mm_loadu_pd(ct##l + t2 + 486);\
                     mc244 = _mm_loadu_pd(ct##l + t2 + 488);\
                     mc245 = _mm_loadu_pd(ct##l + t2 + 490);\
                     mc246 = _mm_loadu_pd(ct##l + t2 + 492);\
                     mc247 = _mm_loadu_pd(ct##l + t2 + 494);\
                     mc248 = _mm_loadu_pd(ct##l + t2 + 496);\
                     mc249 = _mm_loadu_pd(ct##l + t2 + 498);\
                     mc250 = _mm_loadu_pd(ct##l + t2 + 500);\
                     mc251 = _mm_loadu_pd(ct##l + t2 + 502);\
                     mc252 = _mm_loadu_pd(ct##l + t2 + 504);\
                     mc253 = _mm_loadu_pd(ct##l + t2 + 506);\
                     mc254 = _mm_loadu_pd(ct##l + t2 + 508);\
                     mc255 = _mm_loadu_pd(ct##l + t2 + 510);\
                     mc256 = _mm_loadu_pd(ct##l + t2 + 512);\
                     mc257 = _mm_loadu_pd(ct##l + t2 + 514);\
                     mc258 = _mm_loadu_pd(ct##l + t2 + 516);\
                     mc259 = _mm_loadu_pd(ct##l + t2 + 518);\
                     mc260 = _mm_loadu_pd(ct##l + t2 + 520);\
                     mc261 = _mm_loadu_pd(ct##l + t2 + 522);\
                     mc262 = _mm_loadu_pd(ct##l + t2 + 524);\
                     mc263 = _mm_loadu_pd(ct##l + t2 + 526);\
                     mc264 = _mm_loadu_pd(ct##l + t2 + 528);\
                     mc265 = _mm_loadu_pd(ct##l + t2 + 530);\
                     mc266 = _mm_loadu_pd(ct##l + t2 + 532);\
                     mc267 = _mm_loadu_pd(ct##l + t2 + 534);\
                     mc268 = _mm_loadu_pd(ct##l + t2 + 536);\
                     mc269 = _mm_loadu_pd(ct##l + t2 + 538);\
                     mc270 = _mm_loadu_pd(ct##l + t2 + 540);\
                     mc271 = _mm_loadu_pd(ct##l + t2 + 542);\
                     mc272 = _mm_loadu_pd(ct##l + t2 + 544);\
                     mc273 = _mm_loadu_pd(ct##l + t2 + 546);\
                     mc274 = _mm_loadu_pd(ct##l + t2 + 548);\
                     mc275 = _mm_loadu_pd(ct##l + t2 + 550);\
                     mc276 = _mm_loadu_pd(ct##l + t2 + 552);\
                     mc277 = _mm_loadu_pd(ct##l + t2 + 554);\
                     mc278 = _mm_loadu_pd(ct##l + t2 + 556);\
                     mc279 = _mm_loadu_pd(ct##l + t2 + 558);\
                     mc280 = _mm_loadu_pd(ct##l + t2 + 560);\
                     mc281 = _mm_loadu_pd(ct##l + t2 + 562);\
                     mc282 = _mm_loadu_pd(ct##l + t2 + 564);\
                     mc283 = _mm_loadu_pd(ct##l + t2 + 566);\
                     mc284 = _mm_loadu_pd(ct##l + t2 + 568);\
                     mc285 = _mm_loadu_pd(ct##l + t2 + 570);\
                     mc286 = _mm_loadu_pd(ct##l + t2 + 572);\
                     mc287 = _mm_loadu_pd(ct##l + t2 + 574);\
                     mc288 = _mm_loadu_pd(ct##l + t2 + 576);\
                     mc289 = _mm_loadu_pd(ct##l + t2 + 578);\
                     mc290 = _mm_loadu_pd(ct##l + t2 + 580);\
                     mc291 = _mm_loadu_pd(ct##l + t2 + 582);\
                     mc292 = _mm_loadu_pd(ct##l + t2 + 584);\
                     mc293 = _mm_loadu_pd(ct##l + t2 + 586);\
                     mc294 = _mm_loadu_pd(ct##l + t2 + 588);\
                     mc295 = _mm_loadu_pd(ct##l + t2 + 590);\
                     mc296 = _mm_loadu_pd(ct##l + t2 + 592);\
                     mc297 = _mm_loadu_pd(ct##l + t2 + 594);\
                     mc298 = _mm_loadu_pd(ct##l + t2 + 596);\
                     mc299 = _mm_loadu_pd(ct##l + t2 + 598);\
                     mc300 = _mm_loadu_pd(ct##l + t2 + 600);\
                     mc301 = _mm_loadu_pd(ct##l + t2 + 602);\
                     mc302 = _mm_loadu_pd(ct##l + t2 + 604);\
                     mc303 = _mm_loadu_pd(ct##l + t2 + 606);\
                     mc304 = _mm_loadu_pd(ct##l + t2 + 608);\
                     mc305 = _mm_loadu_pd(ct##l + t2 + 610);\
                     mc306 = _mm_loadu_pd(ct##l + t2 + 612);\
                     mc307 = _mm_loadu_pd(ct##l + t2 + 614);\
                     mc308 = _mm_loadu_pd(ct##l + t2 + 616);\
                     mc309 = _mm_loadu_pd(ct##l + t2 + 618);\
                     mc310 = _mm_loadu_pd(ct##l + t2 + 620);\
                     mc311 = _mm_loadu_pd(ct##l + t2 + 622);\
                     mc312 = _mm_loadu_pd(ct##l + t2 + 624);\
                     mc313 = _mm_loadu_pd(ct##l + t2 + 626);\
                     mc314 = _mm_loadu_pd(ct##l + t2 + 628);\
                     mc315 = _mm_loadu_pd(ct##l + t2 + 630);\
                     mc316 = _mm_loadu_pd(ct##l + t2 + 632);\
                     mc317 = _mm_loadu_pd(ct##l + t2 + 634);\
                     mc318 = _mm_loadu_pd(ct##l + t2 + 636);\
                     mc319 = _mm_loadu_pd(ct##l + t2 + 638);\
                     mc320 = _mm_loadu_pd(ct##l + t2 + 640);\
                     mc321 = _mm_loadu_pd(ct##l + t2 + 642);\
                     mc322 = _mm_loadu_pd(ct##l + t2 + 644);\
                     mc323 = _mm_loadu_pd(ct##l + t2 + 646);\
                     mc324 = _mm_loadu_pd(ct##l + t2 + 648);\
                     mc325 = _mm_loadu_pd(ct##l + t2 + 650);\
                     mc326 = _mm_loadu_pd(ct##l + t2 + 652);\
                     mc327 = _mm_loadu_pd(ct##l + t2 + 654);\
                     mc328 = _mm_loadu_pd(ct##l + t2 + 656);\
                     mc329 = _mm_loadu_pd(ct##l + t2 + 658);\
                     mc330 = _mm_loadu_pd(ct##l + t2 + 660);\
                     mc331 = _mm_loadu_pd(ct##l + t2 + 662);\
                     mc332 = _mm_loadu_pd(ct##l + t2 + 664);\
                     mc333 = _mm_loadu_pd(ct##l + t2 + 666);\
                     mc334 = _mm_loadu_pd(ct##l + t2 + 668);\
                     mc335 = _mm_loadu_pd(ct##l + t2 + 670);\
                     mc336 = _mm_loadu_pd(ct##l + t2 + 672);\
                     mc337 = _mm_loadu_pd(ct##l + t2 + 674);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx8, mc8));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx9, mc9));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx10, mc10));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx11, mc11));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx12, mc12));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc13));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc15));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc16));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc17));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc18));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc19));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc20));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx8, mc21));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx9, mc22));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx10, mc23));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx11, mc24));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx12, mc25));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc26));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc27));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc28));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc29));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc30));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc31));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc32));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc33));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx8, mc34));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx9, mc35));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx10, mc36));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx11, mc37));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx12, mc38));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc39));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc40));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc41));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc42));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc43));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc44));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc45));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc46));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx8, mc47));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx9, mc48));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx10, mc49));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx11, mc50));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx12, mc51));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc52));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc53));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc54));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc55));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc56));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc57));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc58));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc59));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx8, mc60));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx9, mc61));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx10, mc62));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx11, mc63));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx12, mc64));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc65));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc66));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc67));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc68));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc69));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc70));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc71));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc72));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx8, mc73));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx9, mc74));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx10, mc75));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx11, mc76));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx12, mc77));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc78));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc79));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc80));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc81));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc82));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc83));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc84));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc85));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx8, mc86));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx9, mc87));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx10, mc88));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx11, mc89));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx12, mc90));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc91));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc92));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc93));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc94));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc95));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc96));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc97));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc98));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx8, mc99));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx9, mc100));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx10, mc101));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx11, mc102));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx12, mc103));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc104));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc105));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc106));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc107));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc108));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc109));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc110));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc111));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx8, mc112));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx9, mc113));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx10, mc114));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx11, mc115));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx12, mc116));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc117));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc118));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc119));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc120));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc121));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc122));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc123));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc124));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx8, mc125));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx9, mc126));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx10, mc127));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx11, mc128));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx12, mc129));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc130));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc131));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc132));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc133));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc134));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc135));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc136));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc137));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx8, mc138));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx9, mc139));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx10, mc140));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx11, mc141));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx12, mc142));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc143));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc144));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc145));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc146));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc147));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc148));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc149));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc150));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx8, mc151));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx9, mc152));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx10, mc153));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx11, mc154));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx12, mc155));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc156));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc157));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc158));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc159));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc160));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc161));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc162));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc163));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx8, mc164));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx9, mc165));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx10, mc166));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx11, mc167));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx12, mc168));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc169));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc170));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc171));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc172));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc173));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc174));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc175));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc176));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx8, mc177));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx9, mc178));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx10, mc179));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx11, mc180));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx12, mc181));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc182));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc183));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc184));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc185));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc186));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc187));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc188));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc189));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx8, mc190));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx9, mc191));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx10, mc192));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx11, mc193));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx12, mc194));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc195));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc196));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc197));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc198));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc199));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc200));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc201));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc202));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx8, mc203));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx9, mc204));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx10, mc205));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx11, mc206));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx12, mc207));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx0, mc208));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx1, mc209));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx2, mc210));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx3, mc211));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx4, mc212));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx5, mc213));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx6, mc214));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx7, mc215));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx8, mc216));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx9, mc217));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx10, mc218));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx11, mc219));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx12, mc220));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx0, mc221));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx1, mc222));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx2, mc223));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx3, mc224));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx4, mc225));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx5, mc226));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx6, mc227));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx7, mc228));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx8, mc229));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx9, mc230));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx10, mc231));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx11, mc232));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx12, mc233));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx0, mc234));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx1, mc235));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx2, mc236));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx3, mc237));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx4, mc238));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx5, mc239));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx6, mc240));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx7, mc241));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx8, mc242));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx9, mc243));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx10, mc244));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx11, mc245));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx12, mc246));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx0, mc247));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx1, mc248));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx2, mc249));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx3, mc250));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx4, mc251));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx5, mc252));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx6, mc253));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx7, mc254));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx8, mc255));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx9, mc256));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx10, mc257));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx11, mc258));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx12, mc259));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx0, mc260));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx1, mc261));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx2, mc262));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx3, mc263));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx4, mc264));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx5, mc265));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx6, mc266));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx7, mc267));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx8, mc268));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx9, mc269));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx10, mc270));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx11, mc271));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx12, mc272));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx0, mc273));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx1, mc274));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx2, mc275));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx3, mc276));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx4, mc277));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx5, mc278));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx6, mc279));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx7, mc280));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx8, mc281));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx9, mc282));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx10, mc283));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx11, mc284));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx12, mc285));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx0, mc286));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx1, mc287));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx2, mc288));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx3, mc289));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx4, mc290));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx5, mc291));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx6, mc292));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx7, mc293));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx8, mc294));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx9, mc295));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx10, mc296));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx11, mc297));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx12, mc298));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx0, mc299));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx1, mc300));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx2, mc301));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx3, mc302));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx4, mc303));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx5, mc304));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx6, mc305));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx7, mc306));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx8, mc307));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx9, mc308));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx10, mc309));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx11, mc310));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx12, mc311));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx0, mc312));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx1, mc313));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx2, mc314));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx3, mc315));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx4, mc316));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx5, mc317));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx6, mc318));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx7, mc319));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx8, mc320));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx9, mc321));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx10, mc322));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx11, mc323));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx12, mc324));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx0, mc325));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx1, mc326));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx2, mc327));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx3, mc328));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx4, mc329));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx5, mc330));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx6, mc331));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx7, mc332));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx8, mc333));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx9, mc334));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx10, mc335));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx11, mc336));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx12, mc337))

#define setup_26(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0);\
                         msum16 = _mm_set_pd(0,0);\
                         msum17 = _mm_set_pd(0,0);\
                         msum18 = _mm_set_pd(0,0);\
                         msum19 = _mm_set_pd(0,0);\
                         msum20 = _mm_set_pd(0,0);\
                         msum21 = _mm_set_pd(0,0);\
                         msum22 = _mm_set_pd(0,0);\
                         msum23 = _mm_set_pd(0,0);\
                         msum24 = _mm_set_pd(0,0);\
                         msum25 = _mm_set_pd(0,0)

#define save_26() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  msum16 = _mm_hadd_pd(msum16, msum17);\
                  msum18 = _mm_hadd_pd(msum18, msum19);\
                  msum20 = _mm_hadd_pd(msum20, msum21);\
                  msum22 = _mm_hadd_pd(msum22, msum23);\
                  msum24 = _mm_hadd_pd(msum24, msum25);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14);\
                  _mm_storeu_pd(y + t1 + 16,msum16);\
                  _mm_storeu_pd(y + t1 + 18,msum18);\
                  _mm_storeu_pd(y + t1 + 20,msum20);\
                  _mm_storeu_pd(y + t1 + 22,msum22);\
                  _mm_storeu_pd(y + t1 + 24,msum24)

PetscErrorCode BSG_MatMult_26(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8, mx9, mx10, mx11, mx12, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, mc128, mc129, mc130, mc131, mc132, mc133, mc134, mc135, mc136, mc137, mc138, mc139, mc140, mc141, mc142, mc143, mc144, mc145, mc146, mc147, mc148, mc149, mc150, mc151, mc152, mc153, mc154, mc155, mc156, mc157, mc158, mc159, mc160, mc161, mc162, mc163, mc164, mc165, mc166, mc167, mc168, mc169, mc170, mc171, mc172, mc173, mc174, mc175, mc176, mc177, mc178, mc179, mc180, mc181, mc182, mc183, mc184, mc185, mc186, mc187, mc188, mc189, mc190, mc191, mc192, mc193, mc194, mc195, mc196, mc197, mc198, mc199, mc200, mc201, mc202, mc203, mc204, mc205, mc206, mc207, mc208, mc209, mc210, mc211, mc212, mc213, mc214, mc215, mc216, mc217, mc218, mc219, mc220, mc221, mc222, mc223, mc224, mc225, mc226, mc227, mc228, mc229, mc230, mc231, mc232, mc233, mc234, mc235, mc236, mc237, mc238, mc239, mc240, mc241, mc242, mc243, mc244, mc245, mc246, mc247, mc248, mc249, mc250, mc251, mc252, mc253, mc254, mc255, mc256, mc257, mc258, mc259, mc260, mc261, mc262, mc263, mc264, mc265, mc266, mc267, mc268, mc269, mc270, mc271, mc272, mc273, mc274, mc275, mc276, mc277, mc278, mc279, mc280, mc281, mc282, mc283, mc284, mc285, mc286, mc287, mc288, mc289, mc290, mc291, mc292, mc293, mc294, mc295, mc296, mc297, mc298, mc299, mc300, mc301, mc302, mc303, mc304, mc305, mc306, mc307, mc308, mc309, mc310, mc311, mc312, mc313, mc314, mc315, mc316, mc317, mc318, mc319, mc320, mc321, mc322, mc323, mc324, mc325, mc326, mc327, mc328, mc329, mc330, mc331, mc332, mc333, mc334, mc335, mc336, mc337, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15, msum16, msum17, msum18, msum19, msum20, msum21, msum22, msum23, msum24, msum25;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_26(0);
            inline_26(3);
            inline_26(4);
            inline_26(5);
            inline_26(6);
            save_26();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_26(1);
            inline_26(2);
            inline_26(3);
            inline_26(4);
            inline_26(5);
            inline_26(6);
            save_26();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_26(lda3);
            inline_26(1);
            inline_26(2);
            inline_26(3);
            inline_26(4);
            inline_26(5);
            inline_26(6);
            save_26();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_26(lda2);
            inline_26(0);
            inline_26(1);
            inline_26(2);
            inline_26(3);
            inline_26(4);
            inline_26(5);
            inline_26(6);
            save_26();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_26(lda1 - lda2);
            inline_26(0);
            inline_26(1);
            inline_26(2);
            inline_26(3);
            inline_26(4);
            inline_26(5);
            save_26();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_26(lda1 - lda3);
            inline_26(0);
            inline_26(1);
            inline_26(2);
            inline_26(3);
            inline_26(4);
            save_26();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_26(lda1 - 1);
            inline_26(0);
            inline_26(1);
            inline_26(2);
            inline_26(3);
            save_26();
        }
    }

PetscFunctionReturn(0);
}

#define inline_28(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mx8 = _mm_loadu_pd(xt##l + t1 + 16);\
                     mx9 = _mm_loadu_pd(xt##l + t1 + 18);\
                     mx10 = _mm_loadu_pd(xt##l + t1 + 20);\
                     mx11 = _mm_loadu_pd(xt##l + t1 + 22);\
                     mx12 = _mm_loadu_pd(xt##l + t1 + 24);\
                     mx13 = _mm_loadu_pd(xt##l + t1 + 26);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     mc128 = _mm_loadu_pd(ct##l + t2 + 256);\
                     mc129 = _mm_loadu_pd(ct##l + t2 + 258);\
                     mc130 = _mm_loadu_pd(ct##l + t2 + 260);\
                     mc131 = _mm_loadu_pd(ct##l + t2 + 262);\
                     mc132 = _mm_loadu_pd(ct##l + t2 + 264);\
                     mc133 = _mm_loadu_pd(ct##l + t2 + 266);\
                     mc134 = _mm_loadu_pd(ct##l + t2 + 268);\
                     mc135 = _mm_loadu_pd(ct##l + t2 + 270);\
                     mc136 = _mm_loadu_pd(ct##l + t2 + 272);\
                     mc137 = _mm_loadu_pd(ct##l + t2 + 274);\
                     mc138 = _mm_loadu_pd(ct##l + t2 + 276);\
                     mc139 = _mm_loadu_pd(ct##l + t2 + 278);\
                     mc140 = _mm_loadu_pd(ct##l + t2 + 280);\
                     mc141 = _mm_loadu_pd(ct##l + t2 + 282);\
                     mc142 = _mm_loadu_pd(ct##l + t2 + 284);\
                     mc143 = _mm_loadu_pd(ct##l + t2 + 286);\
                     mc144 = _mm_loadu_pd(ct##l + t2 + 288);\
                     mc145 = _mm_loadu_pd(ct##l + t2 + 290);\
                     mc146 = _mm_loadu_pd(ct##l + t2 + 292);\
                     mc147 = _mm_loadu_pd(ct##l + t2 + 294);\
                     mc148 = _mm_loadu_pd(ct##l + t2 + 296);\
                     mc149 = _mm_loadu_pd(ct##l + t2 + 298);\
                     mc150 = _mm_loadu_pd(ct##l + t2 + 300);\
                     mc151 = _mm_loadu_pd(ct##l + t2 + 302);\
                     mc152 = _mm_loadu_pd(ct##l + t2 + 304);\
                     mc153 = _mm_loadu_pd(ct##l + t2 + 306);\
                     mc154 = _mm_loadu_pd(ct##l + t2 + 308);\
                     mc155 = _mm_loadu_pd(ct##l + t2 + 310);\
                     mc156 = _mm_loadu_pd(ct##l + t2 + 312);\
                     mc157 = _mm_loadu_pd(ct##l + t2 + 314);\
                     mc158 = _mm_loadu_pd(ct##l + t2 + 316);\
                     mc159 = _mm_loadu_pd(ct##l + t2 + 318);\
                     mc160 = _mm_loadu_pd(ct##l + t2 + 320);\
                     mc161 = _mm_loadu_pd(ct##l + t2 + 322);\
                     mc162 = _mm_loadu_pd(ct##l + t2 + 324);\
                     mc163 = _mm_loadu_pd(ct##l + t2 + 326);\
                     mc164 = _mm_loadu_pd(ct##l + t2 + 328);\
                     mc165 = _mm_loadu_pd(ct##l + t2 + 330);\
                     mc166 = _mm_loadu_pd(ct##l + t2 + 332);\
                     mc167 = _mm_loadu_pd(ct##l + t2 + 334);\
                     mc168 = _mm_loadu_pd(ct##l + t2 + 336);\
                     mc169 = _mm_loadu_pd(ct##l + t2 + 338);\
                     mc170 = _mm_loadu_pd(ct##l + t2 + 340);\
                     mc171 = _mm_loadu_pd(ct##l + t2 + 342);\
                     mc172 = _mm_loadu_pd(ct##l + t2 + 344);\
                     mc173 = _mm_loadu_pd(ct##l + t2 + 346);\
                     mc174 = _mm_loadu_pd(ct##l + t2 + 348);\
                     mc175 = _mm_loadu_pd(ct##l + t2 + 350);\
                     mc176 = _mm_loadu_pd(ct##l + t2 + 352);\
                     mc177 = _mm_loadu_pd(ct##l + t2 + 354);\
                     mc178 = _mm_loadu_pd(ct##l + t2 + 356);\
                     mc179 = _mm_loadu_pd(ct##l + t2 + 358);\
                     mc180 = _mm_loadu_pd(ct##l + t2 + 360);\
                     mc181 = _mm_loadu_pd(ct##l + t2 + 362);\
                     mc182 = _mm_loadu_pd(ct##l + t2 + 364);\
                     mc183 = _mm_loadu_pd(ct##l + t2 + 366);\
                     mc184 = _mm_loadu_pd(ct##l + t2 + 368);\
                     mc185 = _mm_loadu_pd(ct##l + t2 + 370);\
                     mc186 = _mm_loadu_pd(ct##l + t2 + 372);\
                     mc187 = _mm_loadu_pd(ct##l + t2 + 374);\
                     mc188 = _mm_loadu_pd(ct##l + t2 + 376);\
                     mc189 = _mm_loadu_pd(ct##l + t2 + 378);\
                     mc190 = _mm_loadu_pd(ct##l + t2 + 380);\
                     mc191 = _mm_loadu_pd(ct##l + t2 + 382);\
                     mc192 = _mm_loadu_pd(ct##l + t2 + 384);\
                     mc193 = _mm_loadu_pd(ct##l + t2 + 386);\
                     mc194 = _mm_loadu_pd(ct##l + t2 + 388);\
                     mc195 = _mm_loadu_pd(ct##l + t2 + 390);\
                     mc196 = _mm_loadu_pd(ct##l + t2 + 392);\
                     mc197 = _mm_loadu_pd(ct##l + t2 + 394);\
                     mc198 = _mm_loadu_pd(ct##l + t2 + 396);\
                     mc199 = _mm_loadu_pd(ct##l + t2 + 398);\
                     mc200 = _mm_loadu_pd(ct##l + t2 + 400);\
                     mc201 = _mm_loadu_pd(ct##l + t2 + 402);\
                     mc202 = _mm_loadu_pd(ct##l + t2 + 404);\
                     mc203 = _mm_loadu_pd(ct##l + t2 + 406);\
                     mc204 = _mm_loadu_pd(ct##l + t2 + 408);\
                     mc205 = _mm_loadu_pd(ct##l + t2 + 410);\
                     mc206 = _mm_loadu_pd(ct##l + t2 + 412);\
                     mc207 = _mm_loadu_pd(ct##l + t2 + 414);\
                     mc208 = _mm_loadu_pd(ct##l + t2 + 416);\
                     mc209 = _mm_loadu_pd(ct##l + t2 + 418);\
                     mc210 = _mm_loadu_pd(ct##l + t2 + 420);\
                     mc211 = _mm_loadu_pd(ct##l + t2 + 422);\
                     mc212 = _mm_loadu_pd(ct##l + t2 + 424);\
                     mc213 = _mm_loadu_pd(ct##l + t2 + 426);\
                     mc214 = _mm_loadu_pd(ct##l + t2 + 428);\
                     mc215 = _mm_loadu_pd(ct##l + t2 + 430);\
                     mc216 = _mm_loadu_pd(ct##l + t2 + 432);\
                     mc217 = _mm_loadu_pd(ct##l + t2 + 434);\
                     mc218 = _mm_loadu_pd(ct##l + t2 + 436);\
                     mc219 = _mm_loadu_pd(ct##l + t2 + 438);\
                     mc220 = _mm_loadu_pd(ct##l + t2 + 440);\
                     mc221 = _mm_loadu_pd(ct##l + t2 + 442);\
                     mc222 = _mm_loadu_pd(ct##l + t2 + 444);\
                     mc223 = _mm_loadu_pd(ct##l + t2 + 446);\
                     mc224 = _mm_loadu_pd(ct##l + t2 + 448);\
                     mc225 = _mm_loadu_pd(ct##l + t2 + 450);\
                     mc226 = _mm_loadu_pd(ct##l + t2 + 452);\
                     mc227 = _mm_loadu_pd(ct##l + t2 + 454);\
                     mc228 = _mm_loadu_pd(ct##l + t2 + 456);\
                     mc229 = _mm_loadu_pd(ct##l + t2 + 458);\
                     mc230 = _mm_loadu_pd(ct##l + t2 + 460);\
                     mc231 = _mm_loadu_pd(ct##l + t2 + 462);\
                     mc232 = _mm_loadu_pd(ct##l + t2 + 464);\
                     mc233 = _mm_loadu_pd(ct##l + t2 + 466);\
                     mc234 = _mm_loadu_pd(ct##l + t2 + 468);\
                     mc235 = _mm_loadu_pd(ct##l + t2 + 470);\
                     mc236 = _mm_loadu_pd(ct##l + t2 + 472);\
                     mc237 = _mm_loadu_pd(ct##l + t2 + 474);\
                     mc238 = _mm_loadu_pd(ct##l + t2 + 476);\
                     mc239 = _mm_loadu_pd(ct##l + t2 + 478);\
                     mc240 = _mm_loadu_pd(ct##l + t2 + 480);\
                     mc241 = _mm_loadu_pd(ct##l + t2 + 482);\
                     mc242 = _mm_loadu_pd(ct##l + t2 + 484);\
                     mc243 = _mm_loadu_pd(ct##l + t2 + 486);\
                     mc244 = _mm_loadu_pd(ct##l + t2 + 488);\
                     mc245 = _mm_loadu_pd(ct##l + t2 + 490);\
                     mc246 = _mm_loadu_pd(ct##l + t2 + 492);\
                     mc247 = _mm_loadu_pd(ct##l + t2 + 494);\
                     mc248 = _mm_loadu_pd(ct##l + t2 + 496);\
                     mc249 = _mm_loadu_pd(ct##l + t2 + 498);\
                     mc250 = _mm_loadu_pd(ct##l + t2 + 500);\
                     mc251 = _mm_loadu_pd(ct##l + t2 + 502);\
                     mc252 = _mm_loadu_pd(ct##l + t2 + 504);\
                     mc253 = _mm_loadu_pd(ct##l + t2 + 506);\
                     mc254 = _mm_loadu_pd(ct##l + t2 + 508);\
                     mc255 = _mm_loadu_pd(ct##l + t2 + 510);\
                     mc256 = _mm_loadu_pd(ct##l + t2 + 512);\
                     mc257 = _mm_loadu_pd(ct##l + t2 + 514);\
                     mc258 = _mm_loadu_pd(ct##l + t2 + 516);\
                     mc259 = _mm_loadu_pd(ct##l + t2 + 518);\
                     mc260 = _mm_loadu_pd(ct##l + t2 + 520);\
                     mc261 = _mm_loadu_pd(ct##l + t2 + 522);\
                     mc262 = _mm_loadu_pd(ct##l + t2 + 524);\
                     mc263 = _mm_loadu_pd(ct##l + t2 + 526);\
                     mc264 = _mm_loadu_pd(ct##l + t2 + 528);\
                     mc265 = _mm_loadu_pd(ct##l + t2 + 530);\
                     mc266 = _mm_loadu_pd(ct##l + t2 + 532);\
                     mc267 = _mm_loadu_pd(ct##l + t2 + 534);\
                     mc268 = _mm_loadu_pd(ct##l + t2 + 536);\
                     mc269 = _mm_loadu_pd(ct##l + t2 + 538);\
                     mc270 = _mm_loadu_pd(ct##l + t2 + 540);\
                     mc271 = _mm_loadu_pd(ct##l + t2 + 542);\
                     mc272 = _mm_loadu_pd(ct##l + t2 + 544);\
                     mc273 = _mm_loadu_pd(ct##l + t2 + 546);\
                     mc274 = _mm_loadu_pd(ct##l + t2 + 548);\
                     mc275 = _mm_loadu_pd(ct##l + t2 + 550);\
                     mc276 = _mm_loadu_pd(ct##l + t2 + 552);\
                     mc277 = _mm_loadu_pd(ct##l + t2 + 554);\
                     mc278 = _mm_loadu_pd(ct##l + t2 + 556);\
                     mc279 = _mm_loadu_pd(ct##l + t2 + 558);\
                     mc280 = _mm_loadu_pd(ct##l + t2 + 560);\
                     mc281 = _mm_loadu_pd(ct##l + t2 + 562);\
                     mc282 = _mm_loadu_pd(ct##l + t2 + 564);\
                     mc283 = _mm_loadu_pd(ct##l + t2 + 566);\
                     mc284 = _mm_loadu_pd(ct##l + t2 + 568);\
                     mc285 = _mm_loadu_pd(ct##l + t2 + 570);\
                     mc286 = _mm_loadu_pd(ct##l + t2 + 572);\
                     mc287 = _mm_loadu_pd(ct##l + t2 + 574);\
                     mc288 = _mm_loadu_pd(ct##l + t2 + 576);\
                     mc289 = _mm_loadu_pd(ct##l + t2 + 578);\
                     mc290 = _mm_loadu_pd(ct##l + t2 + 580);\
                     mc291 = _mm_loadu_pd(ct##l + t2 + 582);\
                     mc292 = _mm_loadu_pd(ct##l + t2 + 584);\
                     mc293 = _mm_loadu_pd(ct##l + t2 + 586);\
                     mc294 = _mm_loadu_pd(ct##l + t2 + 588);\
                     mc295 = _mm_loadu_pd(ct##l + t2 + 590);\
                     mc296 = _mm_loadu_pd(ct##l + t2 + 592);\
                     mc297 = _mm_loadu_pd(ct##l + t2 + 594);\
                     mc298 = _mm_loadu_pd(ct##l + t2 + 596);\
                     mc299 = _mm_loadu_pd(ct##l + t2 + 598);\
                     mc300 = _mm_loadu_pd(ct##l + t2 + 600);\
                     mc301 = _mm_loadu_pd(ct##l + t2 + 602);\
                     mc302 = _mm_loadu_pd(ct##l + t2 + 604);\
                     mc303 = _mm_loadu_pd(ct##l + t2 + 606);\
                     mc304 = _mm_loadu_pd(ct##l + t2 + 608);\
                     mc305 = _mm_loadu_pd(ct##l + t2 + 610);\
                     mc306 = _mm_loadu_pd(ct##l + t2 + 612);\
                     mc307 = _mm_loadu_pd(ct##l + t2 + 614);\
                     mc308 = _mm_loadu_pd(ct##l + t2 + 616);\
                     mc309 = _mm_loadu_pd(ct##l + t2 + 618);\
                     mc310 = _mm_loadu_pd(ct##l + t2 + 620);\
                     mc311 = _mm_loadu_pd(ct##l + t2 + 622);\
                     mc312 = _mm_loadu_pd(ct##l + t2 + 624);\
                     mc313 = _mm_loadu_pd(ct##l + t2 + 626);\
                     mc314 = _mm_loadu_pd(ct##l + t2 + 628);\
                     mc315 = _mm_loadu_pd(ct##l + t2 + 630);\
                     mc316 = _mm_loadu_pd(ct##l + t2 + 632);\
                     mc317 = _mm_loadu_pd(ct##l + t2 + 634);\
                     mc318 = _mm_loadu_pd(ct##l + t2 + 636);\
                     mc319 = _mm_loadu_pd(ct##l + t2 + 638);\
                     mc320 = _mm_loadu_pd(ct##l + t2 + 640);\
                     mc321 = _mm_loadu_pd(ct##l + t2 + 642);\
                     mc322 = _mm_loadu_pd(ct##l + t2 + 644);\
                     mc323 = _mm_loadu_pd(ct##l + t2 + 646);\
                     mc324 = _mm_loadu_pd(ct##l + t2 + 648);\
                     mc325 = _mm_loadu_pd(ct##l + t2 + 650);\
                     mc326 = _mm_loadu_pd(ct##l + t2 + 652);\
                     mc327 = _mm_loadu_pd(ct##l + t2 + 654);\
                     mc328 = _mm_loadu_pd(ct##l + t2 + 656);\
                     mc329 = _mm_loadu_pd(ct##l + t2 + 658);\
                     mc330 = _mm_loadu_pd(ct##l + t2 + 660);\
                     mc331 = _mm_loadu_pd(ct##l + t2 + 662);\
                     mc332 = _mm_loadu_pd(ct##l + t2 + 664);\
                     mc333 = _mm_loadu_pd(ct##l + t2 + 666);\
                     mc334 = _mm_loadu_pd(ct##l + t2 + 668);\
                     mc335 = _mm_loadu_pd(ct##l + t2 + 670);\
                     mc336 = _mm_loadu_pd(ct##l + t2 + 672);\
                     mc337 = _mm_loadu_pd(ct##l + t2 + 674);\
                     mc338 = _mm_loadu_pd(ct##l + t2 + 676);\
                     mc339 = _mm_loadu_pd(ct##l + t2 + 678);\
                     mc340 = _mm_loadu_pd(ct##l + t2 + 680);\
                     mc341 = _mm_loadu_pd(ct##l + t2 + 682);\
                     mc342 = _mm_loadu_pd(ct##l + t2 + 684);\
                     mc343 = _mm_loadu_pd(ct##l + t2 + 686);\
                     mc344 = _mm_loadu_pd(ct##l + t2 + 688);\
                     mc345 = _mm_loadu_pd(ct##l + t2 + 690);\
                     mc346 = _mm_loadu_pd(ct##l + t2 + 692);\
                     mc347 = _mm_loadu_pd(ct##l + t2 + 694);\
                     mc348 = _mm_loadu_pd(ct##l + t2 + 696);\
                     mc349 = _mm_loadu_pd(ct##l + t2 + 698);\
                     mc350 = _mm_loadu_pd(ct##l + t2 + 700);\
                     mc351 = _mm_loadu_pd(ct##l + t2 + 702);\
                     mc352 = _mm_loadu_pd(ct##l + t2 + 704);\
                     mc353 = _mm_loadu_pd(ct##l + t2 + 706);\
                     mc354 = _mm_loadu_pd(ct##l + t2 + 708);\
                     mc355 = _mm_loadu_pd(ct##l + t2 + 710);\
                     mc356 = _mm_loadu_pd(ct##l + t2 + 712);\
                     mc357 = _mm_loadu_pd(ct##l + t2 + 714);\
                     mc358 = _mm_loadu_pd(ct##l + t2 + 716);\
                     mc359 = _mm_loadu_pd(ct##l + t2 + 718);\
                     mc360 = _mm_loadu_pd(ct##l + t2 + 720);\
                     mc361 = _mm_loadu_pd(ct##l + t2 + 722);\
                     mc362 = _mm_loadu_pd(ct##l + t2 + 724);\
                     mc363 = _mm_loadu_pd(ct##l + t2 + 726);\
                     mc364 = _mm_loadu_pd(ct##l + t2 + 728);\
                     mc365 = _mm_loadu_pd(ct##l + t2 + 730);\
                     mc366 = _mm_loadu_pd(ct##l + t2 + 732);\
                     mc367 = _mm_loadu_pd(ct##l + t2 + 734);\
                     mc368 = _mm_loadu_pd(ct##l + t2 + 736);\
                     mc369 = _mm_loadu_pd(ct##l + t2 + 738);\
                     mc370 = _mm_loadu_pd(ct##l + t2 + 740);\
                     mc371 = _mm_loadu_pd(ct##l + t2 + 742);\
                     mc372 = _mm_loadu_pd(ct##l + t2 + 744);\
                     mc373 = _mm_loadu_pd(ct##l + t2 + 746);\
                     mc374 = _mm_loadu_pd(ct##l + t2 + 748);\
                     mc375 = _mm_loadu_pd(ct##l + t2 + 750);\
                     mc376 = _mm_loadu_pd(ct##l + t2 + 752);\
                     mc377 = _mm_loadu_pd(ct##l + t2 + 754);\
                     mc378 = _mm_loadu_pd(ct##l + t2 + 756);\
                     mc379 = _mm_loadu_pd(ct##l + t2 + 758);\
                     mc380 = _mm_loadu_pd(ct##l + t2 + 760);\
                     mc381 = _mm_loadu_pd(ct##l + t2 + 762);\
                     mc382 = _mm_loadu_pd(ct##l + t2 + 764);\
                     mc383 = _mm_loadu_pd(ct##l + t2 + 766);\
                     mc384 = _mm_loadu_pd(ct##l + t2 + 768);\
                     mc385 = _mm_loadu_pd(ct##l + t2 + 770);\
                     mc386 = _mm_loadu_pd(ct##l + t2 + 772);\
                     mc387 = _mm_loadu_pd(ct##l + t2 + 774);\
                     mc388 = _mm_loadu_pd(ct##l + t2 + 776);\
                     mc389 = _mm_loadu_pd(ct##l + t2 + 778);\
                     mc390 = _mm_loadu_pd(ct##l + t2 + 780);\
                     mc391 = _mm_loadu_pd(ct##l + t2 + 782);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx8, mc8));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx9, mc9));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx10, mc10));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx11, mc11));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx12, mc12));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx13, mc13));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc15));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc16));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc17));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc18));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc19));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc20));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc21));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx8, mc22));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx9, mc23));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx10, mc24));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx11, mc25));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx12, mc26));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx13, mc27));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc28));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc29));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc30));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc31));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc32));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc33));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc34));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc35));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx8, mc36));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx9, mc37));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx10, mc38));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx11, mc39));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx12, mc40));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx13, mc41));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc42));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc43));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc44));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc45));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc46));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc47));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc48));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc49));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx8, mc50));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx9, mc51));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx10, mc52));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx11, mc53));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx12, mc54));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx13, mc55));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc56));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc57));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc58));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc59));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc60));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc61));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc62));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc63));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx8, mc64));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx9, mc65));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx10, mc66));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx11, mc67));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx12, mc68));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx13, mc69));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc70));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc71));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc72));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc73));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc74));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc75));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc76));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc77));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx8, mc78));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx9, mc79));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx10, mc80));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx11, mc81));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx12, mc82));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx13, mc83));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc84));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc85));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc86));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc87));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc88));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc89));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc90));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc91));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx8, mc92));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx9, mc93));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx10, mc94));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx11, mc95));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx12, mc96));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx13, mc97));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc98));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc99));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc100));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc101));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc102));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc103));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc104));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc105));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx8, mc106));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx9, mc107));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx10, mc108));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx11, mc109));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx12, mc110));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx13, mc111));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc112));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc113));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc114));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc115));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc116));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc117));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc118));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc119));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx8, mc120));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx9, mc121));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx10, mc122));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx11, mc123));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx12, mc124));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx13, mc125));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc126));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc127));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc128));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc129));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc130));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc131));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc132));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc133));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx8, mc134));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx9, mc135));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx10, mc136));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx11, mc137));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx12, mc138));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx13, mc139));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc140));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc141));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc142));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc143));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc144));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc145));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc146));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc147));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx8, mc148));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx9, mc149));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx10, mc150));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx11, mc151));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx12, mc152));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx13, mc153));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc154));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc155));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc156));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc157));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc158));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc159));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc160));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc161));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx8, mc162));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx9, mc163));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx10, mc164));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx11, mc165));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx12, mc166));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx13, mc167));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc168));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc169));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc170));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc171));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc172));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc173));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc174));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc175));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx8, mc176));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx9, mc177));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx10, mc178));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx11, mc179));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx12, mc180));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx13, mc181));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc182));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc183));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc184));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc185));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc186));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc187));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc188));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc189));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx8, mc190));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx9, mc191));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx10, mc192));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx11, mc193));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx12, mc194));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx13, mc195));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc196));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc197));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc198));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc199));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc200));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc201));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc202));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc203));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx8, mc204));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx9, mc205));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx10, mc206));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx11, mc207));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx12, mc208));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx13, mc209));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc210));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc211));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc212));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc213));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc214));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc215));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc216));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc217));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx8, mc218));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx9, mc219));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx10, mc220));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx11, mc221));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx12, mc222));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx13, mc223));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx0, mc224));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx1, mc225));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx2, mc226));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx3, mc227));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx4, mc228));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx5, mc229));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx6, mc230));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx7, mc231));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx8, mc232));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx9, mc233));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx10, mc234));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx11, mc235));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx12, mc236));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx13, mc237));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx0, mc238));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx1, mc239));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx2, mc240));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx3, mc241));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx4, mc242));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx5, mc243));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx6, mc244));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx7, mc245));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx8, mc246));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx9, mc247));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx10, mc248));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx11, mc249));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx12, mc250));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx13, mc251));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx0, mc252));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx1, mc253));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx2, mc254));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx3, mc255));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx4, mc256));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx5, mc257));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx6, mc258));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx7, mc259));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx8, mc260));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx9, mc261));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx10, mc262));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx11, mc263));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx12, mc264));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx13, mc265));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx0, mc266));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx1, mc267));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx2, mc268));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx3, mc269));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx4, mc270));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx5, mc271));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx6, mc272));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx7, mc273));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx8, mc274));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx9, mc275));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx10, mc276));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx11, mc277));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx12, mc278));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx13, mc279));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx0, mc280));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx1, mc281));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx2, mc282));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx3, mc283));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx4, mc284));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx5, mc285));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx6, mc286));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx7, mc287));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx8, mc288));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx9, mc289));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx10, mc290));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx11, mc291));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx12, mc292));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx13, mc293));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx0, mc294));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx1, mc295));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx2, mc296));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx3, mc297));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx4, mc298));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx5, mc299));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx6, mc300));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx7, mc301));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx8, mc302));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx9, mc303));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx10, mc304));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx11, mc305));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx12, mc306));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx13, mc307));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx0, mc308));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx1, mc309));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx2, mc310));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx3, mc311));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx4, mc312));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx5, mc313));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx6, mc314));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx7, mc315));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx8, mc316));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx9, mc317));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx10, mc318));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx11, mc319));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx12, mc320));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx13, mc321));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx0, mc322));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx1, mc323));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx2, mc324));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx3, mc325));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx4, mc326));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx5, mc327));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx6, mc328));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx7, mc329));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx8, mc330));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx9, mc331));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx10, mc332));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx11, mc333));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx12, mc334));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx13, mc335));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx0, mc336));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx1, mc337));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx2, mc338));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx3, mc339));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx4, mc340));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx5, mc341));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx6, mc342));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx7, mc343));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx8, mc344));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx9, mc345));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx10, mc346));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx11, mc347));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx12, mc348));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx13, mc349));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx0, mc350));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx1, mc351));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx2, mc352));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx3, mc353));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx4, mc354));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx5, mc355));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx6, mc356));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx7, mc357));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx8, mc358));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx9, mc359));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx10, mc360));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx11, mc361));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx12, mc362));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx13, mc363));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx0, mc364));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx1, mc365));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx2, mc366));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx3, mc367));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx4, mc368));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx5, mc369));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx6, mc370));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx7, mc371));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx8, mc372));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx9, mc373));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx10, mc374));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx11, mc375));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx12, mc376));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx13, mc377));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx0, mc378));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx1, mc379));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx2, mc380));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx3, mc381));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx4, mc382));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx5, mc383));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx6, mc384));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx7, mc385));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx8, mc386));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx9, mc387));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx10, mc388));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx11, mc389));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx12, mc390));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx13, mc391))

#define setup_28(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0);\
                         msum16 = _mm_set_pd(0,0);\
                         msum17 = _mm_set_pd(0,0);\
                         msum18 = _mm_set_pd(0,0);\
                         msum19 = _mm_set_pd(0,0);\
                         msum20 = _mm_set_pd(0,0);\
                         msum21 = _mm_set_pd(0,0);\
                         msum22 = _mm_set_pd(0,0);\
                         msum23 = _mm_set_pd(0,0);\
                         msum24 = _mm_set_pd(0,0);\
                         msum25 = _mm_set_pd(0,0);\
                         msum26 = _mm_set_pd(0,0);\
                         msum27 = _mm_set_pd(0,0)

#define save_28() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  msum16 = _mm_hadd_pd(msum16, msum17);\
                  msum18 = _mm_hadd_pd(msum18, msum19);\
                  msum20 = _mm_hadd_pd(msum20, msum21);\
                  msum22 = _mm_hadd_pd(msum22, msum23);\
                  msum24 = _mm_hadd_pd(msum24, msum25);\
                  msum26 = _mm_hadd_pd(msum26, msum27);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14);\
                  _mm_storeu_pd(y + t1 + 16,msum16);\
                  _mm_storeu_pd(y + t1 + 18,msum18);\
                  _mm_storeu_pd(y + t1 + 20,msum20);\
                  _mm_storeu_pd(y + t1 + 22,msum22);\
                  _mm_storeu_pd(y + t1 + 24,msum24);\
                  _mm_storeu_pd(y + t1 + 26,msum26)

PetscErrorCode BSG_MatMult_28(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8, mx9, mx10, mx11, mx12, mx13, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, mc128, mc129, mc130, mc131, mc132, mc133, mc134, mc135, mc136, mc137, mc138, mc139, mc140, mc141, mc142, mc143, mc144, mc145, mc146, mc147, mc148, mc149, mc150, mc151, mc152, mc153, mc154, mc155, mc156, mc157, mc158, mc159, mc160, mc161, mc162, mc163, mc164, mc165, mc166, mc167, mc168, mc169, mc170, mc171, mc172, mc173, mc174, mc175, mc176, mc177, mc178, mc179, mc180, mc181, mc182, mc183, mc184, mc185, mc186, mc187, mc188, mc189, mc190, mc191, mc192, mc193, mc194, mc195, mc196, mc197, mc198, mc199, mc200, mc201, mc202, mc203, mc204, mc205, mc206, mc207, mc208, mc209, mc210, mc211, mc212, mc213, mc214, mc215, mc216, mc217, mc218, mc219, mc220, mc221, mc222, mc223, mc224, mc225, mc226, mc227, mc228, mc229, mc230, mc231, mc232, mc233, mc234, mc235, mc236, mc237, mc238, mc239, mc240, mc241, mc242, mc243, mc244, mc245, mc246, mc247, mc248, mc249, mc250, mc251, mc252, mc253, mc254, mc255, mc256, mc257, mc258, mc259, mc260, mc261, mc262, mc263, mc264, mc265, mc266, mc267, mc268, mc269, mc270, mc271, mc272, mc273, mc274, mc275, mc276, mc277, mc278, mc279, mc280, mc281, mc282, mc283, mc284, mc285, mc286, mc287, mc288, mc289, mc290, mc291, mc292, mc293, mc294, mc295, mc296, mc297, mc298, mc299, mc300, mc301, mc302, mc303, mc304, mc305, mc306, mc307, mc308, mc309, mc310, mc311, mc312, mc313, mc314, mc315, mc316, mc317, mc318, mc319, mc320, mc321, mc322, mc323, mc324, mc325, mc326, mc327, mc328, mc329, mc330, mc331, mc332, mc333, mc334, mc335, mc336, mc337, mc338, mc339, mc340, mc341, mc342, mc343, mc344, mc345, mc346, mc347, mc348, mc349, mc350, mc351, mc352, mc353, mc354, mc355, mc356, mc357, mc358, mc359, mc360, mc361, mc362, mc363, mc364, mc365, mc366, mc367, mc368, mc369, mc370, mc371, mc372, mc373, mc374, mc375, mc376, mc377, mc378, mc379, mc380, mc381, mc382, mc383, mc384, mc385, mc386, mc387, mc388, mc389, mc390, mc391, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15, msum16, msum17, msum18, msum19, msum20, msum21, msum22, msum23, msum24, msum25, msum26, msum27;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_28(0);
            inline_28(3);
            inline_28(4);
            inline_28(5);
            inline_28(6);
            save_28();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_28(1);
            inline_28(2);
            inline_28(3);
            inline_28(4);
            inline_28(5);
            inline_28(6);
            save_28();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_28(lda3);
            inline_28(1);
            inline_28(2);
            inline_28(3);
            inline_28(4);
            inline_28(5);
            inline_28(6);
            save_28();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_28(lda2);
            inline_28(0);
            inline_28(1);
            inline_28(2);
            inline_28(3);
            inline_28(4);
            inline_28(5);
            inline_28(6);
            save_28();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_28(lda1 - lda2);
            inline_28(0);
            inline_28(1);
            inline_28(2);
            inline_28(3);
            inline_28(4);
            inline_28(5);
            save_28();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_28(lda1 - lda3);
            inline_28(0);
            inline_28(1);
            inline_28(2);
            inline_28(3);
            inline_28(4);
            save_28();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_28(lda1 - 1);
            inline_28(0);
            inline_28(1);
            inline_28(2);
            inline_28(3);
            save_28();
        }
    }

PetscFunctionReturn(0);
}

#define inline_30(l) \
                     mx0 = _mm_loadu_pd(xt##l + t1 + 0);\
                     mx1 = _mm_loadu_pd(xt##l + t1 + 2);\
                     mx2 = _mm_loadu_pd(xt##l + t1 + 4);\
                     mx3 = _mm_loadu_pd(xt##l + t1 + 6);\
                     mx4 = _mm_loadu_pd(xt##l + t1 + 8);\
                     mx5 = _mm_loadu_pd(xt##l + t1 + 10);\
                     mx6 = _mm_loadu_pd(xt##l + t1 + 12);\
                     mx7 = _mm_loadu_pd(xt##l + t1 + 14);\
                     mx8 = _mm_loadu_pd(xt##l + t1 + 16);\
                     mx9 = _mm_loadu_pd(xt##l + t1 + 18);\
                     mx10 = _mm_loadu_pd(xt##l + t1 + 20);\
                     mx11 = _mm_loadu_pd(xt##l + t1 + 22);\
                     mx12 = _mm_loadu_pd(xt##l + t1 + 24);\
                     mx13 = _mm_loadu_pd(xt##l + t1 + 26);\
                     mx14 = _mm_loadu_pd(xt##l + t1 + 28);\
                     mc0 = _mm_loadu_pd(ct##l + t2 + 0);\
                     mc1 = _mm_loadu_pd(ct##l + t2 + 2);\
                     mc2 = _mm_loadu_pd(ct##l + t2 + 4);\
                     mc3 = _mm_loadu_pd(ct##l + t2 + 6);\
                     mc4 = _mm_loadu_pd(ct##l + t2 + 8);\
                     mc5 = _mm_loadu_pd(ct##l + t2 + 10);\
                     mc6 = _mm_loadu_pd(ct##l + t2 + 12);\
                     mc7 = _mm_loadu_pd(ct##l + t2 + 14);\
                     mc8 = _mm_loadu_pd(ct##l + t2 + 16);\
                     mc9 = _mm_loadu_pd(ct##l + t2 + 18);\
                     mc10 = _mm_loadu_pd(ct##l + t2 + 20);\
                     mc11 = _mm_loadu_pd(ct##l + t2 + 22);\
                     mc12 = _mm_loadu_pd(ct##l + t2 + 24);\
                     mc13 = _mm_loadu_pd(ct##l + t2 + 26);\
                     mc14 = _mm_loadu_pd(ct##l + t2 + 28);\
                     mc15 = _mm_loadu_pd(ct##l + t2 + 30);\
                     mc16 = _mm_loadu_pd(ct##l + t2 + 32);\
                     mc17 = _mm_loadu_pd(ct##l + t2 + 34);\
                     mc18 = _mm_loadu_pd(ct##l + t2 + 36);\
                     mc19 = _mm_loadu_pd(ct##l + t2 + 38);\
                     mc20 = _mm_loadu_pd(ct##l + t2 + 40);\
                     mc21 = _mm_loadu_pd(ct##l + t2 + 42);\
                     mc22 = _mm_loadu_pd(ct##l + t2 + 44);\
                     mc23 = _mm_loadu_pd(ct##l + t2 + 46);\
                     mc24 = _mm_loadu_pd(ct##l + t2 + 48);\
                     mc25 = _mm_loadu_pd(ct##l + t2 + 50);\
                     mc26 = _mm_loadu_pd(ct##l + t2 + 52);\
                     mc27 = _mm_loadu_pd(ct##l + t2 + 54);\
                     mc28 = _mm_loadu_pd(ct##l + t2 + 56);\
                     mc29 = _mm_loadu_pd(ct##l + t2 + 58);\
                     mc30 = _mm_loadu_pd(ct##l + t2 + 60);\
                     mc31 = _mm_loadu_pd(ct##l + t2 + 62);\
                     mc32 = _mm_loadu_pd(ct##l + t2 + 64);\
                     mc33 = _mm_loadu_pd(ct##l + t2 + 66);\
                     mc34 = _mm_loadu_pd(ct##l + t2 + 68);\
                     mc35 = _mm_loadu_pd(ct##l + t2 + 70);\
                     mc36 = _mm_loadu_pd(ct##l + t2 + 72);\
                     mc37 = _mm_loadu_pd(ct##l + t2 + 74);\
                     mc38 = _mm_loadu_pd(ct##l + t2 + 76);\
                     mc39 = _mm_loadu_pd(ct##l + t2 + 78);\
                     mc40 = _mm_loadu_pd(ct##l + t2 + 80);\
                     mc41 = _mm_loadu_pd(ct##l + t2 + 82);\
                     mc42 = _mm_loadu_pd(ct##l + t2 + 84);\
                     mc43 = _mm_loadu_pd(ct##l + t2 + 86);\
                     mc44 = _mm_loadu_pd(ct##l + t2 + 88);\
                     mc45 = _mm_loadu_pd(ct##l + t2 + 90);\
                     mc46 = _mm_loadu_pd(ct##l + t2 + 92);\
                     mc47 = _mm_loadu_pd(ct##l + t2 + 94);\
                     mc48 = _mm_loadu_pd(ct##l + t2 + 96);\
                     mc49 = _mm_loadu_pd(ct##l + t2 + 98);\
                     mc50 = _mm_loadu_pd(ct##l + t2 + 100);\
                     mc51 = _mm_loadu_pd(ct##l + t2 + 102);\
                     mc52 = _mm_loadu_pd(ct##l + t2 + 104);\
                     mc53 = _mm_loadu_pd(ct##l + t2 + 106);\
                     mc54 = _mm_loadu_pd(ct##l + t2 + 108);\
                     mc55 = _mm_loadu_pd(ct##l + t2 + 110);\
                     mc56 = _mm_loadu_pd(ct##l + t2 + 112);\
                     mc57 = _mm_loadu_pd(ct##l + t2 + 114);\
                     mc58 = _mm_loadu_pd(ct##l + t2 + 116);\
                     mc59 = _mm_loadu_pd(ct##l + t2 + 118);\
                     mc60 = _mm_loadu_pd(ct##l + t2 + 120);\
                     mc61 = _mm_loadu_pd(ct##l + t2 + 122);\
                     mc62 = _mm_loadu_pd(ct##l + t2 + 124);\
                     mc63 = _mm_loadu_pd(ct##l + t2 + 126);\
                     mc64 = _mm_loadu_pd(ct##l + t2 + 128);\
                     mc65 = _mm_loadu_pd(ct##l + t2 + 130);\
                     mc66 = _mm_loadu_pd(ct##l + t2 + 132);\
                     mc67 = _mm_loadu_pd(ct##l + t2 + 134);\
                     mc68 = _mm_loadu_pd(ct##l + t2 + 136);\
                     mc69 = _mm_loadu_pd(ct##l + t2 + 138);\
                     mc70 = _mm_loadu_pd(ct##l + t2 + 140);\
                     mc71 = _mm_loadu_pd(ct##l + t2 + 142);\
                     mc72 = _mm_loadu_pd(ct##l + t2 + 144);\
                     mc73 = _mm_loadu_pd(ct##l + t2 + 146);\
                     mc74 = _mm_loadu_pd(ct##l + t2 + 148);\
                     mc75 = _mm_loadu_pd(ct##l + t2 + 150);\
                     mc76 = _mm_loadu_pd(ct##l + t2 + 152);\
                     mc77 = _mm_loadu_pd(ct##l + t2 + 154);\
                     mc78 = _mm_loadu_pd(ct##l + t2 + 156);\
                     mc79 = _mm_loadu_pd(ct##l + t2 + 158);\
                     mc80 = _mm_loadu_pd(ct##l + t2 + 160);\
                     mc81 = _mm_loadu_pd(ct##l + t2 + 162);\
                     mc82 = _mm_loadu_pd(ct##l + t2 + 164);\
                     mc83 = _mm_loadu_pd(ct##l + t2 + 166);\
                     mc84 = _mm_loadu_pd(ct##l + t2 + 168);\
                     mc85 = _mm_loadu_pd(ct##l + t2 + 170);\
                     mc86 = _mm_loadu_pd(ct##l + t2 + 172);\
                     mc87 = _mm_loadu_pd(ct##l + t2 + 174);\
                     mc88 = _mm_loadu_pd(ct##l + t2 + 176);\
                     mc89 = _mm_loadu_pd(ct##l + t2 + 178);\
                     mc90 = _mm_loadu_pd(ct##l + t2 + 180);\
                     mc91 = _mm_loadu_pd(ct##l + t2 + 182);\
                     mc92 = _mm_loadu_pd(ct##l + t2 + 184);\
                     mc93 = _mm_loadu_pd(ct##l + t2 + 186);\
                     mc94 = _mm_loadu_pd(ct##l + t2 + 188);\
                     mc95 = _mm_loadu_pd(ct##l + t2 + 190);\
                     mc96 = _mm_loadu_pd(ct##l + t2 + 192);\
                     mc97 = _mm_loadu_pd(ct##l + t2 + 194);\
                     mc98 = _mm_loadu_pd(ct##l + t2 + 196);\
                     mc99 = _mm_loadu_pd(ct##l + t2 + 198);\
                     mc100 = _mm_loadu_pd(ct##l + t2 + 200);\
                     mc101 = _mm_loadu_pd(ct##l + t2 + 202);\
                     mc102 = _mm_loadu_pd(ct##l + t2 + 204);\
                     mc103 = _mm_loadu_pd(ct##l + t2 + 206);\
                     mc104 = _mm_loadu_pd(ct##l + t2 + 208);\
                     mc105 = _mm_loadu_pd(ct##l + t2 + 210);\
                     mc106 = _mm_loadu_pd(ct##l + t2 + 212);\
                     mc107 = _mm_loadu_pd(ct##l + t2 + 214);\
                     mc108 = _mm_loadu_pd(ct##l + t2 + 216);\
                     mc109 = _mm_loadu_pd(ct##l + t2 + 218);\
                     mc110 = _mm_loadu_pd(ct##l + t2 + 220);\
                     mc111 = _mm_loadu_pd(ct##l + t2 + 222);\
                     mc112 = _mm_loadu_pd(ct##l + t2 + 224);\
                     mc113 = _mm_loadu_pd(ct##l + t2 + 226);\
                     mc114 = _mm_loadu_pd(ct##l + t2 + 228);\
                     mc115 = _mm_loadu_pd(ct##l + t2 + 230);\
                     mc116 = _mm_loadu_pd(ct##l + t2 + 232);\
                     mc117 = _mm_loadu_pd(ct##l + t2 + 234);\
                     mc118 = _mm_loadu_pd(ct##l + t2 + 236);\
                     mc119 = _mm_loadu_pd(ct##l + t2 + 238);\
                     mc120 = _mm_loadu_pd(ct##l + t2 + 240);\
                     mc121 = _mm_loadu_pd(ct##l + t2 + 242);\
                     mc122 = _mm_loadu_pd(ct##l + t2 + 244);\
                     mc123 = _mm_loadu_pd(ct##l + t2 + 246);\
                     mc124 = _mm_loadu_pd(ct##l + t2 + 248);\
                     mc125 = _mm_loadu_pd(ct##l + t2 + 250);\
                     mc126 = _mm_loadu_pd(ct##l + t2 + 252);\
                     mc127 = _mm_loadu_pd(ct##l + t2 + 254);\
                     mc128 = _mm_loadu_pd(ct##l + t2 + 256);\
                     mc129 = _mm_loadu_pd(ct##l + t2 + 258);\
                     mc130 = _mm_loadu_pd(ct##l + t2 + 260);\
                     mc131 = _mm_loadu_pd(ct##l + t2 + 262);\
                     mc132 = _mm_loadu_pd(ct##l + t2 + 264);\
                     mc133 = _mm_loadu_pd(ct##l + t2 + 266);\
                     mc134 = _mm_loadu_pd(ct##l + t2 + 268);\
                     mc135 = _mm_loadu_pd(ct##l + t2 + 270);\
                     mc136 = _mm_loadu_pd(ct##l + t2 + 272);\
                     mc137 = _mm_loadu_pd(ct##l + t2 + 274);\
                     mc138 = _mm_loadu_pd(ct##l + t2 + 276);\
                     mc139 = _mm_loadu_pd(ct##l + t2 + 278);\
                     mc140 = _mm_loadu_pd(ct##l + t2 + 280);\
                     mc141 = _mm_loadu_pd(ct##l + t2 + 282);\
                     mc142 = _mm_loadu_pd(ct##l + t2 + 284);\
                     mc143 = _mm_loadu_pd(ct##l + t2 + 286);\
                     mc144 = _mm_loadu_pd(ct##l + t2 + 288);\
                     mc145 = _mm_loadu_pd(ct##l + t2 + 290);\
                     mc146 = _mm_loadu_pd(ct##l + t2 + 292);\
                     mc147 = _mm_loadu_pd(ct##l + t2 + 294);\
                     mc148 = _mm_loadu_pd(ct##l + t2 + 296);\
                     mc149 = _mm_loadu_pd(ct##l + t2 + 298);\
                     mc150 = _mm_loadu_pd(ct##l + t2 + 300);\
                     mc151 = _mm_loadu_pd(ct##l + t2 + 302);\
                     mc152 = _mm_loadu_pd(ct##l + t2 + 304);\
                     mc153 = _mm_loadu_pd(ct##l + t2 + 306);\
                     mc154 = _mm_loadu_pd(ct##l + t2 + 308);\
                     mc155 = _mm_loadu_pd(ct##l + t2 + 310);\
                     mc156 = _mm_loadu_pd(ct##l + t2 + 312);\
                     mc157 = _mm_loadu_pd(ct##l + t2 + 314);\
                     mc158 = _mm_loadu_pd(ct##l + t2 + 316);\
                     mc159 = _mm_loadu_pd(ct##l + t2 + 318);\
                     mc160 = _mm_loadu_pd(ct##l + t2 + 320);\
                     mc161 = _mm_loadu_pd(ct##l + t2 + 322);\
                     mc162 = _mm_loadu_pd(ct##l + t2 + 324);\
                     mc163 = _mm_loadu_pd(ct##l + t2 + 326);\
                     mc164 = _mm_loadu_pd(ct##l + t2 + 328);\
                     mc165 = _mm_loadu_pd(ct##l + t2 + 330);\
                     mc166 = _mm_loadu_pd(ct##l + t2 + 332);\
                     mc167 = _mm_loadu_pd(ct##l + t2 + 334);\
                     mc168 = _mm_loadu_pd(ct##l + t2 + 336);\
                     mc169 = _mm_loadu_pd(ct##l + t2 + 338);\
                     mc170 = _mm_loadu_pd(ct##l + t2 + 340);\
                     mc171 = _mm_loadu_pd(ct##l + t2 + 342);\
                     mc172 = _mm_loadu_pd(ct##l + t2 + 344);\
                     mc173 = _mm_loadu_pd(ct##l + t2 + 346);\
                     mc174 = _mm_loadu_pd(ct##l + t2 + 348);\
                     mc175 = _mm_loadu_pd(ct##l + t2 + 350);\
                     mc176 = _mm_loadu_pd(ct##l + t2 + 352);\
                     mc177 = _mm_loadu_pd(ct##l + t2 + 354);\
                     mc178 = _mm_loadu_pd(ct##l + t2 + 356);\
                     mc179 = _mm_loadu_pd(ct##l + t2 + 358);\
                     mc180 = _mm_loadu_pd(ct##l + t2 + 360);\
                     mc181 = _mm_loadu_pd(ct##l + t2 + 362);\
                     mc182 = _mm_loadu_pd(ct##l + t2 + 364);\
                     mc183 = _mm_loadu_pd(ct##l + t2 + 366);\
                     mc184 = _mm_loadu_pd(ct##l + t2 + 368);\
                     mc185 = _mm_loadu_pd(ct##l + t2 + 370);\
                     mc186 = _mm_loadu_pd(ct##l + t2 + 372);\
                     mc187 = _mm_loadu_pd(ct##l + t2 + 374);\
                     mc188 = _mm_loadu_pd(ct##l + t2 + 376);\
                     mc189 = _mm_loadu_pd(ct##l + t2 + 378);\
                     mc190 = _mm_loadu_pd(ct##l + t2 + 380);\
                     mc191 = _mm_loadu_pd(ct##l + t2 + 382);\
                     mc192 = _mm_loadu_pd(ct##l + t2 + 384);\
                     mc193 = _mm_loadu_pd(ct##l + t2 + 386);\
                     mc194 = _mm_loadu_pd(ct##l + t2 + 388);\
                     mc195 = _mm_loadu_pd(ct##l + t2 + 390);\
                     mc196 = _mm_loadu_pd(ct##l + t2 + 392);\
                     mc197 = _mm_loadu_pd(ct##l + t2 + 394);\
                     mc198 = _mm_loadu_pd(ct##l + t2 + 396);\
                     mc199 = _mm_loadu_pd(ct##l + t2 + 398);\
                     mc200 = _mm_loadu_pd(ct##l + t2 + 400);\
                     mc201 = _mm_loadu_pd(ct##l + t2 + 402);\
                     mc202 = _mm_loadu_pd(ct##l + t2 + 404);\
                     mc203 = _mm_loadu_pd(ct##l + t2 + 406);\
                     mc204 = _mm_loadu_pd(ct##l + t2 + 408);\
                     mc205 = _mm_loadu_pd(ct##l + t2 + 410);\
                     mc206 = _mm_loadu_pd(ct##l + t2 + 412);\
                     mc207 = _mm_loadu_pd(ct##l + t2 + 414);\
                     mc208 = _mm_loadu_pd(ct##l + t2 + 416);\
                     mc209 = _mm_loadu_pd(ct##l + t2 + 418);\
                     mc210 = _mm_loadu_pd(ct##l + t2 + 420);\
                     mc211 = _mm_loadu_pd(ct##l + t2 + 422);\
                     mc212 = _mm_loadu_pd(ct##l + t2 + 424);\
                     mc213 = _mm_loadu_pd(ct##l + t2 + 426);\
                     mc214 = _mm_loadu_pd(ct##l + t2 + 428);\
                     mc215 = _mm_loadu_pd(ct##l + t2 + 430);\
                     mc216 = _mm_loadu_pd(ct##l + t2 + 432);\
                     mc217 = _mm_loadu_pd(ct##l + t2 + 434);\
                     mc218 = _mm_loadu_pd(ct##l + t2 + 436);\
                     mc219 = _mm_loadu_pd(ct##l + t2 + 438);\
                     mc220 = _mm_loadu_pd(ct##l + t2 + 440);\
                     mc221 = _mm_loadu_pd(ct##l + t2 + 442);\
                     mc222 = _mm_loadu_pd(ct##l + t2 + 444);\
                     mc223 = _mm_loadu_pd(ct##l + t2 + 446);\
                     mc224 = _mm_loadu_pd(ct##l + t2 + 448);\
                     mc225 = _mm_loadu_pd(ct##l + t2 + 450);\
                     mc226 = _mm_loadu_pd(ct##l + t2 + 452);\
                     mc227 = _mm_loadu_pd(ct##l + t2 + 454);\
                     mc228 = _mm_loadu_pd(ct##l + t2 + 456);\
                     mc229 = _mm_loadu_pd(ct##l + t2 + 458);\
                     mc230 = _mm_loadu_pd(ct##l + t2 + 460);\
                     mc231 = _mm_loadu_pd(ct##l + t2 + 462);\
                     mc232 = _mm_loadu_pd(ct##l + t2 + 464);\
                     mc233 = _mm_loadu_pd(ct##l + t2 + 466);\
                     mc234 = _mm_loadu_pd(ct##l + t2 + 468);\
                     mc235 = _mm_loadu_pd(ct##l + t2 + 470);\
                     mc236 = _mm_loadu_pd(ct##l + t2 + 472);\
                     mc237 = _mm_loadu_pd(ct##l + t2 + 474);\
                     mc238 = _mm_loadu_pd(ct##l + t2 + 476);\
                     mc239 = _mm_loadu_pd(ct##l + t2 + 478);\
                     mc240 = _mm_loadu_pd(ct##l + t2 + 480);\
                     mc241 = _mm_loadu_pd(ct##l + t2 + 482);\
                     mc242 = _mm_loadu_pd(ct##l + t2 + 484);\
                     mc243 = _mm_loadu_pd(ct##l + t2 + 486);\
                     mc244 = _mm_loadu_pd(ct##l + t2 + 488);\
                     mc245 = _mm_loadu_pd(ct##l + t2 + 490);\
                     mc246 = _mm_loadu_pd(ct##l + t2 + 492);\
                     mc247 = _mm_loadu_pd(ct##l + t2 + 494);\
                     mc248 = _mm_loadu_pd(ct##l + t2 + 496);\
                     mc249 = _mm_loadu_pd(ct##l + t2 + 498);\
                     mc250 = _mm_loadu_pd(ct##l + t2 + 500);\
                     mc251 = _mm_loadu_pd(ct##l + t2 + 502);\
                     mc252 = _mm_loadu_pd(ct##l + t2 + 504);\
                     mc253 = _mm_loadu_pd(ct##l + t2 + 506);\
                     mc254 = _mm_loadu_pd(ct##l + t2 + 508);\
                     mc255 = _mm_loadu_pd(ct##l + t2 + 510);\
                     mc256 = _mm_loadu_pd(ct##l + t2 + 512);\
                     mc257 = _mm_loadu_pd(ct##l + t2 + 514);\
                     mc258 = _mm_loadu_pd(ct##l + t2 + 516);\
                     mc259 = _mm_loadu_pd(ct##l + t2 + 518);\
                     mc260 = _mm_loadu_pd(ct##l + t2 + 520);\
                     mc261 = _mm_loadu_pd(ct##l + t2 + 522);\
                     mc262 = _mm_loadu_pd(ct##l + t2 + 524);\
                     mc263 = _mm_loadu_pd(ct##l + t2 + 526);\
                     mc264 = _mm_loadu_pd(ct##l + t2 + 528);\
                     mc265 = _mm_loadu_pd(ct##l + t2 + 530);\
                     mc266 = _mm_loadu_pd(ct##l + t2 + 532);\
                     mc267 = _mm_loadu_pd(ct##l + t2 + 534);\
                     mc268 = _mm_loadu_pd(ct##l + t2 + 536);\
                     mc269 = _mm_loadu_pd(ct##l + t2 + 538);\
                     mc270 = _mm_loadu_pd(ct##l + t2 + 540);\
                     mc271 = _mm_loadu_pd(ct##l + t2 + 542);\
                     mc272 = _mm_loadu_pd(ct##l + t2 + 544);\
                     mc273 = _mm_loadu_pd(ct##l + t2 + 546);\
                     mc274 = _mm_loadu_pd(ct##l + t2 + 548);\
                     mc275 = _mm_loadu_pd(ct##l + t2 + 550);\
                     mc276 = _mm_loadu_pd(ct##l + t2 + 552);\
                     mc277 = _mm_loadu_pd(ct##l + t2 + 554);\
                     mc278 = _mm_loadu_pd(ct##l + t2 + 556);\
                     mc279 = _mm_loadu_pd(ct##l + t2 + 558);\
                     mc280 = _mm_loadu_pd(ct##l + t2 + 560);\
                     mc281 = _mm_loadu_pd(ct##l + t2 + 562);\
                     mc282 = _mm_loadu_pd(ct##l + t2 + 564);\
                     mc283 = _mm_loadu_pd(ct##l + t2 + 566);\
                     mc284 = _mm_loadu_pd(ct##l + t2 + 568);\
                     mc285 = _mm_loadu_pd(ct##l + t2 + 570);\
                     mc286 = _mm_loadu_pd(ct##l + t2 + 572);\
                     mc287 = _mm_loadu_pd(ct##l + t2 + 574);\
                     mc288 = _mm_loadu_pd(ct##l + t2 + 576);\
                     mc289 = _mm_loadu_pd(ct##l + t2 + 578);\
                     mc290 = _mm_loadu_pd(ct##l + t2 + 580);\
                     mc291 = _mm_loadu_pd(ct##l + t2 + 582);\
                     mc292 = _mm_loadu_pd(ct##l + t2 + 584);\
                     mc293 = _mm_loadu_pd(ct##l + t2 + 586);\
                     mc294 = _mm_loadu_pd(ct##l + t2 + 588);\
                     mc295 = _mm_loadu_pd(ct##l + t2 + 590);\
                     mc296 = _mm_loadu_pd(ct##l + t2 + 592);\
                     mc297 = _mm_loadu_pd(ct##l + t2 + 594);\
                     mc298 = _mm_loadu_pd(ct##l + t2 + 596);\
                     mc299 = _mm_loadu_pd(ct##l + t2 + 598);\
                     mc300 = _mm_loadu_pd(ct##l + t2 + 600);\
                     mc301 = _mm_loadu_pd(ct##l + t2 + 602);\
                     mc302 = _mm_loadu_pd(ct##l + t2 + 604);\
                     mc303 = _mm_loadu_pd(ct##l + t2 + 606);\
                     mc304 = _mm_loadu_pd(ct##l + t2 + 608);\
                     mc305 = _mm_loadu_pd(ct##l + t2 + 610);\
                     mc306 = _mm_loadu_pd(ct##l + t2 + 612);\
                     mc307 = _mm_loadu_pd(ct##l + t2 + 614);\
                     mc308 = _mm_loadu_pd(ct##l + t2 + 616);\
                     mc309 = _mm_loadu_pd(ct##l + t2 + 618);\
                     mc310 = _mm_loadu_pd(ct##l + t2 + 620);\
                     mc311 = _mm_loadu_pd(ct##l + t2 + 622);\
                     mc312 = _mm_loadu_pd(ct##l + t2 + 624);\
                     mc313 = _mm_loadu_pd(ct##l + t2 + 626);\
                     mc314 = _mm_loadu_pd(ct##l + t2 + 628);\
                     mc315 = _mm_loadu_pd(ct##l + t2 + 630);\
                     mc316 = _mm_loadu_pd(ct##l + t2 + 632);\
                     mc317 = _mm_loadu_pd(ct##l + t2 + 634);\
                     mc318 = _mm_loadu_pd(ct##l + t2 + 636);\
                     mc319 = _mm_loadu_pd(ct##l + t2 + 638);\
                     mc320 = _mm_loadu_pd(ct##l + t2 + 640);\
                     mc321 = _mm_loadu_pd(ct##l + t2 + 642);\
                     mc322 = _mm_loadu_pd(ct##l + t2 + 644);\
                     mc323 = _mm_loadu_pd(ct##l + t2 + 646);\
                     mc324 = _mm_loadu_pd(ct##l + t2 + 648);\
                     mc325 = _mm_loadu_pd(ct##l + t2 + 650);\
                     mc326 = _mm_loadu_pd(ct##l + t2 + 652);\
                     mc327 = _mm_loadu_pd(ct##l + t2 + 654);\
                     mc328 = _mm_loadu_pd(ct##l + t2 + 656);\
                     mc329 = _mm_loadu_pd(ct##l + t2 + 658);\
                     mc330 = _mm_loadu_pd(ct##l + t2 + 660);\
                     mc331 = _mm_loadu_pd(ct##l + t2 + 662);\
                     mc332 = _mm_loadu_pd(ct##l + t2 + 664);\
                     mc333 = _mm_loadu_pd(ct##l + t2 + 666);\
                     mc334 = _mm_loadu_pd(ct##l + t2 + 668);\
                     mc335 = _mm_loadu_pd(ct##l + t2 + 670);\
                     mc336 = _mm_loadu_pd(ct##l + t2 + 672);\
                     mc337 = _mm_loadu_pd(ct##l + t2 + 674);\
                     mc338 = _mm_loadu_pd(ct##l + t2 + 676);\
                     mc339 = _mm_loadu_pd(ct##l + t2 + 678);\
                     mc340 = _mm_loadu_pd(ct##l + t2 + 680);\
                     mc341 = _mm_loadu_pd(ct##l + t2 + 682);\
                     mc342 = _mm_loadu_pd(ct##l + t2 + 684);\
                     mc343 = _mm_loadu_pd(ct##l + t2 + 686);\
                     mc344 = _mm_loadu_pd(ct##l + t2 + 688);\
                     mc345 = _mm_loadu_pd(ct##l + t2 + 690);\
                     mc346 = _mm_loadu_pd(ct##l + t2 + 692);\
                     mc347 = _mm_loadu_pd(ct##l + t2 + 694);\
                     mc348 = _mm_loadu_pd(ct##l + t2 + 696);\
                     mc349 = _mm_loadu_pd(ct##l + t2 + 698);\
                     mc350 = _mm_loadu_pd(ct##l + t2 + 700);\
                     mc351 = _mm_loadu_pd(ct##l + t2 + 702);\
                     mc352 = _mm_loadu_pd(ct##l + t2 + 704);\
                     mc353 = _mm_loadu_pd(ct##l + t2 + 706);\
                     mc354 = _mm_loadu_pd(ct##l + t2 + 708);\
                     mc355 = _mm_loadu_pd(ct##l + t2 + 710);\
                     mc356 = _mm_loadu_pd(ct##l + t2 + 712);\
                     mc357 = _mm_loadu_pd(ct##l + t2 + 714);\
                     mc358 = _mm_loadu_pd(ct##l + t2 + 716);\
                     mc359 = _mm_loadu_pd(ct##l + t2 + 718);\
                     mc360 = _mm_loadu_pd(ct##l + t2 + 720);\
                     mc361 = _mm_loadu_pd(ct##l + t2 + 722);\
                     mc362 = _mm_loadu_pd(ct##l + t2 + 724);\
                     mc363 = _mm_loadu_pd(ct##l + t2 + 726);\
                     mc364 = _mm_loadu_pd(ct##l + t2 + 728);\
                     mc365 = _mm_loadu_pd(ct##l + t2 + 730);\
                     mc366 = _mm_loadu_pd(ct##l + t2 + 732);\
                     mc367 = _mm_loadu_pd(ct##l + t2 + 734);\
                     mc368 = _mm_loadu_pd(ct##l + t2 + 736);\
                     mc369 = _mm_loadu_pd(ct##l + t2 + 738);\
                     mc370 = _mm_loadu_pd(ct##l + t2 + 740);\
                     mc371 = _mm_loadu_pd(ct##l + t2 + 742);\
                     mc372 = _mm_loadu_pd(ct##l + t2 + 744);\
                     mc373 = _mm_loadu_pd(ct##l + t2 + 746);\
                     mc374 = _mm_loadu_pd(ct##l + t2 + 748);\
                     mc375 = _mm_loadu_pd(ct##l + t2 + 750);\
                     mc376 = _mm_loadu_pd(ct##l + t2 + 752);\
                     mc377 = _mm_loadu_pd(ct##l + t2 + 754);\
                     mc378 = _mm_loadu_pd(ct##l + t2 + 756);\
                     mc379 = _mm_loadu_pd(ct##l + t2 + 758);\
                     mc380 = _mm_loadu_pd(ct##l + t2 + 760);\
                     mc381 = _mm_loadu_pd(ct##l + t2 + 762);\
                     mc382 = _mm_loadu_pd(ct##l + t2 + 764);\
                     mc383 = _mm_loadu_pd(ct##l + t2 + 766);\
                     mc384 = _mm_loadu_pd(ct##l + t2 + 768);\
                     mc385 = _mm_loadu_pd(ct##l + t2 + 770);\
                     mc386 = _mm_loadu_pd(ct##l + t2 + 772);\
                     mc387 = _mm_loadu_pd(ct##l + t2 + 774);\
                     mc388 = _mm_loadu_pd(ct##l + t2 + 776);\
                     mc389 = _mm_loadu_pd(ct##l + t2 + 778);\
                     mc390 = _mm_loadu_pd(ct##l + t2 + 780);\
                     mc391 = _mm_loadu_pd(ct##l + t2 + 782);\
                     mc392 = _mm_loadu_pd(ct##l + t2 + 784);\
                     mc393 = _mm_loadu_pd(ct##l + t2 + 786);\
                     mc394 = _mm_loadu_pd(ct##l + t2 + 788);\
                     mc395 = _mm_loadu_pd(ct##l + t2 + 790);\
                     mc396 = _mm_loadu_pd(ct##l + t2 + 792);\
                     mc397 = _mm_loadu_pd(ct##l + t2 + 794);\
                     mc398 = _mm_loadu_pd(ct##l + t2 + 796);\
                     mc399 = _mm_loadu_pd(ct##l + t2 + 798);\
                     mc400 = _mm_loadu_pd(ct##l + t2 + 800);\
                     mc401 = _mm_loadu_pd(ct##l + t2 + 802);\
                     mc402 = _mm_loadu_pd(ct##l + t2 + 804);\
                     mc403 = _mm_loadu_pd(ct##l + t2 + 806);\
                     mc404 = _mm_loadu_pd(ct##l + t2 + 808);\
                     mc405 = _mm_loadu_pd(ct##l + t2 + 810);\
                     mc406 = _mm_loadu_pd(ct##l + t2 + 812);\
                     mc407 = _mm_loadu_pd(ct##l + t2 + 814);\
                     mc408 = _mm_loadu_pd(ct##l + t2 + 816);\
                     mc409 = _mm_loadu_pd(ct##l + t2 + 818);\
                     mc410 = _mm_loadu_pd(ct##l + t2 + 820);\
                     mc411 = _mm_loadu_pd(ct##l + t2 + 822);\
                     mc412 = _mm_loadu_pd(ct##l + t2 + 824);\
                     mc413 = _mm_loadu_pd(ct##l + t2 + 826);\
                     mc414 = _mm_loadu_pd(ct##l + t2 + 828);\
                     mc415 = _mm_loadu_pd(ct##l + t2 + 830);\
                     mc416 = _mm_loadu_pd(ct##l + t2 + 832);\
                     mc417 = _mm_loadu_pd(ct##l + t2 + 834);\
                     mc418 = _mm_loadu_pd(ct##l + t2 + 836);\
                     mc419 = _mm_loadu_pd(ct##l + t2 + 838);\
                     mc420 = _mm_loadu_pd(ct##l + t2 + 840);\
                     mc421 = _mm_loadu_pd(ct##l + t2 + 842);\
                     mc422 = _mm_loadu_pd(ct##l + t2 + 844);\
                     mc423 = _mm_loadu_pd(ct##l + t2 + 846);\
                     mc424 = _mm_loadu_pd(ct##l + t2 + 848);\
                     mc425 = _mm_loadu_pd(ct##l + t2 + 850);\
                     mc426 = _mm_loadu_pd(ct##l + t2 + 852);\
                     mc427 = _mm_loadu_pd(ct##l + t2 + 854);\
                     mc428 = _mm_loadu_pd(ct##l + t2 + 856);\
                     mc429 = _mm_loadu_pd(ct##l + t2 + 858);\
                     mc430 = _mm_loadu_pd(ct##l + t2 + 860);\
                     mc431 = _mm_loadu_pd(ct##l + t2 + 862);\
                     mc432 = _mm_loadu_pd(ct##l + t2 + 864);\
                     mc433 = _mm_loadu_pd(ct##l + t2 + 866);\
                     mc434 = _mm_loadu_pd(ct##l + t2 + 868);\
                     mc435 = _mm_loadu_pd(ct##l + t2 + 870);\
                     mc436 = _mm_loadu_pd(ct##l + t2 + 872);\
                     mc437 = _mm_loadu_pd(ct##l + t2 + 874);\
                     mc438 = _mm_loadu_pd(ct##l + t2 + 876);\
                     mc439 = _mm_loadu_pd(ct##l + t2 + 878);\
                     mc440 = _mm_loadu_pd(ct##l + t2 + 880);\
                     mc441 = _mm_loadu_pd(ct##l + t2 + 882);\
                     mc442 = _mm_loadu_pd(ct##l + t2 + 884);\
                     mc443 = _mm_loadu_pd(ct##l + t2 + 886);\
                     mc444 = _mm_loadu_pd(ct##l + t2 + 888);\
                     mc445 = _mm_loadu_pd(ct##l + t2 + 890);\
                     mc446 = _mm_loadu_pd(ct##l + t2 + 892);\
                     mc447 = _mm_loadu_pd(ct##l + t2 + 894);\
                     mc448 = _mm_loadu_pd(ct##l + t2 + 896);\
                     mc449 = _mm_loadu_pd(ct##l + t2 + 898);\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx3, mc3));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx4, mc4));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx5, mc5));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx6, mc6));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx7, mc7));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx8, mc8));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx9, mc9));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx10, mc10));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx11, mc11));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx12, mc12));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx13, mc13));\
                     msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx14, mc14));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc15));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc16));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc17));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx3, mc18));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx4, mc19));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx5, mc20));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx6, mc21));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx7, mc22));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx8, mc23));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx9, mc24));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx10, mc25));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx11, mc26));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx12, mc27));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx13, mc28));\
                     msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx14, mc29));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc30));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc31));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc32));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx3, mc33));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx4, mc34));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx5, mc35));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx6, mc36));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx7, mc37));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx8, mc38));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx9, mc39));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx10, mc40));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx11, mc41));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx12, mc42));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx13, mc43));\
                     msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx14, mc44));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc45));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc46));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc47));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx3, mc48));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx4, mc49));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx5, mc50));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx6, mc51));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx7, mc52));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx8, mc53));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx9, mc54));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx10, mc55));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx11, mc56));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx12, mc57));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx13, mc58));\
                     msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx14, mc59));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc60));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc61));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc62));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx3, mc63));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx4, mc64));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx5, mc65));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx6, mc66));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx7, mc67));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx8, mc68));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx9, mc69));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx10, mc70));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx11, mc71));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx12, mc72));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx13, mc73));\
                     msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx14, mc74));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc75));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc76));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc77));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx3, mc78));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx4, mc79));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx5, mc80));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx6, mc81));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx7, mc82));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx8, mc83));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx9, mc84));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx10, mc85));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx11, mc86));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx12, mc87));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx13, mc88));\
                     msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx14, mc89));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx0, mc90));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx1, mc91));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx2, mc92));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx3, mc93));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx4, mc94));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx5, mc95));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx6, mc96));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx7, mc97));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx8, mc98));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx9, mc99));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx10, mc100));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx11, mc101));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx12, mc102));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx13, mc103));\
                     msum6 = _mm_add_pd(msum6 , _mm_mul_pd(mx14, mc104));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx0, mc105));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx1, mc106));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx2, mc107));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx3, mc108));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx4, mc109));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx5, mc110));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx6, mc111));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx7, mc112));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx8, mc113));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx9, mc114));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx10, mc115));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx11, mc116));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx12, mc117));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx13, mc118));\
                     msum7 = _mm_add_pd(msum7 , _mm_mul_pd(mx14, mc119));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx0, mc120));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx1, mc121));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx2, mc122));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx3, mc123));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx4, mc124));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx5, mc125));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx6, mc126));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx7, mc127));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx8, mc128));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx9, mc129));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx10, mc130));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx11, mc131));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx12, mc132));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx13, mc133));\
                     msum8 = _mm_add_pd(msum8 , _mm_mul_pd(mx14, mc134));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx0, mc135));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx1, mc136));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx2, mc137));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx3, mc138));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx4, mc139));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx5, mc140));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx6, mc141));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx7, mc142));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx8, mc143));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx9, mc144));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx10, mc145));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx11, mc146));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx12, mc147));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx13, mc148));\
                     msum9 = _mm_add_pd(msum9 , _mm_mul_pd(mx14, mc149));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx0, mc150));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx1, mc151));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx2, mc152));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx3, mc153));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx4, mc154));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx5, mc155));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx6, mc156));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx7, mc157));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx8, mc158));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx9, mc159));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx10, mc160));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx11, mc161));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx12, mc162));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx13, mc163));\
                     msum10 = _mm_add_pd(msum10 , _mm_mul_pd(mx14, mc164));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx0, mc165));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx1, mc166));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx2, mc167));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx3, mc168));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx4, mc169));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx5, mc170));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx6, mc171));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx7, mc172));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx8, mc173));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx9, mc174));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx10, mc175));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx11, mc176));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx12, mc177));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx13, mc178));\
                     msum11 = _mm_add_pd(msum11 , _mm_mul_pd(mx14, mc179));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx0, mc180));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx1, mc181));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx2, mc182));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx3, mc183));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx4, mc184));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx5, mc185));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx6, mc186));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx7, mc187));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx8, mc188));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx9, mc189));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx10, mc190));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx11, mc191));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx12, mc192));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx13, mc193));\
                     msum12 = _mm_add_pd(msum12 , _mm_mul_pd(mx14, mc194));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx0, mc195));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx1, mc196));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx2, mc197));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx3, mc198));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx4, mc199));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx5, mc200));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx6, mc201));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx7, mc202));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx8, mc203));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx9, mc204));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx10, mc205));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx11, mc206));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx12, mc207));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx13, mc208));\
                     msum13 = _mm_add_pd(msum13 , _mm_mul_pd(mx14, mc209));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx0, mc210));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx1, mc211));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx2, mc212));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx3, mc213));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx4, mc214));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx5, mc215));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx6, mc216));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx7, mc217));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx8, mc218));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx9, mc219));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx10, mc220));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx11, mc221));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx12, mc222));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx13, mc223));\
                     msum14 = _mm_add_pd(msum14 , _mm_mul_pd(mx14, mc224));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx0, mc225));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx1, mc226));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx2, mc227));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx3, mc228));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx4, mc229));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx5, mc230));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx6, mc231));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx7, mc232));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx8, mc233));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx9, mc234));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx10, mc235));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx11, mc236));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx12, mc237));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx13, mc238));\
                     msum15 = _mm_add_pd(msum15 , _mm_mul_pd(mx14, mc239));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx0, mc240));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx1, mc241));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx2, mc242));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx3, mc243));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx4, mc244));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx5, mc245));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx6, mc246));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx7, mc247));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx8, mc248));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx9, mc249));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx10, mc250));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx11, mc251));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx12, mc252));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx13, mc253));\
                     msum16 = _mm_add_pd(msum16 , _mm_mul_pd(mx14, mc254));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx0, mc255));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx1, mc256));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx2, mc257));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx3, mc258));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx4, mc259));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx5, mc260));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx6, mc261));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx7, mc262));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx8, mc263));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx9, mc264));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx10, mc265));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx11, mc266));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx12, mc267));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx13, mc268));\
                     msum17 = _mm_add_pd(msum17 , _mm_mul_pd(mx14, mc269));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx0, mc270));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx1, mc271));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx2, mc272));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx3, mc273));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx4, mc274));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx5, mc275));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx6, mc276));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx7, mc277));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx8, mc278));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx9, mc279));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx10, mc280));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx11, mc281));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx12, mc282));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx13, mc283));\
                     msum18 = _mm_add_pd(msum18 , _mm_mul_pd(mx14, mc284));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx0, mc285));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx1, mc286));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx2, mc287));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx3, mc288));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx4, mc289));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx5, mc290));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx6, mc291));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx7, mc292));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx8, mc293));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx9, mc294));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx10, mc295));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx11, mc296));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx12, mc297));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx13, mc298));\
                     msum19 = _mm_add_pd(msum19 , _mm_mul_pd(mx14, mc299));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx0, mc300));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx1, mc301));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx2, mc302));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx3, mc303));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx4, mc304));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx5, mc305));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx6, mc306));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx7, mc307));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx8, mc308));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx9, mc309));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx10, mc310));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx11, mc311));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx12, mc312));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx13, mc313));\
                     msum20 = _mm_add_pd(msum20 , _mm_mul_pd(mx14, mc314));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx0, mc315));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx1, mc316));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx2, mc317));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx3, mc318));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx4, mc319));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx5, mc320));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx6, mc321));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx7, mc322));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx8, mc323));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx9, mc324));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx10, mc325));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx11, mc326));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx12, mc327));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx13, mc328));\
                     msum21 = _mm_add_pd(msum21 , _mm_mul_pd(mx14, mc329));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx0, mc330));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx1, mc331));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx2, mc332));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx3, mc333));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx4, mc334));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx5, mc335));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx6, mc336));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx7, mc337));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx8, mc338));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx9, mc339));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx10, mc340));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx11, mc341));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx12, mc342));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx13, mc343));\
                     msum22 = _mm_add_pd(msum22 , _mm_mul_pd(mx14, mc344));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx0, mc345));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx1, mc346));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx2, mc347));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx3, mc348));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx4, mc349));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx5, mc350));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx6, mc351));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx7, mc352));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx8, mc353));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx9, mc354));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx10, mc355));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx11, mc356));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx12, mc357));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx13, mc358));\
                     msum23 = _mm_add_pd(msum23 , _mm_mul_pd(mx14, mc359));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx0, mc360));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx1, mc361));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx2, mc362));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx3, mc363));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx4, mc364));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx5, mc365));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx6, mc366));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx7, mc367));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx8, mc368));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx9, mc369));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx10, mc370));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx11, mc371));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx12, mc372));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx13, mc373));\
                     msum24 = _mm_add_pd(msum24 , _mm_mul_pd(mx14, mc374));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx0, mc375));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx1, mc376));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx2, mc377));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx3, mc378));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx4, mc379));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx5, mc380));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx6, mc381));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx7, mc382));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx8, mc383));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx9, mc384));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx10, mc385));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx11, mc386));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx12, mc387));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx13, mc388));\
                     msum25 = _mm_add_pd(msum25 , _mm_mul_pd(mx14, mc389));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx0, mc390));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx1, mc391));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx2, mc392));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx3, mc393));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx4, mc394));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx5, mc395));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx6, mc396));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx7, mc397));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx8, mc398));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx9, mc399));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx10, mc400));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx11, mc401));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx12, mc402));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx13, mc403));\
                     msum26 = _mm_add_pd(msum26 , _mm_mul_pd(mx14, mc404));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx0, mc405));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx1, mc406));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx2, mc407));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx3, mc408));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx4, mc409));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx5, mc410));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx6, mc411));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx7, mc412));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx8, mc413));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx9, mc414));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx10, mc415));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx11, mc416));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx12, mc417));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx13, mc418));\
                     msum27 = _mm_add_pd(msum27 , _mm_mul_pd(mx14, mc419));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx0, mc420));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx1, mc421));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx2, mc422));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx3, mc423));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx4, mc424));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx5, mc425));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx6, mc426));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx7, mc427));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx8, mc428));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx9, mc429));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx10, mc430));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx11, mc431));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx12, mc432));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx13, mc433));\
                     msum28 = _mm_add_pd(msum28 , _mm_mul_pd(mx14, mc434));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx0, mc435));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx1, mc436));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx2, mc437));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx3, mc438));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx4, mc439));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx5, mc440));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx6, mc441));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx7, mc442));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx8, mc443));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx9, mc444));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx10, mc445));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx11, mc446));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx12, mc447));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx13, mc448));\
                     msum29 = _mm_add_pd(msum29 , _mm_mul_pd(mx14, mc449))

#define setup_30(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0);\
                         msum6 = _mm_set_pd(0,0);\
                         msum7 = _mm_set_pd(0,0);\
                         msum8 = _mm_set_pd(0,0);\
                         msum9 = _mm_set_pd(0,0);\
                         msum10 = _mm_set_pd(0,0);\
                         msum11 = _mm_set_pd(0,0);\
                         msum12 = _mm_set_pd(0,0);\
                         msum13 = _mm_set_pd(0,0);\
                         msum14 = _mm_set_pd(0,0);\
                         msum15 = _mm_set_pd(0,0);\
                         msum16 = _mm_set_pd(0,0);\
                         msum17 = _mm_set_pd(0,0);\
                         msum18 = _mm_set_pd(0,0);\
                         msum19 = _mm_set_pd(0,0);\
                         msum20 = _mm_set_pd(0,0);\
                         msum21 = _mm_set_pd(0,0);\
                         msum22 = _mm_set_pd(0,0);\
                         msum23 = _mm_set_pd(0,0);\
                         msum24 = _mm_set_pd(0,0);\
                         msum25 = _mm_set_pd(0,0);\
                         msum26 = _mm_set_pd(0,0);\
                         msum27 = _mm_set_pd(0,0);\
                         msum28 = _mm_set_pd(0,0);\
                         msum29 = _mm_set_pd(0,0)

#define save_30() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  msum6 = _mm_hadd_pd(msum6, msum7);\
                  msum8 = _mm_hadd_pd(msum8, msum9);\
                  msum10 = _mm_hadd_pd(msum10, msum11);\
                  msum12 = _mm_hadd_pd(msum12, msum13);\
                  msum14 = _mm_hadd_pd(msum14, msum15);\
                  msum16 = _mm_hadd_pd(msum16, msum17);\
                  msum18 = _mm_hadd_pd(msum18, msum19);\
                  msum20 = _mm_hadd_pd(msum20, msum21);\
                  msum22 = _mm_hadd_pd(msum22, msum23);\
                  msum24 = _mm_hadd_pd(msum24, msum25);\
                  msum26 = _mm_hadd_pd(msum26, msum27);\
                  msum28 = _mm_hadd_pd(msum28, msum29);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4);\
                  _mm_storeu_pd(y + t1 + 6,msum6);\
                  _mm_storeu_pd(y + t1 + 8,msum8);\
                  _mm_storeu_pd(y + t1 + 10,msum10);\
                  _mm_storeu_pd(y + t1 + 12,msum12);\
                  _mm_storeu_pd(y + t1 + 14,msum14);\
                  _mm_storeu_pd(y + t1 + 16,msum16);\
                  _mm_storeu_pd(y + t1 + 18,msum18);\
                  _mm_storeu_pd(y + t1 + 20,msum20);\
                  _mm_storeu_pd(y + t1 + 22,msum22);\
                  _mm_storeu_pd(y + t1 + 24,msum24);\
                  _mm_storeu_pd(y + t1 + 26,msum26);\
                  _mm_storeu_pd(y + t1 + 28,msum28)

PetscErrorCode BSG_MatMult_30(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset){
    PetscInt k, k1, it, l, t1, t2;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

     __m128d mx0, mx1, mx2, mx3, mx4, mx5, mx6, mx7, mx8, mx9, mx10, mx11, mx12, mx13, mx14, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, mc18, mc19, mc20, mc21, mc22, mc23, mc24, mc25, mc26, mc27, mc28, mc29, mc30, mc31, mc32, mc33, mc34, mc35, mc36, mc37, mc38, mc39, mc40, mc41, mc42, mc43, mc44, mc45, mc46, mc47, mc48, mc49, mc50, mc51, mc52, mc53, mc54, mc55, mc56, mc57, mc58, mc59, mc60, mc61, mc62, mc63, mc64, mc65, mc66, mc67, mc68, mc69, mc70, mc71, mc72, mc73, mc74, mc75, mc76, mc77, mc78, mc79, mc80, mc81, mc82, mc83, mc84, mc85, mc86, mc87, mc88, mc89, mc90, mc91, mc92, mc93, mc94, mc95, mc96, mc97, mc98, mc99, mc100, mc101, mc102, mc103, mc104, mc105, mc106, mc107, mc108, mc109, mc110, mc111, mc112, mc113, mc114, mc115, mc116, mc117, mc118, mc119, mc120, mc121, mc122, mc123, mc124, mc125, mc126, mc127, mc128, mc129, mc130, mc131, mc132, mc133, mc134, mc135, mc136, mc137, mc138, mc139, mc140, mc141, mc142, mc143, mc144, mc145, mc146, mc147, mc148, mc149, mc150, mc151, mc152, mc153, mc154, mc155, mc156, mc157, mc158, mc159, mc160, mc161, mc162, mc163, mc164, mc165, mc166, mc167, mc168, mc169, mc170, mc171, mc172, mc173, mc174, mc175, mc176, mc177, mc178, mc179, mc180, mc181, mc182, mc183, mc184, mc185, mc186, mc187, mc188, mc189, mc190, mc191, mc192, mc193, mc194, mc195, mc196, mc197, mc198, mc199, mc200, mc201, mc202, mc203, mc204, mc205, mc206, mc207, mc208, mc209, mc210, mc211, mc212, mc213, mc214, mc215, mc216, mc217, mc218, mc219, mc220, mc221, mc222, mc223, mc224, mc225, mc226, mc227, mc228, mc229, mc230, mc231, mc232, mc233, mc234, mc235, mc236, mc237, mc238, mc239, mc240, mc241, mc242, mc243, mc244, mc245, mc246, mc247, mc248, mc249, mc250, mc251, mc252, mc253, mc254, mc255, mc256, mc257, mc258, mc259, mc260, mc261, mc262, mc263, mc264, mc265, mc266, mc267, mc268, mc269, mc270, mc271, mc272, mc273, mc274, mc275, mc276, mc277, mc278, mc279, mc280, mc281, mc282, mc283, mc284, mc285, mc286, mc287, mc288, mc289, mc290, mc291, mc292, mc293, mc294, mc295, mc296, mc297, mc298, mc299, mc300, mc301, mc302, mc303, mc304, mc305, mc306, mc307, mc308, mc309, mc310, mc311, mc312, mc313, mc314, mc315, mc316, mc317, mc318, mc319, mc320, mc321, mc322, mc323, mc324, mc325, mc326, mc327, mc328, mc329, mc330, mc331, mc332, mc333, mc334, mc335, mc336, mc337, mc338, mc339, mc340, mc341, mc342, mc343, mc344, mc345, mc346, mc347, mc348, mc349, mc350, mc351, mc352, mc353, mc354, mc355, mc356, mc357, mc358, mc359, mc360, mc361, mc362, mc363, mc364, mc365, mc366, mc367, mc368, mc369, mc370, mc371, mc372, mc373, mc374, mc375, mc376, mc377, mc378, mc379, mc380, mc381, mc382, mc383, mc384, mc385, mc386, mc387, mc388, mc389, mc390, mc391, mc392, mc393, mc394, mc395, mc396, mc397, mc398, mc399, mc400, mc401, mc402, mc403, mc404, mc405, mc406, mc407, mc408, mc409, mc410, mc411, mc412, mc413, mc414, mc415, mc416, mc417, mc418, mc419, mc420, mc421, mc422, mc423, mc424, mc425, mc426, mc427, mc428, mc429, mc430, mc431, mc432, mc433, mc434, mc435, mc436, mc437, mc438, mc439, mc440, mc441, mc442, mc443, mc444, mc445, mc446, mc447, mc448, mc449, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9, msum10, msum11, msum12, msum13, msum14, msum15, msum16, msum17, msum18, msum19, msum20, msum21, msum22, msum23, msum24, msum25, msum26, msum27, msum28, msum29;

    const PetscScalar *xt0, *ct0, *xt1, *ct1, *xt2, *ct2, *xt3, *ct3, *xt4, *ct4, *xt5, *ct5, *xt6, *ct6;
    xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2) * dof;
    xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2) * dof;
    xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2) * dof;
    xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2) * dof;
    xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2) * dof;
    xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2) * dof;
    xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2) * dof;

    for(k1 = (0) , it = 0; k1 < (1); k1+= l3threshold, it++){
        setupct(0, it);
        endval = min(1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_30(0);
            inline_30(3);
            inline_30(4);
            inline_30(5);
            inline_30(6);
            save_30();
        }
    }

    for(k1 = (1) , it = 0; k1 < (lda3); k1+= l3threshold, it++){
        setupct(1, it);
        endval = min(lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_30(1);
            inline_30(2);
            inline_30(3);
            inline_30(4);
            inline_30(5);
            inline_30(6);
            save_30();
        }
    }

    for(k1 = (lda3) , it = 0; k1 < (lda2); k1+= l3threshold, it++){
        setupct(2, it);
        endval = min(lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_30(lda3);
            inline_30(1);
            inline_30(2);
            inline_30(3);
            inline_30(4);
            inline_30(5);
            inline_30(6);
            save_30();
        }
    }

    for(k1 = (lda2) , it = 0; k1 < (lda1 - lda2); k1+= l3threshold, it++){
        setupct(3, it);
        endval = min(lda1 - lda2,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_30(lda2);
            inline_30(0);
            inline_30(1);
            inline_30(2);
            inline_30(3);
            inline_30(4);
            inline_30(5);
            inline_30(6);
            save_30();
        }
    }

    for(k1 = (lda1 - lda2) , it = 0; k1 < (lda1 - lda3); k1+= l3threshold, it++){
        setupct(4, it);
        endval = min(lda1 - lda3,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_30(lda1 - lda2);
            inline_30(0);
            inline_30(1);
            inline_30(2);
            inline_30(3);
            inline_30(4);
            inline_30(5);
            save_30();
        }
    }

    for(k1 = (lda1 - lda3) , it = 0; k1 < (lda1 - 1); k1+= l3threshold, it++){
        setupct(5, it);
        endval = min(lda1 - 1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_30(lda1 - lda3);
            inline_30(0);
            inline_30(1);
            inline_30(2);
            inline_30(3);
            inline_30(4);
            save_30();
        }
    }

    for(k1 = (lda1 - 1) , it = 0; k1 < (lda1); k1+= l3threshold, it++){
        setupct(6, it);
        endval = min(lda1,k1+l3threshold);
        for(k = k1; k < endval; k++){
            setup_30(lda1 - 1);
            inline_30(0);
            inline_30(1);
            inline_30(2);
            inline_30(3);
            save_30();
        }
    }

PetscFunctionReturn(0);
}


