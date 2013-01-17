#include <string.h>

#define min(a,b) (a)<(b) ? (a) : (b)


#include <omp.h>
int OPENMPB = 0;


#include "../src/mat/impls/blockstructgrid/seq/matblockstructgrid.h"
#include "../src/mat/impls/blockstructgrid/seq/commonfunctions.h"



#define setupct(pos, iter)	count = stpoffset[pos] + iter*nos;\
				ct0 = ctl[count + 0];\
				ct1 = ctl[count + 1];\
				ct2 = ctl[count + 2];\
				ct3 = ctl[count + 3];\
				ct4 = ctl[count + 4];\
				ct5 = ctl[count + 5];\
				ct6 = ctl[count + 6]

#define setup_ct(rno) 	count = stpoffset[rno];\
			for(l = lbeg[rno]; l< lend[rno]; l++){\
				ct[l] = ctl[count+l];\
			}

#ifdef _VEC2




#define ninline_2() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm_loadu_pd(xt[l] + t1 + 0);\
                    mc0 = _mm_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm_loadu_pd(ct[l] + t2 + 2);\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc1));\
                    }

#define nsetup_2(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0)

#define nsave_2() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  _mm_storeu_pd(y + t1 + 0,msum0)

PetscErrorCode BSG_MatMult_2(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m128d mx0, mc0, mc1, msum0, msum1;

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(l,t1,t2, mx0, msum0, msum1, mc0, mc1)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_2(rstart[k1]);
            ninline_2();
            nsave_2();
        }
    }

    PetscFunctionReturn(0);
}



#define ninline_4() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm_loadu_pd(xt[l] + t1 + 2);\
                    mc0 = _mm_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm_loadu_pd(ct[l] + t2 + 2);\
                    mc2 = _mm_loadu_pd(ct[l] + t2 + 4);\
                    mc3 = _mm_loadu_pd(ct[l] + t2 + 6);\
                    mc4 = _mm_loadu_pd(ct[l] + t2 + 8);\
                    mc5 = _mm_loadu_pd(ct[l] + t2 + 10);\
                    mc6 = _mm_loadu_pd(ct[l] + t2 + 12);\
                    mc7 = _mm_loadu_pd(ct[l] + t2 + 14);\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc2));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc3));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc4));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc5));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc6));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc7));\
                    }

#define nsetup_4(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0)

#define nsave_4() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2)

PetscErrorCode BSG_MatMult_4(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m128d mx0, mx1, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, msum0, msum1, msum2, msum3;

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(l,t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_4(rstart[k1]);
            ninline_4();
            nsave_4();
        }
    }

    PetscFunctionReturn(0);
}




#define ninline_6() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm_loadu_pd(xt[l] + t1 + 2);\
                    mx2 = _mm_loadu_pd(xt[l] + t1 + 4);\
                    mc0 = _mm_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm_loadu_pd(ct[l] + t2 + 2);\
                    mc2 = _mm_loadu_pd(ct[l] + t2 + 4);\
                    mc3 = _mm_loadu_pd(ct[l] + t2 + 6);\
                    mc4 = _mm_loadu_pd(ct[l] + t2 + 8);\
                    mc5 = _mm_loadu_pd(ct[l] + t2 + 10);\
                    mc6 = _mm_loadu_pd(ct[l] + t2 + 12);\
                    mc7 = _mm_loadu_pd(ct[l] + t2 + 14);\
                    mc8 = _mm_loadu_pd(ct[l] + t2 + 16);\
                    mc9 = _mm_loadu_pd(ct[l] + t2 + 18);\
                    mc10 = _mm_loadu_pd(ct[l] + t2 + 20);\
                    mc11 = _mm_loadu_pd(ct[l] + t2 + 22);\
                    mc12 = _mm_loadu_pd(ct[l] + t2 + 24);\
                    mc13 = _mm_loadu_pd(ct[l] + t2 + 26);\
                    mc14 = _mm_loadu_pd(ct[l] + t2 + 28);\
                    mc15 = _mm_loadu_pd(ct[l] + t2 + 30);\
                    mc16 = _mm_loadu_pd(ct[l] + t2 + 32);\
                    mc17 = _mm_loadu_pd(ct[l] + t2 + 34);\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx0, mc0));\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx1, mc1));\
                    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mx2, mc2));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx0, mc3));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx1, mc4));\
                    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mx2, mc5));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx0, mc6));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx1, mc7));\
                    msum2 = _mm_add_pd(msum2 , _mm_mul_pd(mx2, mc8));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx0, mc9));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx1, mc10));\
                    msum3 = _mm_add_pd(msum3 , _mm_mul_pd(mx2, mc11));\
                    msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx0, mc12));\
                    msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx1, mc13));\
                    msum4 = _mm_add_pd(msum4 , _mm_mul_pd(mx2, mc14));\
                    msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx0, mc15));\
                    msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx1, mc16));\
                    msum5 = _mm_add_pd(msum5 , _mm_mul_pd(mx2, mc17));\
                    }

#define nsetup_6(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0);\
                         msum3 = _mm_set_pd(0,0);\
                         msum4 = _mm_set_pd(0,0);\
                         msum5 = _mm_set_pd(0,0)

#define nsave_6() \
                  msum0 = _mm_hadd_pd(msum0, msum1);\
                  msum2 = _mm_hadd_pd(msum2, msum3);\
                  msum4 = _mm_hadd_pd(msum4, msum5);\
                  _mm_storeu_pd(y + t1 + 0,msum0);\
                  _mm_storeu_pd(y + t1 + 2,msum2);\
                  _mm_storeu_pd(y + t1 + 4,msum4)

PetscErrorCode BSG_MatMult_6(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m128d mx0, mx1, mx2, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17, msum0, msum1, msum2, msum3, msum4, msum5;

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(t1,t2, mx0, mx1, mx2, msum0, msum1, msum2, msum3, msum4,msum5, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9,  mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_6(rstart[k1]);
            ninline_6();
            nsave_6();
        }
    }

PetscFunctionReturn(0);
}




#define nsetup_Neven()	for(i=0;i<dofby2;i++){\
				 msum[i] = _mm_set_pd(0,0);\
			 }\

#define nsave_Neven() for(i=0;i<dofby2;i++){\
				_mm_storeu_pd(y+k*dof+2*i,msum[i]);\
			 }\

#define nsetup12_Neven(xt)	t1= k*dof+l2;\
			mx0 = _mm_loadu_pd(xt+t1);\
                        mx1 = _mm_loadu_pd(xt+t1+2)

#define nsetup34_Neven(xt)	t1= k*dof+l2;\
			mx0 = _mm_loadu_pd(xt+t1)

#define ninline_stage1_Neven(ct,offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+2);\
                        mc4 = _mm_loadu_pd(ct+t2+2*dof);\
                        mc5 = _mm_loadu_pd(ct+t2+2*dof+2);\
                        mc6 = _mm_loadu_pd(ct+t2+3*dof);\
                        mc7 = _mm_loadu_pd(ct+t2+3*dof+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        mc4 = _mm_add_pd(_mm_mul_pd(mx0,mc4),_mm_mul_pd(mx1,mc5));\
                        mc6 = _mm_add_pd(_mm_mul_pd(mx0,mc6),_mm_mul_pd(mx1,mc7));\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
			msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))


#define ninline_stage3_Neven(ct, offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc4 = _mm_loadu_pd(ct+t2+2*dof);\
                        mc6 = _mm_loadu_pd(ct+t2+3*dof);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        mc4 = _mm_mul_pd(mx0,mc4);\
                        mc6 = _mm_mul_pd(mx0,mc6);\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
			msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))

#define ninline_stage2_Neven(ct, offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define ninline_stage4_Neven(ct, offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define ninline_Neven(offset) \
                    for(l3 = lbeg[k1]; l3 < lend[k1]; l3++){ \
			for(l2 = 0; l2 < dof-2; l2 += 4){\
				nsetup12_Neven(xt[l3]);\
				for(l1=0; l1<dofby2-1; l1+= 2){\
					ninline_stage1_Neven(ct[l3], offset);\
				}\
				for(;l1<dofby2;l1++){\
					ninline_stage2_Neven(ct[l3], offset);\
				}\
			}\
			for(; l2 < dof; l2 += 2){\
				nsetup34_Neven(xt[l3]);\
				for(l1=0; l1<dofby2-1; l1+= 2){\
					ninline_stage3_Neven(ct[l3], offset);\
				}\
				for(;l1<dofby2;l1++){\
					ninline_stage4_Neven(ct[l3], offset);\
				}\
			}\
			}

PetscErrorCode BSG_MatMult_Neven(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
	PetscInt i,k,l,k1,it, t1, t2, l1,l2,l3;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = dof/2;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	__m128d mx0, mx1, msum[dofby2], mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
//#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(t1,t2,l1,l2,l3, mx0, mx1, msum, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7,i)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_Neven();
            ninline_Neven(rstart[k1]);
            nsave_Neven();
        }
    }

PetscFunctionReturn(0);
}




 #define nsetup_1(offset) t1= k*dof; t2 = (k-(offset))*bs;\
  msum0 = 0.0
  
#define nsave_1() y[t1] = msum0

PetscInt BSG_MatMult_1(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt * lend, const PetscInt * rstart)
{
PetscInt k,l,k1, it, t1, t2, l1;
const PetscInt lda3 = m;
const PetscInt lda2 = lda3 * n;
const PetscInt lda1 = lda2 * p;
const PetscInt mnos = 3;
PetscScalar msum0;
PetscInt l3threshold = WORKINGSETSIZE / bs;
PetscInt count, endval;

const PetscScalar *xt[nos],*ct[nos];
for(k1 = 0; k1< nos; k1++)
{
        xt[k1] = x + ioff[k1] * dof;
}


for(k1 = 0; k1 < nregion; k1++){
setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(t1,t2,l, msum0)
for(k = rstart[k1] ; k<rstart[k1+1]; k++){
 nsetup_1(rstart[k1]);
 ninline_1();
 nsave_1();
}
}



#define ninline_3() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm_load1_pd(xt[l] + t1 + 2);\
                    mc0 = _mm_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm_loadu_pd(ct[l] + t2 + 2);\
                    mc2 = _mm_loadu_pd(ct[l] + t2 + 4);\
                    mc0 = _mm_mul_pd(mx0, mc0);\
                    mc1 = _mm_mul_pd(mx0, mc1);\
                    mc2 = _mm_mul_pd(mx0, mc2);\
                    msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc1));\
                    msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc2, mc2));\
                    mc0 = _mm_loadu_pd(ct[l] + t2 + 6);\
                    mc1 = _mm_loadu_pd(ct[l] + t2 + 8);\
                    msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx1, mc0));\
                    msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx1, mc1));\
                    }

#define nsetup_3(offset) \
                         t1 = k*3; t2 = (k-(offset))*9;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0)

#define nsave_3() \
                  _mm_storeu_pd(y + t1 + 0, msum0);\
                  _mm_maskstore_pd(y + t1 + 2, xtemp, msum1)

PetscErrorCode BSG_MatMult_3(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m128d mx0, mx1, mc0, mc1, mc2, msum0, msum1;

    __m128i xtemp = _mm_set_epi32(0,0,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y,xtemp) private(l,t1,t2,mx0, mx1, msum0, msum1, mc0, mc1, mc2)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_3(rstart[k1]);
            ninline_3();
            nsave_3();
        }
    }

PetscFunctionReturn(0);
}





#define nninline_5() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm_loadu_pd(xt[l] + t1 + 2);\
                    mx2 = _mm_load1_pd(xt[l] + t1 + 4);\
                    mc0 = _mm_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm_loadu_pd(ct[l] + t2 + 2);\
                    mc2 = _mm_loadu_pd(ct[l] + t2 + 4);\
                    mc3 = _mm_loadu_pd(ct[l] + t2 + 6);\
                    mc4 = _mm_loadu_pd(ct[l] + t2 + 8);\
                    mc5 = _mm_loadu_pd(ct[l] + t2 + 10);\
                    mc6 = _mm_loadu_pd(ct[l] + t2 + 12);\
                    mc7 = _mm_loadu_pd(ct[l] + t2 + 14);\
                    mc8 = _mm_loadu_pd(ct[l] + t2 + 16);\
                    mc9 = _mm_loadu_pd(ct[l] + t2 + 18);\
                    mc0 = _mm_mul_pd(mx0, mc0);\
                    mc0 = _mm_add_pd(mc0, _mm_mul_pd(mx1, mc1));\
                    mc2 = _mm_mul_pd(mx0, mc2);\
                    mc2 = _mm_add_pd(mc2, _mm_mul_pd(mx1, mc3));\
                    mc4 = _mm_mul_pd(mx0, mc4);\
                    mc4 = _mm_add_pd(mc4, _mm_mul_pd(mx1, mc5));\
                    mc6 = _mm_mul_pd(mx0, mc6);\
                    mc6 = _mm_add_pd(mc6, _mm_mul_pd(mx1, mc7));\
                    mc8 = _mm_mul_pd(mx0, mc8);\
                    mc8 = _mm_add_pd(mc8, _mm_mul_pd(mx1, mc9));\
                    msum0 = _mm_add_pd(msum0, _mm_hadd_pd(mc0, mc2));\
                    msum1 = _mm_add_pd(msum1, _mm_hadd_pd(mc4, mc6));\
                    msum2 = _mm_add_pd(msum2, _mm_hadd_pd(mc8, mc8));\
                    mc0 = _mm_loadu_pd(ct[l] + t2 + 20);\
                    mc1 = _mm_loadu_pd(ct[l] + t2 + 22);\
                    mc2 = _mm_loadu_pd(ct[l] + t2 + 24);\
                    msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx2, mc0));\
                    msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx2, mc1));\
                    msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx2, mc2));\
                    }

#define nnsetup_5(offset) \
                         t1 = k*5; t2 = (k-(offset))*25;\
                         msum0 = _mm_set_pd(0,0);\
                         msum1 = _mm_set_pd(0,0);\
                         msum2 = _mm_set_pd(0,0)

#define nnsave_5() \
                  _mm_storeu_pd(y + t1 + 0, msum0);\
                  _mm_storeu_pd(y + t1 + 2, msum1);\
                  _mm_maskstore_pd(y + t1 + 4, xtemp, msum2)


PetscErrorCode BSG_MatMult_5(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m128d mx0, mx1, mx2, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, msum0, msum1, msum2;

    __m128i xtemp = _mm_set_epi32(0,0,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y,xtemp) private(l,t1,t2, mx0, mx1,mx2, msum0, msum1, msum2, mc0, mc1, mc2, mc3, mc4,mc5, mc6, mc7, mc8, mc9)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nnsetup_5(rstart[k1]);
            nninline_5();
            nnsave_5();
        }
    }

PetscFunctionReturn(0);
}




#define nsetup_Nodd()   for(i=0;i<dofby2;i++){\
                                 msum[i] = _mm_set_pd(0,0);\
                         }\

#define nsave_Nodd() for(i=0;i<dofby2-1;i++){\
                                _mm_storeu_pd(y+k*dof+2*i,msum[i]);\
                         }\
  _mm_maskstore_pd(y+(k+1)*dof-1,xtemp,msum[dofby2-1])

#define nsetup123_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt+t1);\
                        mx1 = _mm_loadu_pd(xt+t1+2)

#define nsetup456_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt+t1)

#define nsetup789_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_load1_pd(xt+t1)

#define ninline_stage1_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+1);\
                        mc4 = _mm_loadu_pd(ct+t2+2*dof-2);\
                        mc5 = _mm_loadu_pd(ct+t2+2*dof);\
                        mc6 = _mm_loadu_pd(ct+t2+3*dof-3);\
                        mc7 = _mm_loadu_pd(ct+t2+3*dof-1);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        mc4 = _mm_add_pd(_mm_mul_pd(mx0,mc4),_mm_mul_pd(mx1,mc5));\
                        mc6 = _mm_add_pd(_mm_mul_pd(mx0,mc6),_mm_mul_pd(mx1,mc7));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))


#define ninline_stage4_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc4 = _mm_loadu_pd(ct+t2+2*dof-2);\
                        mc6 = _mm_loadu_pd(ct+t2+3*dof-3);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        mc4 = _mm_mul_pd(mx0,mc4);\
                        mc6 = _mm_mul_pd(mx0,mc6);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))

#define ninline_stage7_Nodd(ct, offset) t2 = (k-(offset))*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_mul_pd(mx0,mc0));\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_mul_pd(mx0,mc1))

#define ninline_stage2_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+1);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define ninline_stage5_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define ninline_stage8_Nodd(ct, offset) t2 = (k-(offset))*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_mul_pd(mx0,mc0))

#define ninline_stage3_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define ninline_stage6_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define ninline_stage9_Nodd(ct, offset) t2 = (k-(offset))*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_mul_pd(mx0,mc0));\

#define ninline_Nodd(offset)	\
                    for(l3 = lbeg[k1]; l3 < lend[k1]; l3++){ \
			for(l2 = 0; l2 < dof-3; l2 += 4){\
				nsetup123_Nodd(xt[l3]);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					ninline_stage1_Nodd(ct[l3], offset);\
				}\
				for(;l1<dofby2-1;l1++){\
					ninline_stage2_Nodd(ct[l3], offset);\
				}\
				for(;l1<dofby2;l1++){\
					ninline_stage3_Nodd(ct[l3], offset);\
				}\
			}\
			for(; l2 < dof-1; l2 += 2){\
				nsetup456_Nodd(xt[l3]);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					ninline_stage4_Nodd(ct[l3], offset);\
				}\
				for(;l1<dofby2-1;l1++){\
					ninline_stage5_Nodd(ct[l3], offset);\
				}\
				for(;l1<dofby2;l1++){\
					ninline_stage6_Nodd(ct[l3], offset);\
				}\
			}\
			for(; l2 < dof; l2 ++){\
				nsetup789_Nodd(xt[l3]);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					ninline_stage7_Nodd(ct[l3], offset);\
				}\
				for(;l1<dofby2-1;l1++){\
					ninline_stage8_Nodd(ct[l3], offset);\
				}\
				for(;l1<dofby2;l1++){\
					ninline_stage9_Nodd(ct[l3], offset);\
				}\
			}\
		}


PetscErrorCode BSG_MatMult_Nodd(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    
	const PetscInt l3threshold = WORKINGSETSIZE / bs;
   	 PetscInt count, endval;

	PetscInt i,k,l,k1, it,t1, t2, l1,l2, l3;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = (dof+1)/2;
	__m128d mx0, mx1, msum[dofby2], mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
    __m128i xtemp = _mm_set_epi32(0,0,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 <  nregion; k1++){
        setup_ct(k1);
//#pragma omp parallel for if(OPENMPB) shared(xt,ct,y,xtemp) private(t1,t2, i, l1, l2,l3,mx0, mx1, msum, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7 )
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
		nsetup_Nodd();
		ninline_Nodd(rstart[k1]);
		nsave_Nodd();
        }
    }

PetscFunctionReturn(0);
}
 #endif

