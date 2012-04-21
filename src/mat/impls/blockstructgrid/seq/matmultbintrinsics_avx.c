#include <string.h>

//#define SPREFETCH
#define min(a,b) (a)<(b) ? (a) : (b)

//#define _mm_maskstore_pd(a,b,c) _mm_storeu_pd(a,c)

#include <omp.h>
extern int OPENMPB;


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



#ifdef _VEC4

#define ninline_4() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x01);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 4);\
                    mc2 = _mm256_loadu_pd(ct[l] + t2 + 8);\
                    mc3 = _mm256_loadu_pd(ct[l] + t2 + 12);\
                    msum0 = _mm256_add_pd(msum0 , _mm256_mul_pd(mx0, mc0));\
                    msum1 = _mm256_add_pd(msum1 , _mm256_mul_pd(mx1, mc1));\
                    msum2 = _mm256_add_pd(msum2 , _mm256_mul_pd(mx0, mc2));\
                    msum3 = _mm256_add_pd(msum3 , _mm256_mul_pd(mx1, mc3));\
                    }

#define nsetup_4(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm256_set_pd(0,0,0,0);\
                         msum1 = _mm256_set_pd(0,0,0,0);\
                         msum2 = _mm256_set_pd(0,0,0,0);\
                         msum3 = _mm256_set_pd(0,0,0,0)

#define nsave_4() \
                  msum0 = _mm256_hadd_pd(msum0, msum1);\
                  msum2 = _mm256_hadd_pd(msum2, msum3);\
                  msum0 = _mm256_hadd_pd(msum0, msum2);\
                  _mm256_storeu_pd(y + t1 + 0,msum0)

PetscErrorCode BSG_MatMult_4(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m256d mx0, mx1, mc0, mc1, mc2, mc3, msum0, msum1, msum2, msum3;

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(l,t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
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
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x01);\
                    mx2 = _mm256_loadu_pd(xt[l] + t1 + 4);\
                    mx2 = _mm256_permute2f128_pd(mx2, mx2, 0x00);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 4);\
                    mc2 = _mm256_loadu_pd(ct[l] + t2 + 8);\
                    mc3 = _mm256_loadu_pd(ct[l] + t2 + 12);\
                    mc4 = _mm256_loadu_pd(ct[l] + t2 + 16);\
                    mc5 = _mm256_loadu_pd(ct[l] + t2 + 20);\
                    msum0 = _mm256_add_pd(msum0 , _mm256_mul_pd(mx0, mc0));\
                    msum1 = _mm256_add_pd(msum1 , _mm256_mul_pd(mx1, mc1));\
                    msum2 = _mm256_add_pd(msum2 , _mm256_mul_pd(mx0, mc2));\
                    msum3 = _mm256_add_pd(msum3 , _mm256_mul_pd(mx1, mc3));\
                    msum4 = _mm256_add_pd(msum4 , _mm256_mul_pd(mx2, mc4));\
                    msum5 = _mm256_add_pd(msum5 , _mm256_mul_pd(mx2, mc5));\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 24);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 26);\
                    mc2 = _mm256_loadu_pd(ct[l] + t2 + 28);\
                    mc3 = _mm256_loadu_pd(ct[l] + t2 + 30);\
                    mc4 = _mm256_loadu_pd(ct[l] + t2 + 32);\
                    mc5 = _mm256_loadu_pd(ct[l] + t2 + 34);\
                    msum6 = _mm256_add_pd(msum6 , _mm256_mul_pd(mx0, mc0));\
                    msum7 = _mm256_add_pd(msum7 , _mm256_mul_pd(mx1, mc1));\
                    msum8 = _mm256_add_pd(msum8 , _mm256_mul_pd(mx0, mc2));\
                    msum9 = _mm256_add_pd(msum9 , _mm256_mul_pd(mx1, mc3));\
                    msum10 = _mm256_add_pd(msum10 , _mm256_mul_pd(mx2, mc4));\
                    msum11 = _mm256_add_pd(msum11 , _mm256_mul_pd(mx2, mc5));\
                    }

#define nsetup_6(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm256_set_pd(0,0,0,0);\
                         msum1 = _mm256_set_pd(0,0,0,0);\
                         msum2 = _mm256_set_pd(0,0,0,0);\
                         msum3 = _mm256_set_pd(0,0,0,0);\
                         msum4 = _mm256_set_pd(0,0,0,0);\
                         msum5 = _mm256_set_pd(0,0,0,0);\
                         msum6 = _mm256_set_pd(0,0,0,0);\
                         msum7 = _mm256_set_pd(0,0,0,0);\
                         msum8 = _mm256_set_pd(0,0,0,0);\
                         msum9 = _mm256_set_pd(0,0,0,0);\
                         msum10 = _mm256_set_pd(0,0,0,0);\
                         msum11 = _mm256_set_pd(0,0,0,0)

#define nsave_6() \
                  msum0 = _mm256_hadd_pd(msum0, msum1);\
                  msum2 = _mm256_hadd_pd(msum2, msum3);\
                  msum0 = _mm256_hadd_pd(msum0, msum2);\
                  msum4 = _mm256_hadd_pd(msum4, msum5);\
                  msum0 = _mm256_add_pd(msum0, msum4);\
                  msum6 = _mm256_hadd_pd(msum6, msum7);\
                  msum8 = _mm256_hadd_pd(msum8, msum9);\
                  msum6 = _mm256_hadd_pd(msum6, msum8);\
                  msum10 = _mm256_hadd_pd(msum10, msum11);\
                  msum6 = _mm256_add_pd(msum6, msum10);\
                  _mm256_storeu_pd(y + t1 + 0,msum0);\
                  _mm256_maskstore_pd(y + t1 + 4,xtemp,msum6)

PetscErrorCode BSG_MatMult_6(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m256d mx0, mx1, mx2, mc0, mc1, mc2, mc3, mc4, mc5, msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7, msum8, msum9 , msum10, msum11;
    __m256i xtemp = _mm256_set_epi32(0,0,0,0,-1,-1,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y, xtemp) private(l,t1,t2,mx0, mx1, mx2, msum0, msum1, msum2, msum3, msum4, msum5, mc0, mc1, mc2, mc3, mc4, mc5)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_6(rstart[k1]);
            ninline_6();
            nsave_6();
        }
    }

PetscFunctionReturn(0);
}

#define nsetup_Neven()	for(i=0;i<dofby4;i++){\
				 msum[i] = _mm256_set_pd(0,0,0,0);\
			 }\

#define nsave_Neven() for(i=0, j = 0;(j+4)<=dof;i++, j+=4){\
				_mm256_storeu_pd(y+k*dof+4*i,msum[i]);\
			 }\
			for(; (j+2) <=dof; j+=2, i++)\
				_mm256_maskstore_pd(y + k*dof + 4*i,xtemp,msum[i]);\

#define nsetup12_Neven(xt)	t1= k*dof+l2;\
			mx0 = _mm256_loadu_pd(xt+t1);\
                    	mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x01)

#define nsetup34_Neven(xt)	t1= k*dof+l2;\
			mx0 = _mm256_loadu_pd(xt+t1);\
                    	mx0 = _mm256_permute2f128_pd(mx0, mx0, 0x00)

#define ninline_stage1_Neven(ct,offset)	t2 = (k-(offset))*bs+l1*dof+l2*4;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        mc1 = _mm256_loadu_pd(ct+t2+4);\
                        mc2 = _mm256_loadu_pd(ct+t2+8);\
                        mc3 = _mm256_loadu_pd(ct+t2+12);\
                        mc0 = _mm256_mul_pd(mx0,mc0);\
                        mc1 = _mm256_mul_pd(mx1,mc1);\
                        mc2 = _mm256_mul_pd(mx0,mc2);\
                        mc3 = _mm256_mul_pd(mx1,mc3);\
                  	mc0 = _mm256_hadd_pd(mc0, mc1);\
                  	mc2 = _mm256_hadd_pd(mc2, mc3);\
			msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc2));\


#define ninline_stage3_Neven(ct, offset)	t2 = (k-(offset))*bs+l1*dof+l2*4;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        mc1 = _mm256_loadu_pd(ct+t2+4);\
                        mc0 = _mm256_mul_pd(mx0,mc0);\
                        mc1 = _mm256_mul_pd(mx0,mc1);\
			msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc1));\

#define ninline_stage2_Neven(ct, offset)	t2 = (k-(offset))*bs+l1*dof+l2*2;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        mc1 = _mm256_loadu_pd(ct+t2+2);\
                        mc2 = _mm256_loadu_pd(ct+t2+4);\
                        mc3 = _mm256_loadu_pd(ct+t2+6);\
                        mc0 = _mm256_mul_pd(mx0,mc0);\
                        mc1 = _mm256_mul_pd(mx1,mc1);\
                        mc2 = _mm256_mul_pd(mx0,mc2);\
                        mc3 = _mm256_mul_pd(mx1,mc3);\
                  	mc0 = _mm256_hadd_pd(mc0, mc1);\
                  	mc2 = _mm256_hadd_pd(mc2, mc3);\
			msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc2));\

#define ninline_stage4_Neven(ct, offset)	t2 = (k-(offset))*bs+l1*dof+l2*2;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        mc1 = _mm256_loadu_pd(ct+t2+2);\
                        mc0 = _mm256_mul_pd(mx0,mc0);\
                        mc1 = _mm256_mul_pd(mx0,mc1);\
			msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc1));\

#define ninline_Neven(offset) \
                    for(l3 = lbeg[k1]; l3 < lend[k1]; l3++){ \
			for(l2 = 0; (l2+4) <= dof; l2 += 4){\
				nsetup12_Neven(xt[l3]);\
				for(l1=0, l4=0; (l1+4)<=dof; l1+=4,l4++){\
					ninline_stage1_Neven(ct[l3], offset);\
				}\
				for(;(l1+2) <= dof; l1+=2,l4++){\
					ninline_stage2_Neven(ct[l3], offset);\
				}\
			}\
			for(; l2 < dof; l2 += 2){\
				nsetup34_Neven(xt[l3]);\
				for(l1=0, l4=0; (l1+4) <= dof; l1+=4, l4++){\
					ninline_stage3_Neven(ct[l3], offset);\
				}\
				for(;(l1+2) <= dof; l1+=2,l4++){\
					ninline_stage4_Neven(ct[l3], offset);\
				}\
			}\
			}

PetscErrorCode BSG_MatMult_Neven(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
	PetscInt i,j,k,l,k1,it, t1, t2, l1,l2,l3, l4;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby4 = (dof+3/4);
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	__m256d mx0, mx1, msum[dofby4], mc0, mc1, mc2, mc3;
    	__m256i xtemp = _mm256_set_epi32(0,0,0,0,-1,-1,-1,-1);
    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y,xtemp)private(t1,t2,l1,l2,l3,l4, mx0, mx1, msum, mc0, mc1, mc2, mc3, i)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_Neven();
            ninline_Neven(rstart[k1]);
            nsave_Neven();
        }
    }

PetscFunctionReturn(0);
}



#endif

