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

#define ninline_2() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mx0 = _mm256_permute2f128_pd(mx0, mx0, 0x00);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 2);\
                    msum0 = _mm256_add_pd(msum0 , _mm256_mul_pd(mx0, mc0));\
                    msum1 = _mm256_add_pd(msum1 , _mm256_mul_pd(mx0, mc1));\
                    }

#define nsetup_2(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm256_set_pd(0,0,0,0);\
                         msum1 = _mm256_set_pd(0,0,0,0);\
                         msum2 = _mm256_set_pd(0,0,0,0)

#define nsave_2() \
                  msum0 = _mm256_hadd_pd(msum0, msum1);\
                  _mm256_maskstore_pd(y + t1 + 0,xtemp,msum0)

PetscErrorCode BSG_MatMult_2(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m256d mx0, mx1, mc0, mc1, mc2, mc3, msum0, msum1, msum2, msum3;
    __m256i xtemp = _mm256_set_epi32(0,0,0,0,-1,-1,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(l,t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_2(rstart[k1]);
            ninline_2();
            nsave_2();
        }
    }

PetscFunctionReturn(0);
}

#define inline_4(l,m) \
                    mx0 = _mm256_loadu_pd(xt##l + t1 + 0);\
                    mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x01);\
                    mc0 = _mm256_loadu_pd(ct##l + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct##l + t2 + 4);\
                    mc2 = _mm256_loadu_pd(ct##l + t2 + 8);\
                    mc3 = _mm256_loadu_pd(ct##l + t2 + 12);\
                    msum0 = _mm256_add_pd(msum0 , _mm256_mul_pd(mx0, mc0));\
                    msum1 = _mm256_add_pd(msum1 , _mm256_mul_pd(mx1, mc1));\
                    msum2 = _mm256_add_pd(msum2 , _mm256_mul_pd(mx0, mc2));\
                    msum3 = _mm256_add_pd(msum3 , _mm256_mul_pd(mx1, mc3));\

#define setup_4(offset) \
                         t1 = k*dof; t2 = (k-(k1))*bs;\
                         msum0 = _mm256_set_pd(0,0,0,0);\
                         msum1 = _mm256_set_pd(0,0,0,0);\
                         msum2 = _mm256_set_pd(0,0,0,0);\
                         msum3 = _mm256_set_pd(0,0,0,0)

#define save_4() \
                  msum0 = _mm256_hadd_pd(msum0, msum1);\
                  msum2 = _mm256_hadd_pd(msum2, msum3);\
                  msum0 = _mm256_hadd_pd(msum0, msum2);\
                  _mm256_storeu_pd(y + t1 + 0,msum0)

PetscInt BSG_MatMult_4_11(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt k,l,k1, it, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	__m256d mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (ioff[0] )*dof;
		xt1 = x + (ioff[1] )*dof;
		xt2 = x + (ioff[2] )*dof;
		xt3 = x + (ioff[3] )*dof;
		xt4 = x + (ioff[4] )*dof;
		xt5 = x + (ioff[5] )*dof;
		xt6 = x + (ioff[6] )*dof;
	for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
#pragma omp parallel for if(OPENMPB) shared(xt0,xt1,xt2,xt3,xt4,xt5,xt6,ct0,ct1,ct2,ct3,ct4,ct5,ct6,xt7,ct7,y) private(t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
		for(k = k1; k < endval; k++)
		{
		setup_4(0);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
		}
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
#pragma omp parallel for if(OPENMPB) shared(xt0,xt1,xt2,xt3,xt4,xt5,xt6,ct0,ct1,ct2,ct3,ct4,ct5,ct6,xt7,ct7,y) private(t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
		for(k = k1; k < endval; k++)
		{
		setup_4(1);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
		}
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
#pragma omp parallel for if(OPENMPB) shared(xt0,xt1,xt2,xt3,xt4,xt5,xt6,ct0,ct1,ct2,ct3,ct4,ct5,ct6,xt7,ct7,y) private(t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
		for(k = k1; k < endval; k++)
		{
		setup_4(lda3);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
		}
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
#pragma omp parallel for if(OPENMPB) shared(xt0,xt1,xt2,xt3,xt4,xt5,xt6,ct0,ct1,ct2,ct3,ct4,ct5,ct6,xt7,ct7,y) private(t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
		for(k = k1; k < endval; k++)
		{
		setup_4(lda2);
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
		}
		}

		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		{
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
#pragma omp parallel for if(OPENMPB) shared(xt0,xt1,xt2,xt3,xt4,xt5,xt6,ct0,ct1,ct2,ct3,ct4,ct5,ct6,xt7,ct7,y) private(t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
		for(k = k1; k < endval; k++)
		{
		setup_4(lda1-lda2);
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		save_4();
		}
		}

		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		{
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
#pragma omp parallel for if(OPENMPB) shared(xt0,xt1,xt2,xt3,xt4,xt5,xt6,ct0,ct1,ct2,ct3,ct4,ct5,ct6,xt7,ct7,y) private(t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
		for(k = k1; k < endval; k++)
		{
		setup_4(lda1-lda3);
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		save_4();
		}
		}

		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		{
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
#pragma omp parallel for if(OPENMPB) shared(xt0,xt1,xt2,xt3,xt4,xt5,xt6,ct0,ct1,ct2,ct3,ct4,ct5,ct6,xt7,ct7,y) private(t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
		for(k = k1; k < endval; k++)
		{
		setup_4(lda1-1);
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		save_4();
		}
		}
	PetscFunctionReturn(0);
}

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
			for(; (l2+2) <= dof; l2 += 2){\
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
//#pragma omp parallel for if(OPENMPB) shared(xt,ct,y,xtemp)private(t1,t2,l1,l2,l3,l4, mx0, mx1, msum, mc0, mc1, mc2, mc3, i)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_Neven();
            ninline_Neven(rstart[k1]);
            nsave_Neven();
        }
    }

PetscFunctionReturn(0);
}

#define ninline_3() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mx0 = _mm256_permute2f128_pd(mx0, mx0, 0x00);\
                        mx1 = _mm256_loadu_pd(xt[l] + t1 + 2);\
			mx1 = _mm256_permute2f128_pd(mx1, mx1, 0x00);\
			mx1 = _mm256_permute_pd(mx1, 0x0);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 2);\
                    mc2 = _mm256_loadu_pd(ct[l] + t2 + 4);\
                    mc3 = _mm256_loadu_pd(ct[l] + t2 + 6);\
                    mc4 = _mm256_loadu_pd(ct[l] + t2 + 8);\
                    msum0 = _mm256_add_pd(msum0 , _mm256_mul_pd(mx0, mc0));\
                    msum1 = _mm256_add_pd(msum1 , _mm256_mul_pd(mx0, mc1));\
                    msum2 = _mm256_add_pd(msum2 , _mm256_mul_pd(mx1, mc2));\
                    msum3 = _mm256_add_pd(msum3 , _mm256_mul_pd(mx0, mc3));\
                    msum4 = _mm256_add_pd(msum4 , _mm256_mul_pd(mx1, mc4));\
                    }

#define nsetup_3(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm256_set_pd(0,0,0,0);\
                         msum1 = _mm256_set_pd(0,0,0,0);\
                         msum2 = _mm256_set_pd(0,0,0,0);\
                         msum3 = _mm256_set_pd(0,0,0,0);\
                         msum4 = _mm256_set_pd(0,0,0,0)

#define nsave_3() \
                  msum2 = _mm256_add_pd(msum2, _mm256_hadd_pd(msum0, msum1));\
                  msum4 = _mm256_add_pd(msum4, _mm256_hadd_pd(msum3, msum3));\
                  _mm256_maskstore_pd(y + t1 + 0,xtemp1,msum2);\
                  _mm256_maskstore_pd(y + t1 + 2,xtemp2,msum4)

PetscErrorCode BSG_MatMult_3(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m256d mx0, mx1, mc0, mc1, mc2, mc3, mc4, msum0, msum1, msum2, msum3, msum4;
    __m256i xtemp1 = _mm256_set_epi32(0,0,0,0,-1,-1,-1,-1);
    __m256i xtemp2 = _mm256_set_epi32(0,0,0,0,0,0,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(l,t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_3(rstart[k1]);
            ninline_3();
            nsave_3();
        }
    }

PetscFunctionReturn(0);
}

#define ninline_5() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x01);\
                        mx2 = _mm256_loadu_pd(xt[l] + t1 + 4);\
			mx2 = _mm256_permute2f128_pd(mx2, mx2, 0x00);\
			mx2 = _mm256_permute_pd(mx2, 0x0);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 4);\
                    mc2 = _mm256_loadu_pd(ct[l] + t2 + 8);\
                    mc3 = _mm256_loadu_pd(ct[l] + t2 + 12);\
                    mc4 = _mm256_loadu_pd(ct[l] + t2 + 16);\
                    mc5 = _mm256_loadu_pd(ct[l] + t2 + 20);\
                    mc7 = _mm256_loadu_pd(ct[l] + t2 + 24);\
                    msum0 = _mm256_add_pd(msum0 , _mm256_mul_pd(mx0, mc0));\
                    msum1 = _mm256_add_pd(msum1 , _mm256_mul_pd(mx1, mc1));\
                    msum2 = _mm256_add_pd(msum2 , _mm256_mul_pd(mx0, mc2));\
                    msum3 = _mm256_add_pd(msum3 , _mm256_mul_pd(mx1, mc3));\
                    msum4 = _mm256_add_pd(msum4 , _mm256_mul_pd(mx2, mc4));\
                    msum5 = _mm256_add_pd(msum5 , _mm256_mul_pd(mx0, mc5));\
                    msum7 = _mm256_add_pd(msum7 , _mm256_mul_pd(mx2, mc7));\
                    }

#define nsetup_5(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm256_set_pd(0,0,0,0);\
                         msum1 = _mm256_set_pd(0,0,0,0);\
                         msum2 = _mm256_set_pd(0,0,0,0);\
                         msum3 = _mm256_set_pd(0,0,0,0);\
                         msum4 = _mm256_set_pd(0,0,0,0);\
                         msum5 = _mm256_set_pd(0,0,0,0);\
                         msum6 = _mm256_set_pd(0,0,0,0);\
                         msum7 = _mm256_set_pd(0,0,0,0)

#define nsave_5() \
                  	msum0 = _mm256_hadd_pd(msum0,msum1);\
                  	msum2 = _mm256_hadd_pd(msum2,msum3);\
			msum0 = _mm256_hadd_pd(msum0,msum2);\
			msum4 = _mm256_add_pd(msum4,msum0);\
\
                    	msum6 = _mm256_permute2f128_pd(msum5, msum5, 0x01);\
			msum5 = _mm256_hadd_pd(msum5,msum6);\
			msum5 = _mm256_hadd_pd(msum5,msum5);\
                        msum7 = _mm256_add_pd(msum7,msum5);\
                  _mm256_storeu_pd(y + t1 + 0,msum4);\
                  _mm256_maskstore_pd(y + t1 + 4,xtemp2,msum7)

PetscErrorCode BSG_MatMult_5(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m256d mx0, mx1, mx2, mc0, mc1, mc2, mc3, mc4, mc5, mc7,  msum0, msum1, msum2, msum3, msum4, msum5, msum6, msum7;
    __m256i xtemp1 = _mm256_set_epi32(0,0,0,0,-1,-1,-1,-1);
    __m256i xtemp2 = _mm256_set_epi32(0,0,0,0,0,0,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(l,t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_5(rstart[k1]);
            ninline_5();
            nsave_5();
        }
    }

PetscFunctionReturn(0);
}

#define ninline_7() \
                    for(l = lbeg[k1]; l < lend[k1]; l++){ \
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x01);\
			mx2 = _mm256_loadu_pd(xt[l] + t1 +4);\
                    	mx2 = _mm256_permute2f128_pd(mx2, mx2, 0x00);\
                        mx3 = _mm256_loadu_pd(xt[l] + t1 + 6);\
			mx3 = _mm256_permute2f128_pd(mx3, mx3, 0x00);\
			mx3 = _mm256_permute_pd(mx3, 0x0);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 4);\
                    mc2 = _mm256_loadu_pd(ct[l] + t2 + 8);\
                    mc3 = _mm256_loadu_pd(ct[l] + t2 + 12);\
                    mc4 = _mm256_loadu_pd(ct[l] + t2 + 16);\
                    mc5 = _mm256_loadu_pd(ct[l] + t2 + 20);\
                    mc6 = _mm256_loadu_pd(ct[l] + t2 + 24);\
                    mc7 = _mm256_loadu_pd(ct[l] + t2 + 28);\
                    mc8 = _mm256_loadu_pd(ct[l] + t2 + 30);\
                    mc9 = _mm256_loadu_pd(ct[l] + t2 + 32);\
                    mc10 = _mm256_loadu_pd(ct[l] + t2 + 34);\
                    mc11 = _mm256_loadu_pd(ct[l] + t2 + 36);\
                    mc12 = _mm256_loadu_pd(ct[l] + t2 + 38);\
                    mc13 = _mm256_loadu_pd(ct[l] + t2 + 40);\
                    mc14 = _mm256_loadu_pd(ct[l] + t2 + 42);\
                    mc16 = _mm256_loadu_pd(ct[l] + t2 + 46);\
                    mc17 = _mm256_loadu_pd(ct[l] + t2 + 48);\
			mc0 = _mm256_mul_pd(mx0,mc0);\
                        mc1 = _mm256_mul_pd(mx1,mc1);\
                        mc2 = _mm256_mul_pd(mx0,mc2);\
                        mc3 = _mm256_mul_pd(mx1,mc3);\
                        mc4 = _mm256_mul_pd(mx2,mc4);\
                        mc5 = _mm256_mul_pd(mx2,mc5);\
			mc6 = _mm256_mul_pd(mx3,mc6);\
                        mc7 = _mm256_mul_pd(mx0,mc7);\
                        mc8 = _mm256_mul_pd(mx1,mc8);\
                        mc9 = _mm256_mul_pd(mx0,mc9);\
                        mc10 = _mm256_mul_pd(mx1,mc10);\
                        mc11 = _mm256_mul_pd(mx2,mc11);\
                        mc12 = _mm256_mul_pd(mx2,mc12);\
			mc13 = _mm256_mul_pd(mx3,mc13);\
                        mc14 = _mm256_mul_pd(mx0,mc14);\
			mc16 = _mm256_mul_pd(mx2,mc16);\
			mc17 = _mm256_mul_pd(mx3,mc17);\
                  	\
			mc0 = _mm256_hadd_pd(mc0,mc1);\
                  	mc2 = _mm256_hadd_pd(mc2,mc3);\
			mc0 = _mm256_hadd_pd(mc0,mc2);\
	 		mc4 = _mm256_hadd_pd(mc4,mc5);\
                  	mc7 = _mm256_hadd_pd(mc7,mc8);\
                  	mc9 = _mm256_hadd_pd(mc9,mc10);\
			mc7 = _mm256_hadd_pd(mc7,mc9);\
			mc11 = _mm256_hadd_pd(mc11,mc12);\
                    	mc15 = _mm256_permute2f128_pd(mc14, mc14, 0x01);\
			mc14 = _mm256_hadd_pd(mc14,mc15);\
			mc14 = _mm256_hadd_pd(mc14,mc15);\
			mc16 = _mm256_hadd_pd(mc16,mc16);\
			\
			mc0 = _mm256_add_pd(mc0,mc4);\
			mc0 = _mm256_add_pd(mc0,mc6);\
			mc7 = _mm256_add_pd(mc7,mc11);\
			mc7 = _mm256_add_pd(mc7,mc13);\
                        mc14 = _mm256_add_pd(mc14,mc16);\
                        mc14 = _mm256_add_pd(mc14,mc17);\
                        \
			msum0 = _mm256_add_pd(msum0,mc0);\
                        msum1 = _mm256_add_pd(msum1, mc7);\
                        msum2 = _mm256_add_pd(msum2,mc14);\
                    }

#define nsetup_7(offset) \
                         t1 = k*dof; t2 = (k-(offset))*bs;\
                         msum0 = _mm256_set_pd(0,0,0,0);\
                         msum1 = _mm256_set_pd(0,0,0,0);\
                         msum2 = _mm256_set_pd(0,0,0,0)

#define nsave_7() \
                  _mm256_storeu_pd(y + t1 + 0,msum0);\
                  _mm256_maskstore_pd(y + t1 + 4,xtemp1,msum1);\
                  _mm256_maskstore_pd(y + t1 + 6,xtemp2,msum2)

PetscErrorCode BSG_MatMult_7(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    PetscInt k, k1, it, l, t1, t2, l1;
    const PetscInt lda3 = m ;
    const PetscInt lda2 = lda3 * n;
    const PetscInt lda1 = lda2 * p;
    const PetscInt mnos = dim;
    const PetscInt l3threshold = WORKINGSETSIZE / bs;
    PetscInt count, endval;

    __m256d mx0, mx1, mx2, mx3, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9, mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17,  msum0, msum1, msum2, msum3;
    __m256i xtemp1 = _mm256_set_epi32(0,0,0,0,-1,-1,-1,-1);
    __m256i xtemp2 = _mm256_set_epi32(0,0,0,0,0,0,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 < nregion; k1++){
        setup_ct(k1);
#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(l,t1,t2,mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
            nsetup_7(rstart[k1]);
            ninline_7();
            nsave_7();
        }
    }

PetscFunctionReturn(0);
}

#define nsetup_Nodd()   for(i=0;i<dofby4;i++){\
                                 msum[i] = _mm256_set_pd(0,0,0,0);\
                         }\


#define nsave_Nodd() for(i=0, j =0;(j+4)<=dof;i++, j+=4){\
                                _mm256_storeu_pd(y+k*dof+j,msum[i]);\
                         }\
			for(; (j+2) <=dof; j+=2, i++)\
				_mm256_maskstore_pd(y + k*dof + j,xtemp1,msum[i]);\
  			_mm256_maskstore_pd(y+k*dof+j,xtemp2,msum[i])

#define nsetup127_Nodd(xt) t1= k*dof+l2;\
			mx0 = _mm256_loadu_pd(xt+t1);\
                    	mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x01)

#define nsetup348_Nodd(xt) t1= k*dof+l2;\
			mx0 = _mm256_loadu_pd(xt+t1);\
                    	mx0 = _mm256_permute2f128_pd(mx0, mx0, 0x00)

#define nsetup569_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm256_loadu_pd(xt+t1);\
			mx1 = _mm256_permute2f128_pd(mx0, mx0, 0x00);\
			mx0 = _mm256_permute_pd(mx1, 0x0)

#define ninline_stage1_Nodd(ct,offset)	t2 = (k-(offset))*bs+l1*dof+l2*4;\
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


#define ninline_stage3_Nodd(ct, offset)	t2 = (k-(offset))*bs+l1*dof+l2*4;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        mc1 = _mm256_loadu_pd(ct+t2+4);\
                        mc0 = _mm256_mul_pd(mx0,mc0);\
                        mc1 = _mm256_mul_pd(mx0,mc1);\
			msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc1));\

#define ninline_stage5_Nodd(ct, offset)   t2 = (k-(offset))*bs+l1*dof+l2*4;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        msum[l4] = _mm256_add_pd(msum[l4], _mm256_mul_pd(mx0,mc0))

#define ninline_stage2_Nodd(ct, offset)	t2 = (k-(offset))*bs+l1*dof+l2*2;\
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

#define ninline_stage4_Nodd(ct, offset)	t2 = (k-(offset))*bs+l1*dof+l2*2;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        mc1 = _mm256_loadu_pd(ct+t2+2);\
                        mc0 = _mm256_mul_pd(mx0,mc0);\
                        mc1 = _mm256_mul_pd(mx0,mc1);\
			msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc1));\

#define ninline_stage6_Nodd(ct, offset)   t2 = (k-(offset))*bs+l1*dof+l2*2;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        msum[l4] = _mm256_add_pd(msum[l4], _mm256_mul_pd(mx0,mc0))


#define ninline_stage7_Nodd(ct, offset) t2 = (k-(offset))*bs+l1*dof + l2;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        mc0 = _mm256_mul_pd(mx0,mc0);\
                    	mc1 = _mm256_permute2f128_pd(mc0, mc0, 0x01);\
			mc0 = _mm256_hadd_pd(mc0, mc1);\
                        msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc0))



#define ninline_stage8_Nodd(ct, offset) t2 = (k-(offset))*bs+l1*dof+l2;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
			mc0 = _mm256_mul_pd(mx0, mc0);\
                        msum[l4] = _mm256_add_pd(msum[l4], _mm256_hadd_pd(mc0,mc0))



#define ninline_stage9_Nodd(ct, offset) t2 = (k-(offset))*bs+l1*dof + l2;\
                        mc0 = _mm256_loadu_pd(ct+t2);\
                        msum[l4] = _mm256_add_pd(msum[l4], _mm256_mul_pd(mx0,mc0));\

#define ninline_Nodd(offset)	\
                    for(l3 = lbeg[k1]; l3 < lend[k1]; l3++){ \
			for(l2 = 0; (l2+4) <= dof ; l2 += 4){\
				nsetup127_Nodd(xt[l3]);\
				for(l1=0, l4 = 0; (l1+4)<=dof; l1+= 4, l4++){\
					ninline_stage1_Nodd(ct[l3], offset);\
				}\
				for(;(l1+2)<=dof;l1+=2, l4++){\
					ninline_stage2_Nodd(ct[l3], offset);\
				}\
				for(;l1 < dof;l1++, l4++){\
					ninline_stage7_Nodd(ct[l3], offset);\
				}\
			}\
			for(; (l2+2) <= dof; l2 += 2){\
				nsetup348_Nodd(xt[l3]);\
				for(l1=0, l4 = 0; (l1+4)<=dof; l1+= 4, l4++){\
					ninline_stage3_Nodd(ct[l3], offset);\
				}\
				for(;(l1+2)<=dof;l1+=2, l4++){\
					ninline_stage4_Nodd(ct[l3], offset);\
				}\
				for(;l1 < dof;l1++, l4++){\
					ninline_stage8_Nodd(ct[l3], offset);\
				}\
			}\
			for(; l2 < dof; l2 ++){\
				nsetup569_Nodd(xt[l3]);\
				for(l1=0, l4 = 0; (l1+4)<=dof; l1+= 4, l4++){\
					ninline_stage5_Nodd(ct[l3], offset);\
				}\
				for(;(l1+2)<=dof;l1+=2, l4++){\
					ninline_stage6_Nodd(ct[l3], offset);\
				}\
				for(;l1 < dof;l1++, l4++){\
					ninline_stage9_Nodd(ct[l3], offset);\
				}\
			}\
		}

PetscErrorCode BSG_MatMult_Nodd(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt *ioff, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset, PetscInt nregion, const PetscInt * lbeg, const PetscInt *lend , const PetscInt *rstart){
    
	const PetscInt l3threshold = WORKINGSETSIZE / bs;
   	 PetscInt count, endval;

	PetscInt i,j,k,l,k1, it,t1, t2, l1,l2, l3, l4;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby4 = (dof-2)/4 + 2;
	__m256d mx0, mx1, msum[dofby4], mc0, mc1, mc2, mc3;
    __m256i xtemp1 = _mm256_set_epi32(0,0,0,0,-1,-1,-1,-1);
    __m256i xtemp2 = _mm256_set_epi32(0,0,0,0,0,0,-1,-1);

    const PetscScalar *xt[nos], *ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

    for(k1 = 0 ; k1 <  nregion; k1++){
        setup_ct(k1);
//#pragma omp parallel for if(OPENMPB) shared(xt,ct,y,xtemp1, xtemp2) private(t1,t2, i,j,k, l1, l2,l3,mx0, mx1, msum, mc0, mc1, mc2, mc3)
        for(k = rstart[k1]; k < rstart[k1+1]; k++){
		nsetup_Nodd();
		ninline_Nodd(rstart[k1]);
		nsave_Nodd();
        }
    }

PetscFunctionReturn(0);
}
#endif

#define ninline_1_av2() for(l = lbeg[k1] ; l < lend[k1] ; l++){\
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mx1 = _mm256_loadu_pd(xt[l] + t1 + 4);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    mc1 = _mm256_loadu_pd(ct[l] + t2 + 4);\
                    msum00 = _mm256_add_pd(msum00 , _mm256_mul_pd(mx0, mc0));\
                    msum01 = _mm256_add_pd(msum01 , _mm256_mul_pd(mx1, mc1));\
			}

#define nsetup_1_av2(offset) t1= k*dof; t2 = (k-(offset))*bs;\
                         msum00 = _mm256_set_pd(0,0,0,0);\
                         msum01 = _mm256_set_pd(0,0,0,0);\

#define nsave_1_av2()\
                  _mm256_storeu_pd(y + t1 + 0,msum00);\
                  _mm256_storeu_pd(y + t1 + 4,msum01)

#define ninline_1_av() for(l = lbeg[k1] ; l < lend[k1] ; l++){\
                    mx0 = _mm256_loadu_pd(xt[l] + t1 + 0);\
                    mc0 = _mm256_loadu_pd(ct[l] + t2 + 0);\
                    msum00 = _mm256_add_pd(msum00 , _mm256_mul_pd(mx0, mc0));\
			}

#define nsetup_1_av(offset) t1= k*dof; t2 = (k-(offset))*bs;\
                         msum00 = _mm256_set_pd(0,0,0,0);\

#define nsave_1_av()\
                  _mm256_storeu_pd(y + t1 + 0,msum00)

#define ninline_1() for(l = lbeg[k1] ; l < lend[k1] ; l++){\
			msum0 += xt[l][t1] * ct[l][t2];\
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
	
	__m256d msum00, msum01, mx0, mc0, mx1, mc1;	

	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;

	
	const PetscScalar *xt[nos],*ct[nos];
    for(k1 = 0; k1< nos; k1++){
        xt[k1] = x + ioff[k1] * dof;
    }

	
	for(k1 = 0; k1 < nregion; k1++){
		setup_ct(k1);
//#pragma omp parallel for if(OPENMPB) shared(xt,ct,y) private(t1,t2,l, msum00, mx0, mc0)
		for(k = rstart[k1] ; k+7<rstart[k1+1]; k+=8){
			nsetup_1_av2(rstart[k1]);
			ninline_1_av2();
			nsave_1_av2();
			
		}
		for(; k+3<rstart[k1+1]; k+=4){
			nsetup_1_av(rstart[k1]);
			ninline_1_av();
			nsave_1_av();
			
		}
		for(; k<rstart[k1+1]; k++){
			nsetup_1(rstart[k1]);
			ninline_1();
			nsave_1();
			
		}
	}

	PetscFunctionReturn(0);
}


