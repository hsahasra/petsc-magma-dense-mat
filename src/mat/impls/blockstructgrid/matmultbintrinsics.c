#include <string.h>
#include <immintrin.h>

//#define SPREFETCH
#define min(a,b) (a)<(b) ? (a) : (b)


#include "../src/mat/impls/blockstructgrid/matblockstructgrid.h"
#include "../src/mat/impls/blockstructgrid/commonfunctions.h"

/*  -------------------------------------------------------------------- 
     This file implements matrix multiplication for the structgrid data type. 
     The routine employs SSE/AVX intrinsics if they are available on the machine.
     Otherwise, the computations default to normal PetscScalar operations. 
     The instruction for fused addmultiply has not been implemented of date.

     Author: Chekuri S. Choudary, RNET
*/

#define inline_2_unroll_load(l,iter)	mx##iter##0 = _mm_loadu_pd(xt##l+t1);\
			mc##iter##0 = _mm_loadu_pd(ct##l+t2);\
			mc##iter##1 = _mm_loadu_pd(ct##l+t2+2)

#define inline_2_prefetch(m)	PetscPrefetchBlock(xt##m+t1+2,2,0,PETSC_PREFETCH_HINT_T2);\
					PetscPrefetchBlock(ct##m+t2+4,4,0,PETSC_PREFETCH_HINT_T2)

#define inline_2_unroll_compute(iter)	msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx##iter##0,mc##iter##0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx##iter##0,mc##iter##1))

#define inline_2_2(l0,l1) inline_2_unroll_load(l0,0) ;\
				inline_2_unroll_load(l1,1) ;\
				inline_2_unroll_compute(0) ;\
				inline_2_unroll_compute(1)

#define inline_2_3(l0,l1,l2) inline_2_unroll_load(l0,0) ;\
				inline_2_unroll_load(l1,1) ;\
				inline_2_unroll_load(l2,2) ;\
				inline_2_unroll_compute(0) ;\
				inline_2_unroll_compute(1);\
				inline_2_unroll_compute(2)

#define inline_2(l)  mx0 = _mm_loadu_pd(xt##l+t1);\
			mc0 = _mm_loadu_pd(ct##l+t2);\
			mc1 = _mm_loadu_pd(ct##l+t2+2);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc1))

#define setup_2(offset) t1= k*dof; t2 = (k-(offset))*bs;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0)

#define save_2() msum1 = _mm_hadd_pd(msum0,msum1);\
		_mm_storeu_pd(y+t1, msum1)

#define setupct(pos, iter)	count = stpoffset[pos] + iter*nos;\
				ct0 = ctl[count + 0];\
				ct1 = ctl[count + 1];\
				ct2 = ctl[count + 2];\
				ct3 = ctl[count + 3];\
				ct4 = ctl[count + 4];\
				ct5 = ctl[count + 5];\
				ct6 = ctl[count + 6]


PetscInt BSG_MatMult_2(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt k,k1,it,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, msum0, msum1, mc0, mc1;
	__m128d mx00, mc00, mc01;
	__m128d mx10, mc10, mc11;
	__m128d mx20, mc20, mc21;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
		ct5 = ctl[5];
		ct6 = ctl[6];
			
	for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
			setup_2(0);
			inline_2_3(3,4,5);
			inline_2(6);
			save_2();
		}
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_2(xt4,ct4);
		setup_2(1);
		inline_2_3(2,3,4);
		inline_2_2(5,6);
		save_2();
		}
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_2(xt3,ct3);
		setup_2(lda3);
		inline_2_3(1,2,3);
		inline_2_3(4,5,6);
		save_2();
		}
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_2(xt2,ct2);
		setup_2(lda2);
		inline_2_3(0,1,2);
		inline_2_3(3,4,5);
		inline_2(6);
		save_2();
		}
		}

		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		{
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_2(xt2,ct2);
		setup_2(lda1-lda2);
		inline_2_3(0,1,2);
		inline_2_3(3,4,5);
		save_2();
		}
		}

		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		{
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_2(xt2,ct2);
		setup_2(lda1-lda3);
		inline_2_3(0,1,2);
		inline_2_2(3,4);
		save_2();
		}
		}

		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		{
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_2(xt2,ct2);
			setup_2(lda1-1);
			inline_2_2(0,1);
			inline_2_2(2,3);
			save_2();
		}
		}
	PetscFunctionReturn(0);
}

#define inline_4(l,m) mx0 = _mm_loadu_pd(xt##l+t1);\
			mx1 = _mm_loadu_pd(xt##l+t1+2);\
			mc0 = _mm_loadu_pd(ct##l+t2);\
			mc1 = _mm_loadu_pd(ct##l+t2+2);\
			mc2 = _mm_loadu_pd(ct##l+t2+4);\
			mc3 = _mm_loadu_pd(ct##l+t2+6);\
			mc4 = _mm_loadu_pd(ct##l+t2+8);\
			mc5 = _mm_loadu_pd(ct##l+t2+10);\
			mc6 = _mm_loadu_pd(ct##l+t2+12);\
			mc7 = _mm_loadu_pd(ct##l+t2+14);\
			PetscPrefetchBlock(xt##m+t1,4,0,PETSC_PREFETCH_HINT_T2);\
			PetscPrefetchBlock(ct##m+t2,16,0,PETSC_PREFETCH_HINT_T2);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx1,mc1));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc2));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx1,mc3));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx0,mc4));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx1,mc5));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx0,mc6));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx1,mc7))

#define setup_4(offset) t1= k*dof; t2 = (k-(offset))*bs;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0);\
		msum2 =_mm_set_pd(0,0);\
		msum3 =_mm_set_pd(0,0)

#define save_4() msum0 = _mm_hadd_pd(msum0,msum1);\
		msum2 = _mm_hadd_pd(msum2,msum3);\
		_mm_storeu_pd(y+t1, msum0);\
		_mm_storeu_pd(y+t1+2, msum2)

PetscInt BSG_MatMult_4(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt k,l,k1, it, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	__m128d mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
	for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
		ct1 = ctl[1];
		ct2 = ctl[2];
		ct3 = ctl[3];
		ct4 = ctl[4];
		ct5 = ctl[5];
		ct6 = ctl[6];

	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
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
		for(k = k1; k < endval; k++)
		inline_4(xt4,ct4);
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
		for(k = k1; k < endval; k++)
		inline_4(xt3,ct3);
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
		for(k = k1; k < endval; k++)
		inline_4(xt2,ct2);
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
		for(k = k1; k < endval; k++)
		inline_4(xt2,ct2);
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
		for(k = k1; k < endval; k++)
		inline_4(xt2,ct2);
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
		for(k = k1; k < endval; k++)
		inline_4(xt2,ct2);
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


#define inline_6(l,m) 	mx0 = _mm_loadu_pd(xt##l+t1);\
			mx1 = _mm_loadu_pd(xt##l+t1+2);\
			mx2 = _mm_loadu_pd(xt##l+t1+4);\
			mc0 = _mm_loadu_pd(ct##l+t2);\
			mc1 = _mm_loadu_pd(ct##l+t2+2);\
			mc2 = _mm_loadu_pd(ct##l+t2+4);\
			mc3 = _mm_loadu_pd(ct##l+t2+6);\
			mc4 = _mm_loadu_pd(ct##l+t2+8);\
			mc5 = _mm_loadu_pd(ct##l+t2+10);\
			mc6 = _mm_loadu_pd(ct##l+t2+12);\
			mc7 = _mm_loadu_pd(ct##l+t2+14);\
			mc8 = _mm_loadu_pd(ct##l+t2+16);\
			mc9 = _mm_loadu_pd(ct##l+t2+18);\
			mc10 = _mm_loadu_pd(ct##l+t2+20);\
			mc11 = _mm_loadu_pd(ct##l+t2+22);\
			mc12 = _mm_loadu_pd(ct##l+t2+24);\
			mc13 = _mm_loadu_pd(ct##l+t2+26);\
			mc14 = _mm_loadu_pd(ct##l+t2+28);\
			mc15 = _mm_loadu_pd(ct##l+t2+30);\
			mc16 = _mm_loadu_pd(ct##l+t2+32);\
			mc17 = _mm_loadu_pd(ct##l+t2+34);\
			PetscPrefetchBlock(xt##m+t1,6,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m+t2,36,0,PETSC_PREFETCH_HINT_NTA);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx1,mc1));\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx2,mc2));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc3));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx1,mc4));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx2,mc5));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx0,mc6));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx1,mc7));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx2,mc8));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx0,mc9));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx1,mc10));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx2,mc11));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx0,mc12));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx1,mc13));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx2,mc14));\
			msum5 = _mm_add_pd(msum5, _mm_mul_pd(mx0,mc15));\
			msum5 = _mm_add_pd(msum5, _mm_mul_pd(mx1,mc16));\
			msum5 = _mm_add_pd(msum5, _mm_mul_pd(mx2,mc17))

#define setup_6(offset)  t1= k*dof; t2 = (k-(offset))*bs;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0);\
		msum2 =_mm_set_pd(0,0);\
		msum3 =_mm_set_pd(0,0);\
		msum4 =_mm_set_pd(0,0);\
		msum5 =_mm_set_pd(0,0)

#define save_6() msum0 = _mm_hadd_pd(msum0,msum1);\
		msum2 = _mm_hadd_pd(msum2,msum3);\
		msum4 = _mm_hadd_pd(msum4,msum5);\
		_mm_storeu_pd(y+t1, msum0);\
		_mm_storeu_pd(y+t1+2, msum2);\
		_mm_storeu_pd(y+t1+4, msum4)

PetscInt BSG_MatMult_6(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs,const PetscInt * stpoffset)
{
	PetscInt k,l,k1,it, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	__m128d mx0, mx1, mx2, msum0, msum1, msum2, msum3, msum4,msum5, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9,  mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
	for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
		ct1 = ctl[1];
		ct2 = ctl[2];
		ct3 = ctl[3];
		ct4 = ctl[4];
		ct5 = ctl[5];
		ct6 = ctl[6];
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_6(0);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
		}
			inline_6(xt6,ct6);
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_6(xt3,ct3);
		setup_6(1);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
		}
			inline_6(xt4,ct4);
			inline_6(xt5,ct5);
			inline_6(xt6,ct6);
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_6(xt2,ct2);
		setup_6(lda3);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
		}
			inline_6(xt4,ct4);
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_6(xt2,ct2);
		setup_6(lda2);
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
		}
		}

		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		{
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
		setup_6(lda1-lda2);
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
		save_6();
		}
		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		{
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		setup_6(lda1-lda3);
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
		save_6();
		}
		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		{
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		setup_6(lda1-1);
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
		save_6();
		}


#define setup_Neven()	for(i=0;i<dofby2;i++){\
				 msum[i] = _mm_set_pd(0,0);\
			 }\

#define save_Neven() for(i=0;i<dofby2;i++){\
				_mm_storeu_pd(y+k*dof+2*i,msum[i]);\
			 }\

#define setup12_Neven(xt)	t1= k*dof+l2;\
			mx0 = _mm_loadu_pd(xt+t1);\
                        mx1 = _mm_loadu_pd(xt+t1+2)

#define setup34_Neven(xt)	t1= k*dof+l2;\
			mx0 = _mm_loadu_pd(xt+t1)

#define inline_stage1_Neven(ct,offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
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


#define inline_stage3_Neven(ct, offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
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

#define inline_stage2_Neven(ct, offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage4_Neven(ct, offset)	t2 = (k-(offset))*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_Neven(l,m, offset) PetscPrefetchBlock(xt##m+t1,dof,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m+t2,bs,0,PETSC_PREFETCH_HINT_NTA);\
			for(l2 = 0; l2 < dof-2; l2 += 4){\
				setup12_Neven(xt##l);\
				for(l1=0; l1<dofby2-1; l1+= 2){\
					inline_stage1_Neven(ct##l, offset);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage2_Neven(ct##l, offset);\
				}\
			}\
			for(; l2 < dof; l2 += 2){\
				setup34_Neven(xt##l);\
				for(l1=0; l1<dofby2-1; l1+= 2){\
					inline_stage3_Neven(ct##l, offset);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage4_Neven(ct##l, offset);\
				}\
			}\

PetscInt BSG_MatMult_Neven(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt i,k,l,k1,it, t1, t2, l1,l2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = dof/2;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	__m128d mx0, mx1, msum[dofby2], mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
		ct0 = ctl[0];
		ct1 = ctl[1];
		ct2 = ctl[2];
		ct3 = ctl[3];
		ct4 = ctl[4];
		ct5 = ctl[5];
		ct6 = ctl[6];
	for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
				{
					inline_stage1_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct3);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct3);
				}
			}
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
		for(k = k1; k < endval; k++)
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct5);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct5);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct6);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct6);
				}
		{
		setup_Neven();
			inline_Neven(3,4,0);
			inline_Neven(4,5,0);
			inline_Neven(5,6,0);
			inline_Neven(6,7,0);
				{
					inline_stage1_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct2);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt2);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct2);
				}
			}
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct4);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct4);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct5);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct5);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct6);
				}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct6);
		save_Neven();
		}
				{
					inline_stage1_Neven(ct1);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct1);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt1);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct1);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct1);
				}
			}
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct3);
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct3);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct3);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct4);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct4);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct5);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct5);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct6);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct6);
		{
		setup_Neven();
			inline_Neven(2,3,1);
			inline_Neven(3,4,1);
			inline_Neven(4,5,1);
			inline_Neven(5,6,1);
			inline_Neven(6,7,1);
					inline_stage1_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct0);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt0);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct0);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt2);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct2);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt2);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct2);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct3);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct3);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct4);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct4);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct5);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct5);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct6);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt6);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct6);
				}
		save_Neven();
		}
				{
					inline_stage1_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct0);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt0);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct0);
				}
			}
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_Neven();
			inline_Neven(1,2,lda3);
			inline_Neven(2,3,lda3);
			inline_Neven(3,4,lda3);
			inline_Neven(4,5,lda3);
			inline_Neven(5,6,lda3);
			inline_Neven(6,7,lda3);
				{
					inline_stage1_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct2);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt2);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct2);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct3);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct3);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct4);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct4);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct5);
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct5);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt5);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct5);
		save_Neven();
		}
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
				{
					inline_stage1_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct0);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt0);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct0);
				}
			}
		{
		setup_Neven();
			inline_Neven(0,1,lda2);
			inline_Neven(1,2,lda2);
			inline_Neven(2,3,lda2);
			inline_Neven(3,4,lda2);
			inline_Neven(4,5,lda2);
			inline_Neven(5,6,lda2);
			inline_Neven(6,7,lda2);
		save_Neven();
		}
		}

		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		{
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_Neven();
			inline_Neven(0,1,lda1-lda2);
			inline_Neven(1,2,lda1-lda2);
			inline_Neven(2,3,lda1-lda2);
			inline_Neven(3,4,lda1-lda2);
			inline_Neven(4,5,lda1-lda2);
			inline_Neven(5,6,lda1-lda2);
				{
					inline_stage1_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct2);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt2);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct2);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct3);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct3);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct4);
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct4);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt4);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct4);
		save_Neven();
		}
		}

		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		{
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
		for(k = k1; k < endval; k++)
				{
					inline_stage1_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct0);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt0);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct0);
				}
			}
		{
		setup_Neven();
			inline_Neven(0,1,lda1-lda3);
			inline_Neven(1,2,lda1-lda3);
			inline_Neven(2,3,lda1-lda3);
			inline_Neven(3,4,lda1-lda3);
			inline_Neven(4,5,lda1-lda3);
		save_Neven();
		}
		}

		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		{
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_Neven();
			inline_Neven(0,1,lda1-1);
			inline_Neven(1,2,lda1-1);
			inline_Neven(2,3,lda1-1);
			inline_Neven(3,4,lda1-1);
				setup12_Neven(xt2);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct2);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt2);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct2);
				}
			}
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven(ct3);
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven(xt3);
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven(ct3);
				}
		save_Neven();
		}
		}
	PetscFunctionReturn(0);
}

#define inline_1(l) msum0 += xt##l[t1] * ct##l[t2]
#define inline_1_2(l1,l2) msum0 += xt##l1[t1] * ct##l1[t2] + xt##l2[t1] * ct##l2[t2]
#define inline_1_3(l1,l2,l3) msum0 += xt##l1[t1] * ct##l1[t2] + xt##l2[t1] * ct##l2[t2] + xt##l3[t1] * ct##l3[t2] 

#define setup_1(offset) t1= k*dof; t2 = (k-(offset))*bs;\
		msum0 = 0.0

#define save_1() y[t1] = msum0

PetscInt BSG_MatMult_1(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt k,l,k1, it, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	PetscScalar msum0;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt count, endval;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
		ct2 = ctl[2];
		ct3 = ctl[3];
		ct4 = ctl[4];
		ct5 = ctl[5];
		ct6 = ctl[6];
	
	
	for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_1(0);
			inline_1_2(3,4);
			inline_1_2(5,6);
		save_1();
		}
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_1(xt4,ct4);
		setup_1(1);
			inline_1_2(2,3);
			inline_1_3(4,5,6);
		save_1();
		}
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_1(xt3,ct3);
		setup_1(lda3);
			inline_1_3(1,2,3);
			inline_1_3(4,5,6);
		save_1();
		}
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_1(xt2,ct2);
		setup_1(lda2);
			inline_1_2(0,1);
			inline_1_2(2,3);
			inline_1_3(4,5,6);
		save_1();
		}
		}

		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		{
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_1(xt2,ct2);
		setup_1(lda1-lda2);
			inline_1_3(0,1,2);
			inline_1_3(3,4,5);
		save_1();
		}
		}

		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		{
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_1(xt2,ct2);
		setup_1(lda1-lda3);
			inline_1_2(0,1);
			inline_1_3(2,3,4);
		save_1();
		}
		}

		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		{
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_1(xt2,ct2);
		setup_1(lda1-1);
			inline_1_2(0,1);
			inline_1_2(2,3);
		save_1();
		}
		}
	PetscFunctionReturn(0);
}

#define inline_3(l) 	mx0 = _mm_loadu_pd(xt##l+t1);\
			mx1 = _mm_load1_pd(xt##l+t1+2);\
			mc0 = _mm_loadu_pd(ct##l+t2);\
			mc1 = _mm_loadu_pd(ct##l+t2+2);\
			mc2 = _mm_loadu_pd(ct##l+t2+4);\
			mc3 = _mm_loadu_pd(ct##l+t2+6);\
			mc4 = _mm_loadu_pd(ct##l+t2+8);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc1));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx0,mc2));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx1,mc3));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx1,mc4))

#define setup_3(offset) t1= k*dof; t2 = (k-(offset))*bs;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0);\
		msum2 =_mm_set_pd(0,0);\
		msum3 =_mm_set_pd(0,0);\
		msum4 =_mm_set_pd(0,0)

#define save_3() msum0 = _mm_hadd_pd(msum0,msum1);      \
  msum0 = _mm_add_pd(msum0,msum3);                      \
  msum2 = _mm_hadd_pd(msum2,msum2);                     \
  msum2 = _mm_add_pd(msum2,msum4);                      \
  _mm_storeu_pd(y+t1, msum0);                           \
  _mm_maskstore_pd(y+t1+2,(__m128d)xtemp,msum2)

PetscInt BSG_MatMult_3(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt k,k1,l,it, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, mx1, msum0, msum1, msum2,msum3, msum4, mc0, mc1, mc2, mc3, mc4;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
	PetscInt count = 0;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt endval;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
		ct3 = ctl[3];
		ct4 = ctl[4];
		ct5 = ctl[5];
		ct6 = ctl[6];
	
		for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_3(0);
			inline_3(3);
			inline_3(4);
			inline_3(5);
			inline_3(6);
		save_3();
		}
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_3(xt4,ct4);
		setup_3(1);
			inline_3(2);
			inline_3(3);
			inline_3(4);
			inline_3(5);
			inline_3(6);
		save_3();
		}
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_3(xt3,ct3);
		setup_3(lda3);
			inline_3(1);
			inline_3(2);
			inline_3(3);
			inline_3(4);
			inline_3(5);
			inline_3(6);
		save_3();
		}
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_3(xt2,ct2);
		setup_3(lda2);
			inline_3(0);
			inline_3(1);
			inline_3(2);
			inline_3(3);
			inline_3(4);
			inline_3(5);
			inline_3(6);
		save_3();
		}
		}

		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		{
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_3(xt2,ct2);
		setup_3(lda1-lda2);
			inline_3(0);
			inline_3(1);
			inline_3(2);
			inline_3(3);
			inline_3(4);
			inline_3(5);
		save_3();
		}
		}

		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		{
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_3(xt2,ct2);
		setup_3(lda1-lda3);
			inline_3(0);
			inline_3(1);
			inline_3(2);
			inline_3(3);
			inline_3(4);
		save_3();
		}
		}

		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		{
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
		for(k = k1; k < endval; k++)
			inline_3(xt2,ct2);
		setup_3(lda1-1);
			inline_3(0);
			inline_3(1);
			inline_3(2);
			inline_3(3);
		save_3();
		}
		}
	PetscFunctionReturn(0);
}

#define inline_stage1_5(l,m)	mx0 = _mm_loadu_pd(xt##l+t1);\
				mx1 = _mm_loadu_pd(xt##l+t1+2);\
				mx2 = _mm_load1_pd(xt##l+t1+4);\
			PetscPrefetchBlock(xt##m+t1,5,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m+t2,25,0,PETSC_PREFETCH_HINT_NTA);\
			mc0 = _mm_loadu_pd(ct##l+t2);\
                        mc1 = _mm_loadu_pd(ct##l+t2+2);\
			mc2 = _mm_loadu_pd(ct##l+t2+4);\
			mc3 = _mm_loadu_pd(ct##l+t2+6);\
			mc4 = _mm_loadu_pd(ct##l+t2+20);\
                        msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
                        msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx1,mc1));\
                        msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc2));\
                        msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx1,mc3));\
                        msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx2,mc4))

#define setup_stage1_5(offset) t1= k*dof; t2 = (k-(offset))*bs;\
                msum0 =_mm_set_pd(0,0);\
                msum1 =_mm_set_pd(0,0);\
                msum2 =_mm_set_pd(0,0)


#define save_stage1_5() msum0 = _mm_hadd_pd(msum0,msum1);\
		msum0 = _mm_add_pd(msum0,msum2);\
                _mm_storeu_pd(y+t1, msum0)

#define inline_stage2_5(l) 	mc0 = _mm_loadu_pd(ct##l+t2+8);\
			mc1 = _mm_loadu_pd(ct##l+t2+10);\
			mc2 = _mm_loadu_pd(ct##l+t2+12);\
			mc3 = _mm_loadu_pd(ct##l+t2+14);\
			mc4 = _mm_loadu_pd(ct##l+t2+16);\
			mc5 = _mm_loadu_pd(ct##l+t2+18);\
			mc6 = _mm_loadu_pd(ct##l+t2+22);\
			mc7 = _mm_loadu_pd(ct##l+t2+24);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx1,mc1));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc2));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx1,mc3));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx0,mc4));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx1,mc5));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx2,mc6));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx2,mc7))

#define setup_stage2_5(offset) t1= k*dof; t2 = (k-(offset))*bs;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0);\
		msum2 =_mm_set_pd(0,0);\
		msum3 =_mm_set_pd(0,0);\
		msum4 =_mm_set_pd(0,0)

#define save_stage2_5() msum0 = _mm_hadd_pd(msum0,msum1);\
		msum0 = _mm_add_pd(msum0,msum3);\
		msum2 = _mm_hadd_pd(msum2,msum2);\
		msum2 = _mm_add_pd(msum2,msum4);\
		_mm_storeu_pd(y+t1+2, msum0);\
		_mm_maskstore_pd(y+t1+4,(__m128d)xtemp,msum2)

PetscInt BSG_MatMult_5(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt k,l,k1,it, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, mx1,mx2, msum0, msum1, msum2,msum3, msum4, mc0, mc1, mc2, mc3, mc4,mc5, mc6, mc7;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
	PetscInt count = 0;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt endval;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
		ct3 = ctl[3];
		ct4 = ctl[4];
		ct5 = ctl[5];
		ct6 = ctl[6];
	
		
		for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_stage1_5(0);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5(0);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
			inline_stage2_5(ct6);
		save_stage2_5();
		}
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
			inline_stage1_5(xt6,ct6);
			inline_stage2_5(ct5);
			inline_stage2_5(ct6);
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_stage1_5(1);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
			inline_stage2_5(ct5);
			inline_stage2_5(ct6);
		save_stage2_5();
		}
			inline_stage1_5(xt4,ct4);
			inline_stage1_5(xt5,ct5);
			inline_stage2_5(ct4);
			inline_stage2_5(ct5);
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_stage1_5(lda3);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5(lda3);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
		save_stage2_5();
		}
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_stage1_5(lda2);
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5(lda2);
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
		save_stage2_5();
		}
		}
		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_stage1_5(lda1-lda2);
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
		save_stage1_5();
		setup_stage2_5(lda1-lda2);
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
		save_stage2_5();
		}
		}
		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_stage1_5(lda1-lda3);
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
		save_stage1_5();
		setup_stage2_5(lda1-lda3);
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
		save_stage2_5();
		}
		}
		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_stage1_5(lda1-1);
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
		save_stage1_5();
		setup_stage2_5(lda1-1);
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
		save_stage2_5();
		}
		}

#define setup_Nodd()   for(i=0;i<dofby2;i++){\
                                 msum[i] = _mm_set_pd(0,0);\
                         }\

#define save_Nodd() for(i=0;i<dofby2-1;i++){\
                                _mm_storeu_pd(y+k*dof+2*i,msum[i]);\
                         }\
  _mm_maskstore_pd(y+(k+1)*dof-1,(__m128d)xtemp,msum[dofby2-1])

#define setup123_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt+t1);\
                        mx1 = _mm_loadu_pd(xt+t1+2)

#define setup456_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt+t1)

#define setup789_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_load1_pd(xt+t1)

#define inline_stage1_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
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


#define inline_stage4_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
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

#define inline_stage7_Nodd(ct, offset) t2 = (k-(offset))*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0);\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], mc1)

#define inline_stage2_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+1);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage5_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage8_Nodd(ct, offset) t2 = (k-(offset))*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0)

#define inline_stage3_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define inline_stage6_Nodd(ct, offset)   t2 = (k-(offset))*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define inline_stage9_Nodd(ct, offset) t2 = (k-(offset))*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0);\

#define inline_Nodd(l,m, offset)	PetscPrefetchBlock(xt##m+t1,dof,0,PETSC_PREFETCH_HINT_NTA);\
				PetscPrefetchBlock(ct##m+t2,bs,0,PETSC_PREFETCH_HINT_NTA);\
				for(l2 = 0; l2 < dof-3; l2 += 4){\
				setup123_Nodd(xt##l);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					inline_stage1_Nodd(ct##l, offset);\
				}\
				for(;l1<dofby2-1;l1++){\
					inline_stage2_Nodd(ct##l, offset);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage3_Nodd(ct##l, offset);\
				}\
			}\
			for(; l2 < dof-1; l2 += 2){\
				setup456_Nodd(xt##l);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					inline_stage4_Nodd(ct##l, offset);\
				}\
				for(;l1<dofby2-1;l1++){\
					inline_stage5_Nodd(ct##l, offset);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage6_Nodd(ct##l, offset);\
				}\
			}\
			for(; l2 < dof; l2 ++){\
				setup789_Nodd(xt##l);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					inline_stage7_Nodd(ct##l, offset);\
				}\
				for(;l1<dofby2-1;l1++){\
					inline_stage8_Nodd(ct##l, offset);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage9_Nodd(ct##l, offset);\
				}\
			}\

PetscInt BSG_MatMult_Nodd(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs, const PetscInt * stpoffset)
{
	PetscInt i,k,l,k1, it,t1, t2, l1,l2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = (dof+1)/2;
	__m128d mx0, mx1, msum[dofby2], mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
	PetscInt count = 0;
	PetscInt l3threshold = WORKINGSETSIZE / bs;
	PetscInt endval;
	const PetscScalar *xt0,*xt1,*xt2,*xt3,*xt4,*xt5,*xt6,*ct0,*ct1,*ct2,*ct3,*ct4,*ct5,*ct6,*xt7,*ct7;
		xt0 = x + (idx[0] + idy[0]*lda3 + idz[0]*lda2)*dof;
		xt1 = x + (idx[1] + idy[1]*lda3 + idz[1]*lda2)*dof;
		xt2 = x + (idx[2] + idy[2]*lda3 + idz[2]*lda2)*dof;
		xt3 = x + (idx[3] + idy[3]*lda3 + idz[3]*lda2)*dof;
		xt4 = x + (idx[4] + idy[4]*lda3 + idz[4]*lda2)*dof;
		xt5 = x + (idx[5] + idy[5]*lda3 + idz[5]*lda2)*dof;
		xt6 = x + (idx[6] + idy[6]*lda3 + idz[6]*lda2)*dof;
		ct3 = ctl[3];
		ct4 = ctl[4];
		ct5 = ctl[5];
		ct6 = ctl[6];
		
		for(k1 = 0, it =0; k1 < 1; k1+=l3threshold, it++)
	{
		setupct(0,it);
		endval = min(1,k1+l3threshold);
		for(k = k1; k < endval; k++)
				{
					inline_stage1_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct3);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct3);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct3);
				}
			}
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct5);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct5);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct5);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct6);
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct6);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct6);
				}
		{
		setup_Nodd();
			inline_Nodd(3,4,0);
			inline_Nodd(4,5,0);
			inline_Nodd(5,6,0);
			inline_Nodd(6,7,0);
				{
					inline_stage1_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct2);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct2);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct2);
				}
			}
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct4);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct4);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct4);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct5);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct5);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct5);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct6);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct6);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct6);
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct6);
		save_Nodd();
		}
				{
					inline_stage1_Nodd(ct1);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct1);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct1);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt1);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct1);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct1);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct1);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt1);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct1);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct1);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct1);
				}
			}
	}
		for(k1 = 1, it =0; k1 < lda3; k1+=l3threshold, it++)
		{
		setupct(1,it);
		endval = min(lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct3);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct3);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct3);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct4);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct4);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct4);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct5);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct5);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct5);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct6);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct6);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct6);
				}
		setup_Nodd();
			inline_Nodd(2,3,1);
			inline_Nodd(3,4,1);
			inline_Nodd(4,5,1);
			inline_Nodd(5,6,1);
			inline_Nodd(6,7,1);
					inline_stage1_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct0);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct0);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct0);
				}
			}
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct2);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct2);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct2);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct3);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct3);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct3);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct4);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct4);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct4);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct5);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct5);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct5);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct6);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct6);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct6);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt6);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct6);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct6);
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct6);
		save_Nodd();
		}
				{
					inline_stage1_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct0);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct0);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct0);
				}
			}
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct2);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct2);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct2);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct3);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct3);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct3);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct4);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct4);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct4);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct5);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct5);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct5);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt5);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct5);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct5);
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct5);
		}

		for(k1 = lda3, it =0; k1 < lda2; k1+=l3threshold, it++)
		{
		setupct(2,it);
		endval = min(lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
				{
					inline_stage1_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct0);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct0);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct0);
				}
			}
		{
		setup_Nodd();
			inline_Nodd(1,2,lda3);
			inline_Nodd(2,3,lda3);
			inline_Nodd(3,4,lda3);
			inline_Nodd(4,5,lda3);
			inline_Nodd(5,6,lda3);
			inline_Nodd(6,7,lda3);
		save_Nodd();
		}
		}

		for(k1 = lda2, it =0; k1 < lda1-lda2; k1+=l3threshold, it++)
		{
		setupct(3,it);
		endval = min(lda1-lda2,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_Nodd();
			inline_Nodd(0,1,lda2);
			inline_Nodd(1,2,lda2);
			inline_Nodd(2,3,lda2);
			inline_Nodd(3,4,lda2);
			inline_Nodd(4,5,lda2);
			inline_Nodd(5,6,lda2);
			inline_Nodd(6,7,lda2);
					inline_stage1_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct2);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct2);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct2);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct3);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct3);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct3);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct4);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct4);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct4);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt4);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct4);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct4);
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct4);
		save_Nodd();
		}
				{
					inline_stage1_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct0);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct0);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt0);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct0);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct0);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct0);
				}
			}
		}

		for(k1 = lda1-lda2, it =0; k1 < lda1-lda3; k1+=l3threshold, it++)
		{
		setupct(4,it);
		endval = min(lda1-lda3,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_Nodd();
			inline_Nodd(0,1,lda1-lda2);
			inline_Nodd(1,2,lda1-lda2);
			inline_Nodd(2,3,lda1-lda2);
			inline_Nodd(3,4,lda1-lda2);
			inline_Nodd(4,5,lda1-lda2);
			inline_Nodd(5,6,lda1-lda2);
		save_Nodd();
		}
		}

		for(k1 = lda1-lda3, it =0; k1 < lda1-1; k1+=l3threshold, it++)
		{
		setupct(5,it);
		endval = min(lda1-1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_Nodd();
			inline_Nodd(0,1,lda1-lda3);
			inline_Nodd(1,2,lda1-lda3);
			inline_Nodd(2,3,lda1-lda3);
			inline_Nodd(3,4,lda1-lda3);
			inline_Nodd(4,5,lda1-lda3);
		save_Nodd();
		}
		}

		for(k1 = lda1-1, it =0; k1 < lda1; k1+=l3threshold, it++)
		{
		setupct(6,it);
		endval = min(lda1,k1+l3threshold);
		for(k = k1; k < endval; k++)
		{
		setup_Nodd();
			inline_Nodd(0,1, lda1-1);
			inline_Nodd(1,2, lda1-1);
			inline_Nodd(2,3, lda1-1);
			inline_Nodd(3,4, lda1-1);
				setup123_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct2);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct2);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct2);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt2);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct2);
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct2);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct2);
				}
			}
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd(ct3);
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd(ct3);
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd(xt3);
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd(ct3);
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd(ct3);
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd(ct3);
		save_Nodd();
		}
		}
	PetscFunctionReturn(0);
}
