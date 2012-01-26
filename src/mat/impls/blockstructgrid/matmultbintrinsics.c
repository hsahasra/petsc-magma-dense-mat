#include <string.h>
#include <immintrin.h>

//#define SPREFETCH



#include "../src/mat/impls/blockstructgrid/matblockstructgrid.h"
#include "../src/mat/impls/blockstructgrid/commonfunctions.h"

/*  -------------------------------------------------------------------- 
     This file implements matrix multiplication for the structgrid data type. 
     The routine employs SSE/AVX intrinsics if they are available on the machine.
     Otherwise, the computations default to normal PetscScalar operations. 
     The instruction for fused addmultiply has not been implemented of date.

     Author: Chekuri S. Choudary, RNET
*/



#define inline_2(l,m)  mx0 = _mm_loadu_pd(xt##l+t1);\
			mc0 = _mm_loadu_pd(ct##l+t2);\
			mc1 = _mm_loadu_pd(ct##l+t2+2);\
			PetscPrefetchBlock(xt##m,2,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m,4,0,PETSC_PREFETCH_HINT_NTA);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc1))

#define setup_2() t1= k*dof; t2 = t1*dof;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0)

#define save_2() msum1 = _mm_hadd_pd(msum0,msum1);\
		_mm_storeu_pd(y+t1, msum1)

PetscInt BSG_MatMult_2(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, msum0, msum1, mc0, mc1;
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
			
	for(k = 0; k < 1; k++)
	{
		setup_2();
			inline_2(3,4);
			inline_2(4,5);
			inline_2(5,6);
			inline_2(6,7);
		save_2();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_2();
			inline_2(2,3);
			inline_2(3,4);
			inline_2(4,5);
			inline_2(5,6);
			inline_2(6,7);
		save_2();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_2();
			inline_2(1,2);
			inline_2(2,3);
			inline_2(3,4);
			inline_2(4,5);
			inline_2(5,6);
			inline_2(6,7);
		save_2();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_2();
			inline_2(0,1);
			inline_2(1,2);
			inline_2(2,3);
			inline_2(3,4);
			inline_2(4,5);
			inline_2(5,6);
			inline_2(6,7);
		save_2();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_2();
			inline_2(0,1);
			inline_2(1,2);
			inline_2(2,3);
			inline_2(3,4);
			inline_2(4,5);
			inline_2(5,6);
		save_2();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_2();
			inline_2(0,1);
			inline_2(1,2);
			inline_2(2,3);
			inline_2(3,4);
			inline_2(4,5);
		save_2();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_2();
			inline_2(0,1);
			inline_2(1,2);
			inline_2(2,3);
			inline_2(3,4);
		save_2();
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
			PetscPrefetchBlock(xt##m,4,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m,16,0,PETSC_PREFETCH_HINT_NTA);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx1,mc1));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc2));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx1,mc3));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx0,mc4));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx1,mc5));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx0,mc6));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx1,mc7))

#define setup_4() t1= k*dof; t2 = t1*dof;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0);\
		msum2 =_mm_set_pd(0,0);\
		msum3 =_mm_set_pd(0,0)

#define save_4() msum0 = _mm_hadd_pd(msum0,msum1);\
		msum2 = _mm_hadd_pd(msum2,msum3);\
		_mm_storeu_pd(y+t1, msum0);\
		_mm_storeu_pd(y+t1+2, msum2)

PetscInt BSG_MatMult_4(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, mx1, msum0, msum1, msum2, msum3, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
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

	
	for(k = 0; k < 1; k++)
	{
		setup_4();
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_4();
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_4();
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_4();
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		inline_4(6,7);
		save_4();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_4();
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		inline_4(5,6);
		save_4();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_4();
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		inline_4(4,5);
		save_4();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_4();
		inline_4(0,1);
		inline_4(1,2);
		inline_4(2,3);
		inline_4(3,4);
		save_4();
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
			PetscPrefetchBlock(xt##m,6,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m,36,0,PETSC_PREFETCH_HINT_NTA);\
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

#define setup_6()  t1= k*dof; t2 = t1*dof;\
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

PetscInt BSG_MatMult_6(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, mx1, mx2, msum0, msum1, msum2, msum3, msum4,msum5, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9,  mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17;
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
	
	for(k = 0; k < 1; k++)
	{
		setup_6();
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_6();
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_6();
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_6();
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
			inline_6(6,7);
		save_6();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_6();
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
			inline_6(5,6);
		save_6();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_6();
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
			inline_6(4,5);
		save_6();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_6();
			inline_6(0,1);
			inline_6(1,2);
			inline_6(2,3);
			inline_6(3,4);
		save_6();
	}
	PetscFunctionReturn(0);
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

#define inline_stage1_Neven(ct)	t2 = k*bs+2*l1*dof+l2;\
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


#define inline_stage3_Neven(ct)	t2 = k*bs+2*l1*dof+l2;\
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

#define inline_stage2_Neven(ct)	t2 = k*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage4_Neven(ct)	t2 = k*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_Neven(l,m) PetscPrefetchBlock(xt##m,dof,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m,bs,0,PETSC_PREFETCH_HINT_NTA);\
			for(l2 = 0; l2 < dof-2; l2 += 4){\
				setup12_Neven(xt##l);\
				for(l1=0; l1<dofby2-1; l1+= 2){\
					inline_stage1_Neven(ct##l);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage2_Neven(ct##l);\
				}\
			}\
			for(; l2 < dof; l2 += 2){\
				setup34_Neven(xt##l);\
				for(l1=0; l1<dofby2-1; l1+= 2){\
					inline_stage3_Neven(ct##l);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage4_Neven(ct##l);\
				}\
			}\

PetscInt BSG_MatMult_Neven(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt i,k,l, t1, t2, l1,l2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = dof/2;
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
	
	for(k = 0; k < 1; k++)
	{
		setup_Neven();
			inline_Neven(3,4);
			inline_Neven(4,5);
			inline_Neven(5,6);
			inline_Neven(6,7);
		save_Neven();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_Neven();
			inline_Neven(2,3);
			inline_Neven(3,4);
			inline_Neven(4,5);
			inline_Neven(5,6);
			inline_Neven(6,7);
		save_Neven();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_Neven();
			inline_Neven(1,2);
			inline_Neven(2,3);
			inline_Neven(3,4);
			inline_Neven(4,5);
			inline_Neven(5,6);
			inline_Neven(6,7);
		save_Neven();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_Neven();
			inline_Neven(0,1);
			inline_Neven(1,2);
			inline_Neven(2,3);
			inline_Neven(3,4);
			inline_Neven(4,5);
			inline_Neven(5,6);
			inline_Neven(6,7);
		save_Neven();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_Neven();
			inline_Neven(0,1);
			inline_Neven(1,2);
			inline_Neven(2,3);
			inline_Neven(3,4);
			inline_Neven(4,5);
			inline_Neven(5,6);
		save_Neven();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_Neven();
			inline_Neven(0,1);
			inline_Neven(1,2);
			inline_Neven(2,3);
			inline_Neven(3,4);
			inline_Neven(4,5);
		save_Neven();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_Neven();
			inline_Neven(0,1);
			inline_Neven(1,2);
			inline_Neven(2,3);
			inline_Neven(3,4);
		save_Neven();
	}
	PetscFunctionReturn(0);
}

#define inline_1(l) msum0 += xt##l[t1] * ct##l[t2]

#define setup_1() t1= k*dof; t2 = t1*dof;\
		msum0 = 0.0

#define save_1() y[t1] = msum0

PetscInt BSG_MatMult_1(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	PetscScalar msum0;
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
	
	
	for(k = 0; k < 1; k++)
	{
		setup_1();
			inline_1(3);
			inline_1(4);
			inline_1(5);
			inline_1(6);
		save_1();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_1();
			inline_1(2);
			inline_1(3);
			inline_1(4);
			inline_1(5);
			inline_1(6);
		save_1();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_1();
			inline_1(1);
			inline_1(2);
			inline_1(3);
			inline_1(4);
			inline_1(5);
			inline_1(6);
		save_1();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_1();
			inline_1(0);
			inline_1(1);
			inline_1(2);
			inline_1(3);
			inline_1(4);
			inline_1(5);
			inline_1(6);
		save_1();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_1();
			inline_1(0);
			inline_1(1);
			inline_1(2);
			inline_1(3);
			inline_1(4);
			inline_1(5);
		save_1();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_1();
			inline_1(0);
			inline_1(1);
			inline_1(2);
			inline_1(3);
			inline_1(4);
		save_1();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_1();
			inline_1(0);
			inline_1(1);
			inline_1(2);
			inline_1(3);
		save_1();
	}
	PetscFunctionReturn(0);
}

#define inline_3(l,m) 	mx0 = _mm_loadu_pd(xt##l+t1);\
			mx1 = _mm_load1_pd(xt##l+t1+2);\
			mc0 = _mm_loadu_pd(ct##l+t2);\
			mc1 = _mm_loadu_pd(ct##l+t2+2);\
			mc2 = _mm_loadu_pd(ct##l+t2+4);\
			mc3 = _mm_loadu_pd(ct##l+t2+6);\
			mc4 = _mm_loadu_pd(ct##l+t2+8);\
			PetscPrefetchBlock(xt##m,3,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m,9,0,PETSC_PREFETCH_HINT_NTA);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc1));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx0,mc2));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx1,mc3));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx1,mc4))

#define setup_3() t1= k*dof; t2 = t1*dof;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0);\
		msum2 =_mm_set_pd(0,0);\
		msum3 =_mm_set_pd(0,0);\
		msum4 =_mm_set_pd(0,0)

#define save_3() msum0 = _mm_hadd_pd(msum0,msum1);\
		msum0 = _mm_add_pd(msum0,msum3);\
		msum2 = _mm_hadd_pd(msum2,msum2);\
		msum2 = _mm_add_pd(msum2,msum4);\
		_mm_storeu_pd(y+t1, msum0);\
		_mm_maskstore_pd(y+t1+2,xtemp,msum2)

PetscInt BSG_MatMult_3(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, mx1, msum0, msum1, msum2,msum3, msum4, mc0, mc1, mc2, mc3, mc4;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
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
	
	for(k = 0; k < 1; k++)
	{
		setup_3();
			inline_3(3,4);
			inline_3(4,5);
			inline_3(5,6);
			inline_3(6,7);
		save_3();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_3();
			inline_3(2,3);
			inline_3(3,4);
			inline_3(4,5);
			inline_3(5,6);
			inline_3(6,7);
		save_3();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_3();
			inline_3(1,2);
			inline_3(2,3);
			inline_3(3,4);
			inline_3(4,5);
			inline_3(5,6);
			inline_3(6,7);
		save_3();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_3();
			inline_3(0,1);
			inline_3(1,2);
			inline_3(2,3);
			inline_3(3,4);
			inline_3(4,5);
			inline_3(5,6);
			inline_3(6,7);
		save_3();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_3();
			inline_3(0,1);
			inline_3(1,2);
			inline_3(2,3);
			inline_3(3,4);
			inline_3(4,5);
			inline_3(5,6);
		save_3();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_3();
			inline_3(0,1);
			inline_3(1,2);
			inline_3(2,3);
			inline_3(3,4);
			inline_3(4,5);
		save_3();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_3();
			inline_3(0,1);
			inline_3(1,2);
			inline_3(2,3);
			inline_3(3,4);
		save_3();
	}
	PetscFunctionReturn(0);
}

#define inline_stage1_5(l,m)	mx0 = _mm_loadu_pd(xt##l+t1);\
				mx1 = _mm_loadu_pd(xt##l+t1+2);\
				mx2 = _mm_load1_pd(xt##l+t1+4);\
			PetscPrefetchBlock(xt##m,5,0,PETSC_PREFETCH_HINT_NTA);\
			PetscPrefetchBlock(ct##m,25,0,PETSC_PREFETCH_HINT_NTA);\
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

#define setup_stage1_5() t1= k*dof; t2 = t1*dof;\
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

#define setup_stage2_5() t1= k*dof; t2 = t1*dof;\
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
		_mm_maskstore_pd(y+t1+4,xtemp,msum2)

PetscInt BSG_MatMult_5(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	__m128d mx0, mx1,mx2, msum0, msum1, msum2,msum3, msum4, mc0, mc1, mc2, mc3, mc4,mc5, mc6, mc7;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
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
	
	for(k = 0; k < 1; k++)
	{
		setup_stage1_5();
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5();
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
		save_stage2_5();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_stage1_5();
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5();
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
		save_stage2_5();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_stage1_5();
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5();
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
		save_stage2_5();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_stage1_5();
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
			inline_stage1_5(6,7);
		save_stage1_5();
		setup_stage2_5();
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
			inline_stage2_5(6);
		save_stage2_5();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_stage1_5();
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
			inline_stage1_5(5,6);
		save_stage1_5();
		setup_stage2_5();
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
			inline_stage2_5(5);
		save_stage2_5();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_stage1_5();
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
			inline_stage1_5(4,5);
		save_stage1_5();
		setup_stage2_5();
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
			inline_stage2_5(4);
		save_stage2_5();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_stage1_5();
			inline_stage1_5(0,1);
			inline_stage1_5(1,2);
			inline_stage1_5(2,3);
			inline_stage1_5(3,4);
		save_stage1_5();
		setup_stage2_5();
			inline_stage2_5(0);
			inline_stage2_5(1);
			inline_stage2_5(2);
			inline_stage2_5(3);
		save_stage2_5();
	}
	
	PetscFunctionReturn(0);
}


#define setup_Nodd()   for(i=0;i<dofby2;i++){\
                                 msum[i] = _mm_set_pd(0,0);\
                         }\

#define save_Nodd() for(i=0;i<dofby2-1;i++){\
                                _mm_storeu_pd(y+k*dof+2*i,msum[i]);\
                         }\
			_mm_maskstore_pd(y+(k+1)*dof-1,xtemp,msum[dofby2-1])

#define setup123_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt+t1);\
                        mx1 = _mm_loadu_pd(xt+t1+2)

#define setup456_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt+t1)

#define setup789_Nodd(xt) t1= k*dof+l2;\
                        mx0 = _mm_load1_pd(xt+t1)

#define inline_stage1_Nodd(ct)   t2 = k*bs+2*l1*(dof-1)+l2;\
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


#define inline_stage4_Nodd(ct)   t2 = k*bs+2*l1*(dof-1)+l2;\
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

#define inline_stage7_Nodd(ct) t2 = k*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0);\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], mc1)

#define inline_stage2_Nodd(ct)   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc3 = _mm_loadu_pd(ct+t2+dof+1);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage5_Nodd(ct)   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc2 = _mm_loadu_pd(ct+t2+dof-1);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage8_Nodd(ct) t2 = k*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0)

#define inline_stage3_Nodd(ct)   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc1 = _mm_loadu_pd(ct+t2+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define inline_stage6_Nodd(ct)   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define inline_stage9_Nodd(ct) t2 = k*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0);\

#define inline_Nodd(l,m)	PetscPrefetchBlock(xt##m,dof,0,PETSC_PREFETCH_HINT_NTA);\
				PetscPrefetchBlock(ct##m,bs,0,PETSC_PREFETCH_HINT_NTA);\
				for(l2 = 0; l2 < dof-3; l2 += 4){\
				setup123_Nodd(xt##l);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					inline_stage1_Nodd(ct##l);\
				}\
				for(;l1<dofby2-1;l1++){\
					inline_stage2_Nodd(ct##l);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage3_Nodd(ct##l);\
				}\
			}\
			for(; l2 < dof-1; l2 += 2){\
				setup456_Nodd(xt##l);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					inline_stage4_Nodd(ct##l);\
				}\
				for(;l1<dofby2-1;l1++){\
					inline_stage5_Nodd(ct##l);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage6_Nodd(ct##l);\
				}\
			}\
			for(; l2 < dof; l2 ++){\
				setup789_Nodd(xt##l);\
				for(l1=0; l1<dofby2-2; l1+= 2){\
					inline_stage7_Nodd(ct##l);\
				}\
				for(;l1<dofby2-1;l1++){\
					inline_stage8_Nodd(ct##l);\
				}\
				for(;l1<dofby2;l1++){\
					inline_stage9_Nodd(ct##l);\
				}\
			}\

PetscInt BSG_MatMult_Nodd(PetscScalar ** ctl,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt i,k,l, t1, t2, l1,l2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = (dof+1)/2;
	__m128d mx0, mx1, msum[dofby2], mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
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
	
	for(k = 0; k < 1; k++)
	{
		setup_Nodd();
			inline_Nodd(3,4);
			inline_Nodd(4,5);
			inline_Nodd(5,6);
			inline_Nodd(6,7);
		save_Nodd();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_Nodd();
			inline_Nodd(2,3);
			inline_Nodd(3,4);
			inline_Nodd(4,5);
			inline_Nodd(5,6);
			inline_Nodd(6,7);
		save_Nodd();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_Nodd();
			inline_Nodd(1,2);
			inline_Nodd(2,3);
			inline_Nodd(3,4);
			inline_Nodd(4,5);
			inline_Nodd(5,6);
			inline_Nodd(6,7);
		save_Nodd();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_Nodd();
			inline_Nodd(0,1);
			inline_Nodd(1,2);
			inline_Nodd(2,3);
			inline_Nodd(3,4);
			inline_Nodd(4,5);
			inline_Nodd(5,6);
			inline_Nodd(6,7);
		save_Nodd();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_Nodd();
			inline_Nodd(0,1);
			inline_Nodd(1,2);
			inline_Nodd(2,3);
			inline_Nodd(3,4);
			inline_Nodd(4,5);
			inline_Nodd(5,6);
		save_Nodd();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_Nodd();
			inline_Nodd(0,1);
			inline_Nodd(1,2);
			inline_Nodd(2,3);
			inline_Nodd(3,4);
			inline_Nodd(4,5);
		save_Nodd();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_Nodd();
			inline_Nodd(0,1);
			inline_Nodd(1,2);
			inline_Nodd(2,3);
			inline_Nodd(3,4);
		save_Nodd();
	}
	PetscFunctionReturn(0);
}
