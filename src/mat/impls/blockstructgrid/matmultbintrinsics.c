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


#define inline_2()  mx0 = _mm_loadu_pd(xt[l]+t1);\
			mc0 = _mm_loadu_pd(ct[l]+t2);\
			mc1 = _mm_loadu_pd(ct[l]+t2+2);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx0,mc0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx0,mc1))

#define setup_2() t1= k*dof; t2 = t1*dof;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0)

#define save_2() msum1 = _mm_hadd_pd(msum0,msum1);\
		_mm_storeu_pd(y+t1, msum1)

PetscInt BSG_MatMult_2(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscScalar *xt[7];
	__m128d mx0, msum0, msum1, mc0, mc1;
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_2();
                for(l=mnos;l<nos;l++)
                {
			inline_2();
                }
		save_2();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_2();
		for(l=mnos-1;l<nos;l++)
		{
			inline_2();
		}
		save_2();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_2();
		for(l=mnos-2;l<nos;l++)
		{
			inline_2();
		}
		save_2();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_2();
		for(l=0;l<nos;l++)
		{
			inline_2();
		}
		save_2();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_2();
		for(l=0;l<nos-1;l++)
		{
			inline_2();
		}
		save_2();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_2();
		for(l=0;l<nos-2;l++)
		{
			inline_2();
		}
		save_2();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_2();
		for(l=0;l<=mnos;l++)
		{
			inline_2();
		}
		save_2();
	}
	PetscFunctionReturn(0);
}

#define inline_4() mx0 = _mm_loadu_pd(xt[l]+t1);\
			mx1 = _mm_loadu_pd(xt[l]+t1+2);\
			mc0 = _mm_loadu_pd(ct[l]+t2);\
			mc1 = _mm_loadu_pd(ct[l]+t2+2);\
			mc2 = _mm_loadu_pd(ct[l]+t2+4);\
			mc3 = _mm_loadu_pd(ct[l]+t2+6);\
			mc4 = _mm_loadu_pd(ct[l]+t2+8);\
			mc5 = _mm_loadu_pd(ct[l]+t2+10);\
			mc6 = _mm_loadu_pd(ct[l]+t2+12);\
			mc7 = _mm_loadu_pd(ct[l]+t2+14);\
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

PetscInt BSG_MatMult_4(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscScalar *xt[7];
	__m128d mx0, mx1, msum0, msum1, msum2, msum3, msumf1,msumf2, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_4();
                for(l=mnos;l<nos;l++)
                {
			inline_4();
                }
		save_4();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_4();
		for(l=mnos-1;l<nos;l++)
		{
			inline_4();
		}
		save_4();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_4();
		for(l=mnos-2;l<nos;l++)
		{
			inline_4();
		}
		save_4();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_4();
		for(l=0;l<nos;l++)
		{
			inline_4();
		}
		save_4();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_4();
		for(l=0;l<nos-1;l++)
		{
			inline_4();
		}
		save_4();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_4();
		for(l=0;l<nos-2;l++)
		{
			inline_4();
		}
		save_4();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_4();
		for(l=0;l<=mnos;l++)
		{
			inline_4();
		}
		save_4();
	}
	PetscFunctionReturn(0);
}


#define inline_6() 	mx0 = _mm_loadu_pd(xt[l]+t1);\
			mx1 = _mm_loadu_pd(xt[l]+t1+2);\
			mx2 = _mm_loadu_pd(xt[l]+t1+4);\
			mc0 = _mm_loadu_pd(ct[l]+t2);\
			mc1 = _mm_loadu_pd(ct[l]+t2+2);\
			mc2 = _mm_loadu_pd(ct[l]+t2+4);\
			mc3 = _mm_loadu_pd(ct[l]+t2+6);\
			mc4 = _mm_loadu_pd(ct[l]+t2+8);\
			mc5 = _mm_loadu_pd(ct[l]+t2+10);\
			mc6 = _mm_loadu_pd(ct[l]+t2+12);\
			mc7 = _mm_loadu_pd(ct[l]+t2+14);\
			mc8 = _mm_loadu_pd(ct[l]+t2+16);\
			mc9 = _mm_loadu_pd(ct[l]+t2+18);\
			mc10 = _mm_loadu_pd(ct[l]+t2+20);\
			mc11 = _mm_loadu_pd(ct[l]+t2+22);\
			mc12 = _mm_loadu_pd(ct[l]+t2+24);\
			mc13 = _mm_loadu_pd(ct[l]+t2+26);\
			mc14 = _mm_loadu_pd(ct[l]+t2+28);\
			mc15 = _mm_loadu_pd(ct[l]+t2+30);\
			mc16 = _mm_loadu_pd(ct[l]+t2+32);\
			mc17 = _mm_loadu_pd(ct[l]+t2+34);\
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

PetscInt BSG_MatMult_6(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscScalar *xt[7];
	__m128d mx0, mx1, mx2, msum0, msum1, msum2, msum3, msum4,msum5, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, mc8, mc9,  mc10, mc11, mc12, mc13, mc14, mc15, mc16, mc17;
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_6();
                for(l=mnos;l<nos;l++)
                {
			inline_6();
                }
		save_6();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_6();
		for(l=mnos-1;l<nos;l++)
		{
			inline_6();
		}
		save_6();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_6();
		for(l=mnos-2;l<nos;l++)
		{
			inline_6();
		}
		save_6();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_6();
		for(l=0;l<nos;l++)
		{
			inline_6();
		}
		save_6();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_6();
		for(l=0;l<nos-1;l++)
		{
			inline_6();
		}
		save_6();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_6();
		for(l=0;l<nos-2;l++)
		{
			inline_6();
		}
		save_6();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_6();
		for(l=0;l<=mnos;l++)
		{
			inline_6();
		}
		save_6();
	}
	PetscFunctionReturn(0);
}

#define inline_Neven_old()  for(i=0;i<dofby2;i++){\
				mx[i] = _mm_loadu_pd(xt[l]+t1+2*i);\
			}\
			PetscPrefetchBlock(xt[l+1]+t1,dof,0,PETSC_PREFETCH_HINT_NTA);\
			\
			for(j=0;j<dof*dofby2;j++){\
				mc[j] = _mm_loadu_pd(ct[l]+t2+2*j);\
			}\
			PetscPrefetchBlock(ct[l+1]+t2,bs,0,PETSC_PREFETCH_HINT_NTA);\
			\
			for(i=0;i<dof;i+=2){\
				t3 = i*dofby2;\
				t4 = (i+1)*dofby2;\
				for(j=0;j<dofby2;j++){\
					msum[i] = _mm_add_pd(msum[i], _mm_mul_pd(mx[j],mc[t3+j]));\
					msum[i+1] = _mm_add_pd(msum[i+1], _mm_mul_pd(mx[j],mc[t4+j]));\
				}\
			}

#define setup_Neven_old()   t1= k*dof; t2 = t1*dof;\
			for(i=0;i<dof;i+=2){\
				msum[i] =_mm_set_pd(0,0);\
				msum[i+1] =_mm_set_pd(0,0);\
			}

#define save_Neven_old()    for(i=0;i<dof;i+=2){\
				msum[i] = _mm_hadd_pd(msum[i],msum[i+1]);\
				_mm_storeu_pd(y+t1+i, msum[i]);\
			}


PetscInt BSG_MatMult_Neven_old(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt i,j,k,l, t1, t2, t3, t4;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = dof/2;
	const PetscScalar *xt[7];
	__m128d mx[dof/2], msum[dof], mc[dof*dof/2];
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_Neven_old();
                for(l=mnos;l<nos;l++)
                {
			inline_Neven_old();
                }
		save_Neven_old();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_Neven_old();
		for(l=mnos-1;l<nos;l++)
		{
			inline_Neven_old();
		}
		save_Neven_old();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_Neven_old();
		for(l=mnos-2;l<nos;l++)
		{
			inline_Neven_old();
		}
		save_Neven_old();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_Neven_old();
		for(l=0;l<nos;l++)
		{
			inline_Neven_old();
		}
		save_Neven_old();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_Neven_old();
		for(l=0;l<nos-1;l++)
		{
			inline_Neven_old();
		}
		save_Neven_old();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_Neven_old();
		for(l=0;l<nos-2;l++)
		{
			inline_Neven_old();
		}
		save_Neven_old();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_Neven_old();
		for(l=0;l<=mnos;l++)
		{
			inline_Neven_old();
		}
		save_Neven_old();
	}
	PetscFunctionReturn(0);
}

#define setup_Neven()	for(i=0;i<dofby2;i++){\
				 msum[i] = _mm_set_pd(0,0);\
			 }\

#define save_Neven() for(i=0;i<dofby2;i++){\
				_mm_storeu_pd(y+k*dof+2*i,msum[i]);\
			 }\

#define setup12_Neven()	t1= k*dof+l2;\
			mx0 = _mm_loadu_pd(xt[l]+t1);\
                        mx1 = _mm_loadu_pd(xt[l]+t1+2)

#define setup34_Neven()	t1= k*dof+l2;\
			mx0 = _mm_loadu_pd(xt[l]+t1)

#define inline_stage1_Neven()	t2 = k*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc1 = _mm_loadu_pd(ct[l]+t2+2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof);\
                        mc3 = _mm_loadu_pd(ct[l]+t2+dof+2);\
                        mc4 = _mm_loadu_pd(ct[l]+t2+2*dof);\
                        mc5 = _mm_loadu_pd(ct[l]+t2+2*dof+2);\
                        mc6 = _mm_loadu_pd(ct[l]+t2+3*dof);\
                        mc7 = _mm_loadu_pd(ct[l]+t2+3*dof+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        mc4 = _mm_add_pd(_mm_mul_pd(mx0,mc4),_mm_mul_pd(mx1,mc5));\
                        mc6 = _mm_add_pd(_mm_mul_pd(mx0,mc6),_mm_mul_pd(mx1,mc7));\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
			msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))


#define inline_stage3_Neven()	t2 = k*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof);\
                        mc4 = _mm_loadu_pd(ct[l]+t2+2*dof);\
                        mc6 = _mm_loadu_pd(ct[l]+t2+3*dof);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        mc4 = _mm_mul_pd(mx0,mc4);\
                        mc6 = _mm_mul_pd(mx0,mc6);\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
			msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))

#define inline_stage2_Neven()	t2 = k*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc1 = _mm_loadu_pd(ct[l]+t2+2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof);\
                        mc3 = _mm_loadu_pd(ct[l]+t2+dof+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage4_Neven()	t2 = k*bs+2*l1*dof+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
			msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

PetscInt BSG_MatMult_Neven(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt i,j,k,l, t1, t2, t3, t4, l1,l2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = dof/2;
	const PetscScalar *xt[7];
	__m128d mx0, mx1, msum[dofby2], mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_Neven();
                for(l=mnos;l<nos;l++)
                {
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven();
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven();
				}
			}
                }
		save_Neven();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_Neven();
		for(l=mnos-1;l<nos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven();
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven();
				}
			}
		}
		save_Neven();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_Neven();
		for(l=mnos-2;l<nos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven();
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven();
				}
			}
		}
		save_Neven();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_Neven();
		for(l=0;l<nos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven();
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven();
				}
			}
		}
		save_Neven();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_Neven();
		for(l=0;l<nos-1;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven();
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven();
				}
			}
		}
		save_Neven();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_Neven();
		for(l=0;l<nos-2;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven();
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven();
				}
			}
		}
		save_Neven();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_Neven();
		for(l=0;l<=mnos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-2; l2 += 4)
			{
				setup12_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage1_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage2_Neven();
				}
			}
			for(; l2 < dof; l2 += 2)
			{
				setup34_Neven();
				for(l1=0; l1<dofby2-1; l1+= 2)
				{
					inline_stage3_Neven();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage4_Neven();
				}
			}
		}
		save_Neven();
	}
	PetscFunctionReturn(0);
}

#define inline_1() msum0 += xt[l][t1] * ct[l][t2]

#define setup_1() t1= k*dof; t2 = t1*dof;\
		msum0 = 0.0

#define save_1() y[t1] = msum0

PetscInt BSG_MatMult_1(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscScalar *xt[7];
	PetscScalar msum0;
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_1();
                for(l=mnos;l<nos;l++)
                {
			inline_1();
                }
		save_1();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_1();
		for(l=mnos-1;l<nos;l++)
		{
			inline_1();
		}
		save_1();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_1();
		for(l=mnos-2;l<nos;l++)
		{
			inline_1();
		}
		save_1();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_1();
		for(l=0;l<nos;l++)
		{
			inline_1();
		}
		save_1();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_1();
		for(l=0;l<nos-1;l++)
		{
			inline_1();
		}
		save_1();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_1();
		for(l=0;l<nos-2;l++)
		{
			inline_1();
		}
		save_1();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_1();
		for(l=0;l<=mnos;l++)
		{
			inline_1();
		}
		save_1();
	}
	PetscFunctionReturn(0);
}

#define inline_3() 	mx0 = _mm_loadu_pd(xt[l]+t1);\
			mx1 = _mm_load1_pd(xt[l]+t1+2);\
			mc0 = _mm_loadu_pd(ct[l]+t2);\
			mc1 = _mm_loadu_pd(ct[l]+t2+2);\
			mc2 = _mm_loadu_pd(ct[l]+t2+4);\
			mc3 = _mm_loadu_pd(ct[l]+t2+6);\
			mc4 = _mm_loadu_pd(ct[l]+t2+8);\
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

PetscInt BSG_MatMult_3(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscScalar *xt[7];
	__m128d mx0, mx1, msum0, msum1, msum2,msum3, msum4, mc0, mc1, mc2, mc3, mc4;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_3();
                for(l=mnos;l<nos;l++)
                {
			inline_3();
                }
		save_3();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_3();
		for(l=mnos-1;l<nos;l++)
		{
			inline_3();
		}
		save_3();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_3();
		for(l=mnos-2;l<nos;l++)
		{
			inline_3();
		}
		save_3();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_3();
		for(l=0;l<nos;l++)
		{
			inline_3();
		}
		save_3();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_3();
		for(l=0;l<nos-1;l++)
		{
			inline_3();
		}
		save_3();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_3();
		for(l=0;l<nos-2;l++)
		{
			inline_3();
		}
		save_3();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_3();
		for(l=0;l<=mnos;l++)
		{
			inline_3();
		}
		save_3();
	}
	PetscFunctionReturn(0);
}

#define inline_stage1_5()  mx0 = _mm_loadu_pd(xt[l]+t1);\
		 	mx1 = _mm_loadu_pd(xt[l]+t1+2);\
			mx2 = _mm_load1_pd(xt[l]+t1+4);\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc1 = _mm_loadu_pd(ct[l]+t2+2);\
			mc2 = _mm_loadu_pd(ct[l]+t2+4);\
			mc3 = _mm_loadu_pd(ct[l]+t2+6);\
			mc4 = _mm_loadu_pd(ct[l]+t2+20);\
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

#define inline_stage2_5() 	mc0 = _mm_loadu_pd(ct[l]+t2+8);\
			mc1 = _mm_loadu_pd(ct[l]+t2+10);\
			mc2 = _mm_loadu_pd(ct[l]+t2+12);\
			mc3 = _mm_loadu_pd(ct[l]+t2+14);\
			mc4 = _mm_loadu_pd(ct[l]+t2+16);\
			mc5 = _mm_loadu_pd(ct[l]+t2+18);\
			mc6 = _mm_loadu_pd(ct[l]+t2+22);\
			mc7 = _mm_loadu_pd(ct[l]+t2+24);\
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

PetscInt BSG_MatMult_5(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscScalar *xt[7];
	__m128d mx0, mx1,mx2, msum0, msum1, msum2,msum3, msum4, mc0, mc1, mc2, mc3, mc4,mc5, mc6, mc7;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	for(k = 0; k < 1; k++)
	{
		setup_stage1_5();
                for(l=mnos;l<nos;l++)
                {
			inline_stage1_5();
                }
		save_stage1_5();
		setup_stage2_5();
                for(l=mnos;l<nos;l++)
                {
			inline_stage2_5();
                }
		save_stage2_5();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_stage1_5();
                for(l=mnos-1;l<nos;l++)
                {
			inline_stage1_5();
                }
		save_stage1_5();
		setup_stage2_5();
		for(l=mnos-1;l<nos;l++)
		{
			inline_stage2_5();
		}
		save_stage2_5();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_stage1_5();
                for(l=mnos-2;l<nos;l++)
                {
			inline_stage1_5();
                }
		save_stage1_5();
		setup_stage2_5();
		for(l=mnos-2;l<nos;l++)
		{
			inline_stage2_5();
		}
		save_stage2_5();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_stage1_5();
                for(l=0;l<nos;l++)
                {
			inline_stage1_5();
                }
		save_stage1_5();
		setup_stage2_5();
		for(l=0;l<nos;l++)
		{
			inline_stage2_5();
		}
		save_stage2_5();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_stage1_5();
                for(l=0;l<nos-1;l++)
                {
			inline_stage1_5();
                }
		save_stage1_5();
		setup_stage2_5();
		for(l=0;l<nos-1;l++)
		{
			inline_stage2_5();
		}
		save_stage2_5();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_stage1_5();
                for(l=0;l<nos-2;l++)
                {
			inline_stage1_5();
                }
		save_stage1_5();
		setup_stage2_5();
		for(l=0;l<nos-2;l++)
		{
			inline_stage2_5();
		}
		save_stage2_5();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_stage1_5();
                for(l=0;l<=mnos;l++)
                {
			inline_stage1_5();
                }
		save_stage1_5();
		setup_stage2_5();
		for(l=0;l<=mnos;l++)
		{
			inline_stage2_5();
		}
		save_stage2_5();
	}
	
	PetscFunctionReturn(0);
}

#define setup_5_ver1() t1= k*dof; t2 = t1*dof;\
		msum0 =_mm_set_pd(0,0);\
		msum1 =_mm_set_pd(0,0);\
		msum2 =_mm_set_pd(0,0);\
		msum3 =_mm_set_pd(0,0);\
		msum4 =_mm_set_pd(0,0);\

#define inline_1ststage_5_ver1() 	mx = _mm_loadu_pd(xt[l]+t1);\
			mc0 = _mm_loadu_pd(ct[l]+t2);\
			mc1 = _mm_loadu_pd(ct[l]+t2+4);\
			mc2 = _mm_loadu_pd(ct[l]+t2+8);\
			mc3 = _mm_loadu_pd(ct[l]+t2+12);\
			mc4 = _mm_loadu_pd(ct[l]+t2+16);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx,mc0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx,mc1));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx,mc2));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx,mc3));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx,mc4));\
		 	mx = _mm_loadu_pd(xt[l]+t1+2);\
			mc0 = _mm_loadu_pd(ct[l]+t2+2);\
			mc1 = _mm_loadu_pd(ct[l]+t2+6);\
			mc2 = _mm_loadu_pd(ct[l]+t2+10);\
			mc3 = _mm_loadu_pd(ct[l]+t2+14);\
			mc4 = _mm_loadu_pd(ct[l]+t2+18);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx,mc0));\
			msum1 = _mm_add_pd(msum1, _mm_mul_pd(mx,mc1));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx,mc2));\
			msum3 = _mm_add_pd(msum3, _mm_mul_pd(mx,mc3));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx,mc4));\

#define save_5_ver1() msum0 = _mm_hadd_pd(msum0,msum1);\
		msum2 = _mm_hadd_pd(msum2,msum3);\
		msum4 = _mm_hadd_pd(msum4,msum4)

#define inline_2ndstage_5_ver1() mx = _mm_load1_pd(xt[l]+t1+4);\
			mc0 = _mm_loadu_pd(ct[l]+t2+20);\
			mc2 = _mm_loadu_pd(ct[l]+t2+22);\
			mc4 = _mm_loadu_pd(ct[l]+t2+24);\
			msum0 = _mm_add_pd(msum0, _mm_mul_pd(mx,mc0));\
			msum2 = _mm_add_pd(msum2, _mm_mul_pd(mx,mc2));\
			msum4 = _mm_add_pd(msum4, _mm_mul_pd(mx,mc4));\
			_mm_storeu_pd(y+t1, msum0);\
			_mm_storeu_pd(y+t1+2, msum2);\
			_mm_maskstore_pd(y+t1+4,xtemp,msum4)

PetscInt BSG_MatMult_5_ver1(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt k,l, t1, t2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscScalar *xt[7];
	__m128d mx, msum0,msum1, msum2,msum3, msum4, mc0, mc1, mc2, mc3, mc4;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	for(k = 0; k < 1; k++)
	{
		setup_5_ver1();
                for(l=mnos;l<nos;l++)
                {
			inline_1ststage_5_ver1();
                }
		save_5_ver1();
                for(l=mnos;l<nos;l++)
                {
			inline_2ndstage_5_ver1();
                }
	}

	for(k = 1; k < lda3; k++)
	{
		setup_5_ver1();
		for(l=mnos-1;l<nos;l++)
		{
			inline_1ststage_5_ver1();
		}
		save_5_ver1();
                for(l=mnos-1;l<nos;l++)
                {
			inline_2ndstage_5_ver1();
                }
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_5_ver1();
		for(l=mnos-2;l<nos;l++)
		{
			inline_1ststage_5_ver1();
		}
		save_5_ver1();
                for(l=mnos-2;l<nos;l++)
                {
			inline_2ndstage_5_ver1();
                }
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_5_ver1();
		for(l=0;l<nos;l++)
		{
			inline_1ststage_5_ver1();
		}
		save_5_ver1();
                for(l=0;l<nos;l++)
                {
			inline_2ndstage_5_ver1();
                }
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_5_ver1();
		for(l=0;l<nos-1;l++)
		{
			inline_1ststage_5_ver1();
		}
		save_5_ver1();
                for(l=0;l<nos-1;l++)
                {
			inline_2ndstage_5_ver1();
                }
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_5_ver1();
		for(l=0;l<nos-2;l++)
		{
			inline_1ststage_5_ver1();
		}
		save_5_ver1();
                for(l=0;l<nos-2;l++)
                {
			inline_2ndstage_5_ver1();
                }
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_5_ver1();
		for(l=0;l<=mnos;l++)
		{
			inline_1ststage_5_ver1();
		}
		save_5_ver1();
                for(l=0;l<=mnos;l++)
                {
			inline_2ndstage_5_ver1();
                }
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

#define setup123_Nodd() t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt[l]+t1);\
                        mx1 = _mm_loadu_pd(xt[l]+t1+2)

#define setup456_Nodd() t1= k*dof+l2;\
                        mx0 = _mm_loadu_pd(xt[l]+t1)

#define setup789_Nodd() t1= k*dof+l2;\
                        mx0 = _mm_load1_pd(xt[l]+t1)

#define inline_stage1_Nodd()   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc1 = _mm_loadu_pd(ct[l]+t2+2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof-1);\
                        mc3 = _mm_loadu_pd(ct[l]+t2+dof+1);\
                        mc4 = _mm_loadu_pd(ct[l]+t2+2*dof-2);\
                        mc5 = _mm_loadu_pd(ct[l]+t2+2*dof);\
                        mc6 = _mm_loadu_pd(ct[l]+t2+3*dof-3);\
                        mc7 = _mm_loadu_pd(ct[l]+t2+3*dof-1);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        mc4 = _mm_add_pd(_mm_mul_pd(mx0,mc4),_mm_mul_pd(mx1,mc5));\
                        mc6 = _mm_add_pd(_mm_mul_pd(mx0,mc6),_mm_mul_pd(mx1,mc7));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))


#define inline_stage4_Nodd()   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof-1);\
                        mc4 = _mm_loadu_pd(ct[l]+t2+2*dof-2);\
                        mc6 = _mm_loadu_pd(ct[l]+t2+3*dof-3);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        mc4 = _mm_mul_pd(mx0,mc4);\
                        mc6 = _mm_mul_pd(mx0,mc6);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2));\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], _mm_hadd_pd(mc4,mc6))

#define inline_stage7_Nodd() t2 = k*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc1 = _mm_loadu_pd(ct[l]+t2+2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0);\
                        msum[l1+1] = _mm_add_pd(msum[l1+1], mc1)

#define inline_stage2_Nodd()   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc1 = _mm_loadu_pd(ct[l]+t2+2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof-1);\
                        mc3 = _mm_loadu_pd(ct[l]+t2+dof+1);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        mc2 = _mm_add_pd(_mm_mul_pd(mx0,mc2),_mm_mul_pd(mx1,mc3));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage5_Nodd()   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc2 = _mm_loadu_pd(ct[l]+t2+dof-1);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        mc2 = _mm_mul_pd(mx0,mc2);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc2))

#define inline_stage8_Nodd() t2 = k*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0)

#define inline_stage3_Nodd()   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc1 = _mm_loadu_pd(ct[l]+t2+2);\
                        mc0 = _mm_add_pd(_mm_mul_pd(mx0,mc0),_mm_mul_pd(mx1,mc1));\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define inline_stage6_Nodd()   t2 = k*bs+2*l1*(dof-1)+l2;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        mc0 = _mm_mul_pd(mx0,mc0);\
                        msum[l1] = _mm_add_pd(msum[l1], _mm_hadd_pd(mc0,mc0))

#define inline_stage9_Nodd() t2 = k*bs+dof*(dof-1)+2*l1;\
                        mc0 = _mm_loadu_pd(ct[l]+t2);\
                        msum[l1] = _mm_add_pd(msum[l1], mc0);\

PetscInt BSG_MatMult_Nodd(PetscScalar ** ct,const PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt i,j,k,l, t1, t2, t3, t4, l1,l2;
	const PetscInt lda3 = m;
	const PetscInt lda2 = lda3 * n;
	const PetscInt lda1 = lda2 * p;
	const PetscInt mnos = 3;
	const PetscInt dofby2 = (dof+1)/2;
	const PetscScalar *xt[7];
	__m128d mx0, mx1, msum[dofby2], mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7;
	__m128i xtemp = _mm_set_epi32(0,0,-1,-1);
	for(l=0;l<7;l++)
		xt[l] = x + (idx[l] + idy[l]*lda3 + idz[l]*lda2)*dof;
	
	for(k = 0; k < 1; k++)
	{
		setup_Nodd();
                for(l=mnos;l<nos;l++)
                {
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd();
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd();
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd();
				}
			}
                }
		save_Nodd();
	}

	for(k = 1; k < lda3; k++)
	{
		setup_Nodd();
		for(l=mnos-1;l<nos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd();
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd();
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd();
				}
			}
		}
		save_Nodd();
	}

	for(k = lda3; k < lda2; k++)
	{
		setup_Nodd();
		for(l=mnos-2;l<nos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd();
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd();
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd();
				}
			}
		}
		save_Nodd();
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		setup_Nodd();
		for(l=0;l<nos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd();
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd();
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd();
				}
			}
		}
		save_Nodd();
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		setup_Nodd();
		for(l=0;l<nos-1;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd();
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd();
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd();
				}
			}
		}
		save_Nodd();
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		setup_Nodd();
		for(l=0;l<nos-2;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd();
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd();
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd();
				}
			}
		}
		save_Nodd();
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
		setup_Nodd();
		for(l=0;l<=mnos;l++)
		{
			PetscPrefetchBlock(xt[l+1]+k*dof,dof,0,PETSC_PREFETCH_HINT_NTA);
			PetscPrefetchBlock(ct[l+1]+k*bs,bs,0,PETSC_PREFETCH_HINT_NTA);
			for(l2 = 0; l2 < dof-3; l2 += 4)
			{
				setup123_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage1_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage2_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage3_Nodd();
				}
			}
			for(; l2 < dof-1; l2 += 2)
			{
				setup456_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage4_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage5_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage6_Nodd();
				}
			}
			for(; l2 < dof; l2 ++)
			{
				setup789_Nodd();
				for(l1=0; l1<dofby2-2; l1+= 2)
				{
					inline_stage7_Nodd();
				}
				for(;l1<dofby2-1;l1++)
				{
					inline_stage8_Nodd();
				}
				for(;l1<dofby2;l1++)
				{
					inline_stage9_Nodd();
				}
			}
		}
		save_Nodd();
	}
	PetscFunctionReturn(0);
}
