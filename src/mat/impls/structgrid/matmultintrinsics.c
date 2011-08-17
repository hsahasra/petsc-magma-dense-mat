#include <string.h>
#include<omp.h>
#define NUM_THREADS 4
int OPENMP;



#include "../src/mat/impls/structgrid/matstructgrid.h"

/*  -------------------------------------------------------------------- 
     This file implements matrix multiplication for the structgrid data type. 
     The routine employs SSE/AVX intrinsics if they are available on the machine.
     Otherwise, the computations default to normal PetscScalar operations. 
     The instruction for fused addmultiply has not been implemented of date.

     Author: Chekuri S. Choudary, RNET
*/



/** transparent short vector at compile time**/

#ifdef __AVX__ //Use 256 AVX intrinsics
#include <immintrin.h>
#define SV_DOUBLE_WIDTH 4
#define __sv_dtype __m256d
#define _sv_loadu_pd(a)  _mm256_loadu_pd(a)
#define _sv_add_pd(a,b) _mm256_add_pd(a,b)
#define _sv_mul_pd(a,b) _mm256_mul_pd(a,b)
#define _sv_storeu_pd(a,b) _mm256_storeu_pd(a,b)
#elif defined(__SSE2__) //Use 128 bit SSE intrinsics
#include <emmintrin.h>
#define SV_DOUBLE_WIDTH 2
#define __sv_dtype __m128d
#define _sv_loadu_pd(a)  _mm_loadu_pd(a)
#define _sv_add_pd(a,b) _mm_add_pd(a,b)
#define _sv_mul_pd(a,b) _mm_mul_pd(a,b)
#define _sv_storeu_pd(a,b) _mm_storeu_pd(a,b)
#else
#define SV_DOUBLE_WIDTH 1
#define __sv_dtype PetscScalar
#define _sv_loadu_pd(a) (*(a))
#define _sv_add_pd(a,b) ((a)+(b))
#define _sv_mul_pd(a,b) ((a)*(b))
#define _sv_storeu_pd(a,b) (*(a)=(b))
#endif


PetscInt SG_MatMult(PetscScalar * coeff, PetscScalar * xi, PetscScalar * y,PetscScalar * x, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos )
{
  //printf("Start of SG_MatMult\n");
	PetscInt i,j,k,l,xdisp,ydisp,zdisp;
	PetscInt lda1 = m*n*p*dof;
	PetscInt lda2 = m*n*dof;
	PetscInt lda3 = m*dof;
	__sv_dtype yv, xv, coeffv,xv1,coeffv1;
#ifdef __AVX__
	__sv_dtype xv2, coeffv2, xv3, coeffv3;
#endif

	PetscInt xval[nos], offset[nos];
	memcpy(x+lda2,xi,sizeof(PetscScalar)*lda1);
	for(l=0;l<nos;l++)
	{
		xdisp = idx[l]; ydisp = idy[l] ; zdisp = idz[l]; offset[l] = l*lda1;
                //printf("offset[%d]: %d\n",l,offset[l]);
	 	xval[l] = xdisp + ydisp*lda3 + zdisp*lda2;
	}
if(OPENMP){
	#pragma omp parallel private(coeff,coeffv1,coeffv2,coeffv3,x,xv1,xv2,xv3) shared(yv)
	#pragma omp for schedule(static)
}	for(k=0;(k+SV_DOUBLE_WIDTH)<lda1;k+=SV_DOUBLE_WIDTH)
	{
		yv = _sv_loadu_pd((PetscScalar *)(y+k));
		for(l=0;(l+SV_DOUBLE_WIDTH)<nos;l+=SV_DOUBLE_WIDTH)
		{
			xv = _sv_loadu_pd((PetscScalar *)(x+lda2+xval[l]+k));

#if (SV_DOUBLE_WIDTH > 1)
			xv1 = _sv_loadu_pd((PetscScalar *)(x+lda2+xval[l+1]+k));
#endif

#if (SV_DOUBLE_WIDTH > 2)
			xv2 = _sv_loadu_pd((PetscScalar *)(x+lda2+xval[l+2]+k));
			xv3 = _sv_loadu_pd((PetscScalar *)(x+lda2+xval[l+3]+k));
#endif
			
			coeffv = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k));

#if (SV_DOUBLE_WIDTH > 1)
			coeffv1 = _sv_loadu_pd((PetscScalar *)(coeff+offset[l+1]+k));
#endif

#if (SV_DOUBLE_WIDTH > 2) 
			coeffv2 = _sv_loadu_pd((PetscScalar *)(coeff+offset[l+2]+k));
			coeffv3 = _sv_loadu_pd((PetscScalar *)(coeff+offset[l+3]+k));
#endif

			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv,xv));

#if (SV_DOUBLE_WIDTH > 1)
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv1,xv1));
#endif

#if (SV_DOUBLE_WIDTH > 2)
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv2,xv2));
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv3,xv3));
#endif

		}
		for(;l<nos;l++)
		{
			xv = _sv_loadu_pd((PetscScalar *)(x+lda2+xval[l]+k));
			coeffv = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k));
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv,xv));
		}
		_sv_storeu_pd((PetscScalar *)(y+k),yv);	
	}
	for(;k<lda1;k++)
		for(l=0;l<nos;l++)
	{
			y[k] += (coeff[offset[l]+k] * x[lda2+(xval[l]+k)]);		
	}

	PetscFunctionReturn(0);
}


