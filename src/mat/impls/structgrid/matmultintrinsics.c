#include <string.h>

#include<omp.h>
int OPENMP=0;
//#define SPREFETCH



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

/*

//MatMult With Padding

PetscInt SG_MatMult(PetscScalar * coeff, PetscScalar * xi, PetscScalar * y,PetscScalar * x, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos )
{
	PetscInt i,j,k,l,xdisp,ydisp,zdisp;
	PetscInt lda1 = m*n*p*dof;
	PetscInt lda2 = m*n*dof;
	PetscInt lda3 = m*dof;

	PetscInt xval[nos], offset[nos];
	memcpy(x+lda2,xi,sizeof(PetscScalar)*lda1);
	for(l=0;l<nos;l++)
	{
		xdisp = idx[l]; ydisp = idy[l] ; zdisp = idz[l]; offset[l] = l*lda1;
	 	xval[l] = xdisp + ydisp*lda3 + zdisp*lda2;
	}

#pragma omp parallel if(OPENMP) firstprivate(lda1,lda2,xval,offset,nos,x,coeff) shared(y) default(none)
{
       	//printf("Thread=%d\n",omp_get_thread_num(),k,l);
	
	__sv_dtype yv, xv, coeffv,xv1,coeffv1;
#ifdef __AVX__
	__sv_dtype xv2, coeffv2, xv3, coeffv3;
#endif

	#pragma omp for nowait private(l,xv,coeffv,xv1,coeffv1,yv) 
	for(k=0;k<=(lda1-SV_DOUBLE_WIDTH);k+=SV_DOUBLE_WIDTH)
	{
		yv = _sv_loadu_pd((PetscScalar *)(y+k));
		for(l=0;l<=(nos-SV_DOUBLE_WIDTH);l+=SV_DOUBLE_WIDTH)
		{
		//	printf("Tiled thread=%d k=%d l=%d\n",omp_get_thread_num(),k,l);
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
		//	printf("Tiled rest l:thread=%d k=%d l=%d\n",omp_get_thread_num(),k,l);
			xv = _sv_loadu_pd((PetscScalar *)(x+lda2+xval[l]+k));
			coeffv = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k));
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv,xv));
		}
		//#pragma omp critical
		//{
		_sv_storeu_pd((PetscScalar *)(y+k),yv);	
		//}
	}
}
//	#pragma omp for
	for(k=(lda1-(lda1%SV_DOUBLE_WIDTH));k<lda1;k++){
		for(l=0;l<nos;l++)
		{
		//	printf("Rest: thread=%d k=%d l=%d\n",omp_get_thread_num(),k,l);
			y[k] += (coeff[offset[l]+k] * x[lda2+(xval[l]+k)]);		
		}
	}

	PetscFunctionReturn(0);
}
*/

//MatMult Without Padding, with openmp, with software prefetching

PetscInt SG_MatMult(PetscScalar * coeff, PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim )
{
	PetscInt i,j,k,l,xdisp,ydisp,zdisp;
	PetscInt lda1 = m*n*p*dof;
	PetscInt lda2  = 0;
	if(dim == 3)
	 lda2 = m*n*dof;
	else
	 lda2 = m*dof;
	PetscInt lda3 = m*dof;
	
	PetscInt _smallval, _largeval, _startval;
	if((lda2+(dof-1)) < (lda1-lda2-(dof-1)))
	{
		_smallval = (lda2+(dof-1));
		_largeval = (lda1-lda2-(dof-1));
	}
	else
	{
		_largeval = (lda2+(dof-1));
		_smallval = (lda1-lda2-(dof-1));
	}

	PetscInt xval[nos], offset[nos], vbeg[nos], vend[nos];
	for(l=0;l<nos;l++)
	{
		xdisp = idx[l]; ydisp = idy[l] ; zdisp = idz[l]; offset[l] = l*lda1;
	 	xval[l] = xdisp + ydisp*lda3 + zdisp*lda2;
		vbeg[l] = 0; vend[l] = m*dof*n*p; // by default the boundaries are 0 to Nz
		if(xval[l] > 0) vend[l] -= xval[l]; // for the positive stencils vend = Nz- xval
		else
		{
			vbeg[l] -= xval[l]; // for the negative stencils vbeg is +-xval
			xval[l] = 0; // xval, the starting column index is zero for negative stencils. 
		}
	}
	//printf("OPENMP=%d\n",OPENMP);
       	//printf("Thread=%d\n",omp_get_thread_num());
	
	__sv_dtype yv, xv, coeffv,xv1,coeffv1;
#ifdef __AVX__
	__sv_dtype xv2, coeffv2, xv3, coeffv3;
#endif


/* A Part */	
	//#pragma omp for nowait private(i,k) 
	for(l=0;l<nos;l+=2*(2*dof-1))// for the stencils 000 and positive ones(100,010 etc)
	{
		for(i=0;i<(2*dof-1);i++)// for all degrees of freedom
		{
			for(k=vbeg[l+i];k<dof; k++) //vbeg[l+i] = starting y(row) index, few vbegs are non zeros when dof>1
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
/* F Part */
	//#pragma omp for nowait private(i,k) 
	for(l=(2*dof-1);l<nos;l+=2*(2*dof-1)) // for the negative stencils ie. 1, 3, 6 etc 
	{
		for(i=0;i<(2*dof-1);i++) 
		{
			for(k=vbeg[l+i];k<(lda2+(dof-1)); k++)//starting indices vary whereas ending index would be m*n*dof+(dof-1)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
/*   G part */
	//#pragma omp for nowait private(i,k) 
	for(l=0;l<nos;l+=2*(2*dof-1))
	{
		for(i=0;i<(2*dof-1);i++)
		{
			for(k=(lda1-lda2-(dof-1));k<vend[l+i]; k++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
/*   C part */
	//#pragma omp for nowait private(i,k) 
	for(k=_largeval;k<lda1; k++)
	{
		for(l=(2*dof-1);l<nos;l+=2*(2*dof-1))
		{
			for(i=0;i<(2*dof-1);i++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
/*   B part */
	//#pragma omp for nowait private(l,i) 
	for(k=dof;k<_smallval; k++)
	{
		for(l=0;l<nos;l+=2*(2*dof-1))
		{
			for(i=0;i<(2*dof-1);i++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);
			}
		}
	}
//   E part 
/*#pragma omp parallel if(OPENMP) firstprivate(lda1,lda2,xval,vbeg,vend,offset,dof,nos,x,coeff) shared(y) private(yv,xv,coeffv,xv1,coeffv1) default(none)
{

#if (SV_DOUBLE_WIDTH > 2)
	#pragma omp for nowait private(l,xv,coeffv,xv1,coeffv1,yv,xv2,xv3,coeffv2,coeffv3) 
#elif (SV_DOUBLE_WIDTH > 1)
	#pragma omp for nowait private(l,xv,coeffv,xv1,coeffv1,yv) 
#endif
*/
#if (SV_DOUBLE_WIDTH > 2)
	#pragma omp parallel if(OPENMP) private(yv,xv,coeffv,xv1,coeffv1,xv2,xv3,coeffv2,coeffv3)
#elif (SV_DOUBLE_WIDTH > 1)
#pragma omp parallel if(OPENMP) private(yv,xv,coeffv,xv1,coeffv1)
#endif
{
#pragma omp for nowait private(l)
	for(k=lda2+(dof-1);k<=lda1-lda2-(dof-1)-SV_DOUBLE_WIDTH; k+=SV_DOUBLE_WIDTH)
	{
       	//printf("Thread=%d\n,k=%d",omp_get_thread_num(),k);
		yv = _sv_loadu_pd((PetscScalar *)(y+k));
		for(l=0;l<=(nos-SV_DOUBLE_WIDTH);l+=SV_DOUBLE_WIDTH)
		{
#ifdef SPREFETCH
		//_mm_prefetch( (coeff+offset[l]+k+80),_MM_HINT_NTA);
    		PetscPrefetchBlock(coeff+offset[l]+k+8,8,0,PETSC_PREFETCH_HINT_NTA);    // Prefetch the next cache line
	#if (SV_DOUBLE_WIDTH > 1)
    		PetscPrefetchBlock(coeff+offset[l+1]+k+8,8,0,PETSC_PREFETCH_HINT_NTA);    
	#endif

	#if (SV_DOUBLE_WIDTH > 2) 
    		PetscPrefetchBlock(coeff+offset[l+2]+k+8,8,0,PETSC_PREFETCH_HINT_NTA);
    		PetscPrefetchBlock(coeff+offset[l+3]+k+8,8,0,PETSC_PREFETCH_HINT_NTA);
	#endif
//(Note:Same cache line being prefetched redundantly to avoid if statements, See MatMultv3 for alternate implementation)
#endif

			xv = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k-vbeg[l]));

#if (SV_DOUBLE_WIDTH > 1)
			xv1 = _sv_loadu_pd((PetscScalar *)(x+xval[l+1]+k-vbeg[l+1]));
#endif


#if (SV_DOUBLE_WIDTH > 2)
			xv2 = _sv_loadu_pd((PetscScalar *)(x+xval[l+2]+k-vbeg[l+2]));
			xv3 = _sv_loadu_pd((PetscScalar *)(x+xval[l+3]+k-vbeg[l+3]));
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
			xv = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k-vbeg[l]));
			coeffv = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k));
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv,xv));
		}
		_sv_storeu_pd((PetscScalar *)(y+k),yv);	
	}
}
	for(k=(lda1-lda2-(dof-1))-((lda1-lda2-(dof-1)-(lda2+(dof-1)))%SV_DOUBLE_WIDTH);k<lda1-lda2-(dof-1);k++){//k has to be initialized for the present platform(probably due to gcc 1.2.4. It can probably be removed when GCC is upgraded.
		for(l=0;l<nos;l++)
		{
			y[k] += (coeff[offset[l]+k] * x[(xval[l]+k-vbeg[l])]);		
		}
	}
	PetscFunctionReturn(0);
}

//MatMult Without Padding, with openmp, with software prefetching( optimized for software prefetching for cache line size=64 bytes)
//unroll inner loop for all the elements in a cache line and prefetch the next cache line just once. Note: Implementation incomplete.
/*
PetscInt SG_MatMultv3(PetscScalar * coeff, PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos )
{
	PetscInt i,j,k,l,xdisp,ydisp,zdisp;
	PetscInt lda1 = m*n*p*dof;
	PetscInt lda2 = m*n*dof;
	PetscInt lda3 = m*dof;

	PetscInt xval[nos], offset[nos], vbeg[nos], vend[nos];
	for(l=0;l<nos;l++)
	{
		xdisp = idx[l]; ydisp = idy[l] ; zdisp = idz[l]; offset[l] = l*lda1;
	 	xval[l] = xdisp + ydisp*lda3 + zdisp*lda2;
		vbeg[l] = 0; vend[l] = m*dof*n*p;
		if(xval[l] > 0) vend[l] -= xval[l];
		else
		{
			vbeg[l] -= xval[l];
			xval[l] = 0;
		}
	}

	//printf("OPENMP=%d\n",OPENMP);

// A Part //	
	//#pragma omp for nowait private(i,k) 
	for(l=0;l<nos;l+=2*(2*dof-1))
	{
		for(i=0;i<(2*dof-1);i++)
		{
			for(k=vbeg[l+i];k<dof; k++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
// F Part //
	//#pragma omp for nowait private(i,k) 
	for(l=(2*dof-1);l<nos;l+=2*(2*dof-1))
	{
		for(i=0;i<(2*dof-1);i++)
		{
			for(k=vbeg[l+i];k<(lda2+(dof-1)); k++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
//   G part //
	//#pragma omp for nowait private(i,k) 
	for(l=0;l<nos;l+=2*(2*dof-1))
	{
		for(i=0;i<(2*dof-1);i++)
		{
			for(k=(lda1-lda2-(dof-1));k<vend[l+i]; k++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
//   D part //
	//#pragma omp for nowait private(i,k) 
	for(l=(2*dof-1);l<nos;l+=2*(2*dof-1))
	{
		for(i=0;i<(2*dof-1);i++)
		{
			for(k=lda1-1-(dof-1);k<vend[l+i]; k++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
//   C part //
	//#pragma omp for nowait private(i,k) 
	for(l=(2*dof-1);l<nos;l+=2*(2*dof-1))
	{
		for(i=0;i<(2*dof-1);i++)
		{
			for(k=start(lda1-lda2-(dof-1),vbeg[l+i]);k<lda1-1-(dof-1); k++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);	
			}
		}
	}
//   B part //
	//#pragma omp for nowait private(l,i) 
	for(k=dof;k<lda2+(dof-1); k++)
	{
		for(l=0;l<nos;l+=2*(2*dof-1))
		{
			for(i=0;i<(2*dof-1);i++)
			{
				y[k] += (coeff[offset[l+i]+k] * x[(xval[l+i])+(k-vbeg[l+i])]);
	
			}
		}
	}
//   E part //
	__sv_dtype yv, xv, coeffv,xv1,coeffv1,yv1;
	__sv_dtype yv2, xv2, coeffv2,xv3,coeffv3,yv3;
#ifdef __AVX__
	__sv_dtype xv2, coeffv2, xv3, coeffv3;
#endif
#pragma omp parallel if(OPENMP) firstprivate(lda1,lda2,xval,vbeg,vend,offset,dof,nos,x,coeff) shared(y) private(yv,xv,coeffv,xv1,coeffv1,yv1,yv2,xv2,coeffv2,yv3,xv3,coeffv3) default(none)
{
       	//printf("Thread=%d\n",omp_get_thread_num());
	
#if (SV_DOUBLE_WIDTH > 2)
	#pragma omp for nowait private(l,xv,coeffv,xv1,coeffv1,yv,yv1,xv2,xv3,coeffv2,coeffv3) 
#endif
#if (SV_DOUBLE_WIDTH > 1)
	#pragma omp for nowait private(l,xv,coeffv,xv1,coeffv1,yv,yv1) 
#endif
	for(k=lda2+(dof-1);k<lda1-lda2-(dof-1)-8; k+=8)
	{
	printf("1st:k=%d\n",k);
		yv = _sv_loadu_pd((PetscScalar *)(y+k));
		yv1 = _sv_loadu_pd((PetscScalar *)(y+k+4));
		for(l=0;l<=(nos-1);l++)
		{
			xv = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k-vbeg[l]));
			coeffv = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k));
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv,xv));
			
			xv1 = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k+2-vbeg[l]));
			coeffv1 = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k+2));
			yv1 = _sv_add_pd(yv1,_sv_mul_pd(coeffv1,xv1));
			
			xv2 = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k+4-vbeg[l]));
			coeffv2 = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k+4));
			yv2 = _sv_add_pd(yv2,_sv_mul_pd(coeffv2,xv2));

			xv3 = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k+6-vbeg[l]));
			coeffv3 = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k+6));
			yv3 = _sv_add_pd(yv3,_sv_mul_pd(coeffv3,xv3));
		}
		_sv_storeu_pd((PetscScalar *)(y+k),yv);	
		_sv_storeu_pd((PetscScalar *)(y+k+2),yv1);	
		_sv_storeu_pd((PetscScalar *)(y+k+4),yv2);	
		_sv_storeu_pd((PetscScalar *)(y+k+6),yv3);	
	}
}
	for(k=(lda1-lda2)-((lda1-lda2)%8);k<lda1-lda2-(dof-1)-SV_DOUBLE_WIDTH; k+=SV_DOUBLE_WIDTH)
	{
	printf("2nd:k=%d\n",k);
		yv = _sv_loadu_pd((PetscScalar *)(y+k));
		for(l=0;l<=(nos-SV_DOUBLE_WIDTH);l+=SV_DOUBLE_WIDTH)
		{

			xv = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k-vbeg[l]));

#if (SV_DOUBLE_WIDTH > 1)
			xv1 = _sv_loadu_pd((PetscScalar *)(x+xval[l+1]+k-vbeg[l+1]));
#endif

#if (SV_DOUBLE_WIDTH > 2)
			xv2 = _sv_loadu_pd((PetscScalar *)(x+xval[l+2]+k-vbeg[l+2]));
			xv3 = _sv_loadu_pd((PetscScalar *)(x+xval[l+3]+k-vbeg[l+3]));
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
			xv = _sv_loadu_pd((PetscScalar *)(x+xval[l]+k-vbeg[l]));
			coeffv = _sv_loadu_pd((PetscScalar *)(coeff+offset[l]+k));
			yv = _sv_add_pd(yv,_sv_mul_pd(coeffv,xv));
		}
		_sv_storeu_pd((PetscScalar *)(y+k),yv);	
	}
	for(;k<lda1-lda2-(dof-1);k++){
	printf("3rd:k=%d\n",k);
		for(l=0;l<nos;l++)
		{
			y[k] += (coeff[offset[l]+k] * x[(xval[l]+k-vbeg[l])]);		
		}
	}
	PetscFunctionReturn(0);
}
*/
