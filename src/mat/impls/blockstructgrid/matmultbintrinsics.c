#include <string.h>

//#define SPREFETCH



#include "../src/mat/impls/blockstructgrid/matblockstructgrid.h"

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
#define _sv_set_pd(a,b,c,d) _mm256_set_pd(a,b,c,d)
#elif defined(__SSE2__) //Use 128 bit SSE intrinsics
#include <emmintrin.h>
#define SV_DOUBLE_WIDTH 2
#define __sv_dtype __m128d
#define _sv_loadu_pd(a)  _mm_loadu_pd(a)
#define _sv_add_pd(a,b) _mm_add_pd(a,b)
#define _sv_mul_pd(a,b) _mm_mul_pd(a,b)
#define _sv_storeu_pd(a,b) _mm_storeu_pd(a,b)
#define _sv_set_pd(a,b,c,d) _mm_set_pd(c,d)
#else
#define SV_DOUBLE_WIDTH 1
#define __sv_dtype PetscScalar
#define _sv_loadu_pd(a) (*(a))
#define _sv_add_pd(a,b) ((a)+(b))
#define _sv_mul_pd(a,b) ((a)*(b))
#define _sv_storeu_pd(a,b) (*(a)=(b))
#define _sv_set_pd(a,b,c,d) (d)
#endif


//MatMult Without Padding, with openmp, with software prefetching

PetscInt BSG_MatMult(PetscScalar * coeff, PetscScalar * x, PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt * negRange, PetscInt *coeffOffset, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos, PetscInt dim , PetscInt bs)
{
	PetscInt i,j,k,l,mt;
	PetscInt lda1 = m*n*p;
	PetscInt lda2 = m*n;
	PetscInt lda3 = m;	
	PetscInt xshift[nos];
	PetscInt mnos = dim + 1;
	register PetscInt tempnegRange, tempOffset, tempshift;
	PetscScalar sttemp[SV_DOUBLE_WIDTH];
	//printf("OPENMP=%d\n",OPENMP);
       	//printf("Thread=%d\n",omp_get_thread_num());
	
	__sv_dtype yv[dof], xv, cv, mask_xv;
	PetscInt onexdiff = dof % SV_DOUBLE_WIDTH ;
	if(onexdiff == 0)
		mask_xv = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH > 1)
	else if (onexdiff == 1)
		mask_xv = _sv_set_pd(0,0,0,1);
#endif
#if (SV_DOUBLE_WIDTH > 2)
	else if (onexdiff == 2)
		mask_xv = _sv_set_pd(0,0,1,1);
	else
		mask_xv = _sv_set_pd(0,1,1,1);
#endif
	for(j=0; j< nos; j++)
	{
		xshift[j] = idx[j] + idy[j]*lda3 + idz[j] * lda2;
	}
	for(j=0; j< dof; j++)
		yv[j] = _sv_set_pd(0,0,0,0);
	for(k = 0; k < 1; k++)
	{
		for(l=0;l<mnos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for (j=0; j< dof; j++)
		{
			_sv_storeu_pd(&sttemp[0], yv[j]);
			yv[j] = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH == 1)
				y[k*dof+j] += sttemp[0];
#elif (SV_DOUBLE_WIDTH == 2)
				y[k*dof+j] += sttemp[0] + sttemp[1];
#else
				y[k*dof+j] += sttemp[0] + sttemp[1] + sttemp[2] + sttemp[3];
#endif
		}
	}

	for(k = 1; k < lda3; k++)
	{
		for(l=0;l<mnos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
//		for(l=mnos;l<mnos+1;l++)
//		{
			l=mnos;
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			tempnegRange = negRange[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

//		}
		for (j=0; j< dof; j++)
		{
			_sv_storeu_pd(&sttemp[0], yv[j]);
			yv[j] = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH == 1)
				y[k*dof+j] += sttemp[0];
#elif (SV_DOUBLE_WIDTH == 2)
				y[k*dof+j] += sttemp[0] + sttemp[1];
#else
				y[k*dof+j] += sttemp[0] + sttemp[1] + sttemp[2] + sttemp[3];
#endif
		}
	}

	for(k = lda3; k < lda2; k++)
	{
		for(l=0;l<mnos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for(l=mnos;l<mnos+2;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			tempnegRange = negRange[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for (j=0; j< dof; j++)
		{
			_sv_storeu_pd(&sttemp[0], yv[j]);
			yv[j] = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH == 1)
				y[k*dof+j] += sttemp[0];
#elif (SV_DOUBLE_WIDTH == 2)
				y[k*dof+j] += sttemp[0] + sttemp[1];
#else
				y[k*dof+j] += sttemp[0] + sttemp[1] + sttemp[2] + sttemp[3];
#endif
		}
	}

	for(k = lda2; k < (lda1- lda2); k++)
	{
		for(l=0;l<mnos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for(l=mnos;l<nos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			tempnegRange = negRange[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for (j=0; j< dof; j++)
		{
			_sv_storeu_pd(&sttemp[0], yv[j]);
			yv[j] = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH == 1)
				y[k*dof+j] += sttemp[0];
#elif (SV_DOUBLE_WIDTH == 2)
				y[k*dof+j] += sttemp[0] + sttemp[1];
#else
				y[k*dof+j] += sttemp[0] + sttemp[1] + sttemp[2] + sttemp[3];
#endif
		}
	}

	for(k = (lda1 - lda2); k < (lda1 - lda3); k++)
	{
		for(l=0;l<mnos-1;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for(l=mnos;l<nos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			tempnegRange = negRange[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for (j=0; j< dof; j++)
		{
			_sv_storeu_pd(&sttemp[0], yv[j]);
			yv[j] = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH == 1)
				y[k*dof+j] += sttemp[0];
#elif (SV_DOUBLE_WIDTH == 2)
				y[k*dof+j] += sttemp[0] + sttemp[1];
#else
				y[k*dof+j] += sttemp[0] + sttemp[1] + sttemp[2] + sttemp[3];
#endif
		}
	}

	for(k = (lda1 - lda3); k < (lda1 - 1); k++)
	{
		for(l=0;l<mnos-2;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for(l=mnos;l<nos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			tempnegRange = negRange[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for (j=0; j< dof; j++)
		{
			_sv_storeu_pd(&sttemp[0], yv[j]);
			yv[j] = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH == 1)
				y[k*dof+j] += sttemp[0];
#elif (SV_DOUBLE_WIDTH == 2)
				y[k*dof+j] += sttemp[0] + sttemp[1];
#else
				y[k*dof+j] += sttemp[0] + sttemp[1] + sttemp[2] + sttemp[3];
#endif
		}
	}

	for(k = (lda1 - 1); k < (lda1); k++)
	{
//		for(l=0;l<mnos-3;l++)
//		{
			l=0;
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

//		}
		for(l=mnos;l<nos;l++)
		{
			tempOffset = coeffOffset[l];
			tempshift = xshift[l];
			tempnegRange = negRange[l];
			for(i=0; (i+SV_DOUBLE_WIDTH) <= dof; i+=SV_DOUBLE_WIDTH)
			{
				xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
				for (j=0; j<dof; j++)
				{
					cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
					yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
				}
			}
			xv = _sv_loadu_pd(x+(k+tempshift)*dof + i);
			xv = _sv_mul_pd(xv,mask_xv);
			for (j=0; j<dof; j++)
			{
				cv = _sv_loadu_pd (coeff+tempOffset*bs + k*dof*dof+ j*dof+ i-tempnegRange);
				yv[j] = _sv_add_pd(yv[j], _sv_mul_pd(xv,cv));
			}

		}
		for (j=0; j< dof; j++)
		{
			_sv_storeu_pd(&sttemp[0], yv[j]);
			yv[j] = _sv_set_pd(0,0,0,0);
#if (SV_DOUBLE_WIDTH == 1)
				y[k*dof+j] += sttemp[0];
#elif (SV_DOUBLE_WIDTH == 2)
				y[k*dof+j] += sttemp[0] + sttemp[1];
#else
				y[k*dof+j] += sttemp[0] + sttemp[1] + sttemp[2] + sttemp[3];
#endif
		}
	}

	PetscFunctionReturn(0);
}
