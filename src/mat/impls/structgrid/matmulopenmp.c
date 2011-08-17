#include "../src/mat/impls/structgrid/matstructgrid.h"
#include<omp.h>
#define NUM_THREADS 2
#include<stdio.h>
#include<string.h>

PetscInt SG_MatMultOpenmp(PetscScalar * coeff, PetscScalar * xi, PetscScalar * y,PetscScalar * x, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos )
{
	//printf("In OpenMP Matmult\n");
//	fflush(stdout);

	PetscInt k,l;
	PetscInt lda1 = m*n*p*dof;
	PetscInt lda2 = m*n*dof;
	PetscInt lda3 = m*dof;

	PetscInt xval[nos], offset[nos];
	memcpy(x+lda2,xi,sizeof(PetscScalar)*lda1);
	for(l=0;l<nos;l++)
	{
		offset[l] = l*lda1;
	 	xval[l] = idx[l] + idy[l]*lda3 + idz[l]*lda2;
	}
	
	#pragma omp parallel
	#pragma omp for schedule(static)
	for(k=0;k<lda1;k++)
		for(l=0;l<nos;l++)
	{
			y[k] += (coeff[offset[l]+k] * x[lda2+(xval[l]+k)]);		
	}

	PetscFunctionReturn(0);
}

