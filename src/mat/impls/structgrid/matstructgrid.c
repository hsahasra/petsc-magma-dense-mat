#define PETSCMAT_DLL


#include "../src/mat/impls/structgrid/matstructgrid.h"
#include "petscblaslapack.h"
#include "petscbt.h"
#include "petscmat.h"
#include "../src/mat/impls/aij/seq/aij.h"  


//#include <immintrin.h>

#include <stdio.h>

static struct _MatOps MatOps_Values = {
/*0*/ MatSetValues_SeqSG,MatGetRow_SeqSG,MatRestoreRow_SeqSG,MatMult_SeqSG,0,
/*5*/0,0,0,0,0,
/*10*/0,0,0,0,0,
/*15*/0,0,MatGetDiagonal_SeqSG,MatDiagonalScale_SeqSG,0,
/*20*/0,0,0,MatZeroEntries_SeqSG,0,
/*25*/0,0,0,0,MatSetUpPreallocation_SeqSG,
/*30*/0,0,0,0,0,
/*35*/0,0,0,0,0,
/*40*/0,0,0,0,0,
/*45*/0,0,0,0,0,
/*50*/0,0,0,0,0,
/*55*/0,0,0,MatSetValuesBlocked_SeqSG,0,
/*60*/MatDestroy_SeqSG,MatView_SeqSG,0,0,0,
/*65*/0,0,MatSetValues_SeqSG,0,MatGetRowMaxAbs_SeqSG,
/*70*/0,0,0,0,0,
/*75*/0,0,0,0,0,
/*80*/0,0,0,0,0,
/*85*/0,0,MatSetValuesBlocked_SeqSG,0,0,
/*90*/0,0,0,0,0,
/*95*/0,0,0,0,0,
/*100*/0,0,0,0,0,
/*105*/0,0,0,0,0,
/*110*/0,0,0,0,0,
/*115*/MatCreate_SeqSG,0,0,0,0,
/*120*/0,0,0,0,0,
/*125*/0,0,0,0,0,
/*130*/0,0,0,0,0,
/*135*/0,0,0,0,MatSetStencil_SeqSG,
/*140*/MatSetGrid_SeqSG
};


/** MatSetGrid_SeqSG : Sets the 3d physical grid information
Added by Deepan */
#undef __FUNCT__ 
#define __FUNCT__ "MatSetGrid_SeqSG"

PetscErrorCode MatSetGrid_SeqSG(Mat B, PetscInt m, PetscInt n, PetscInt p)
{
	Mat_SeqSG * b = (Mat_SeqSG*) B->data;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	if(m <= 0 || n <= 0 || p <= 0 ) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Grid Dimension should be atleast 1");

	b->m = m;
	b->n = n;
	b->p = p;
	b->nz = m*n*p*b->dof;
	PetscFunctionReturn(0);
}

/** MatCreate_SeqSG : Creates the struct grid representation 
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatCreate_SeqSG"

PetscErrorCode MatCreate_SeqSG(Mat B)
{
	Mat_SeqSG * b;
	PetscErrorCode ierr;
	PetscMPIInt size;
	PetscFunctionBegin;

	ierr = MPI_Comm_size(((PetscObject)B)->comm, &size); CHKERRQ(ierr);
	if (size > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Comm must be size 1");

	ierr = PetscMalloc(sizeof(Mat_SeqSG),&b);CHKERRQ(ierr);CHKERRQ(ierr);
	B->data = (void *)b;
	memcpy(B->ops,&MatOps_Values,sizeof(struct _MatOps));
	B->same_nonzero= PETSC_FALSE;
	B->spptr = 0;
	b->a = 0;
	b->m = 0;
	b->n = 0;
	b->p = 0;
	b->dof = 0;
	b->nz = 0;
	b->xt = 0;
	b->idx = PETSC_NULL;
	b->idy = PETSC_NULL;
	b->idz = PETSC_NULL;
	ierr = PetscObjectChangeTypeName((PetscObject)B, MATSTRUCTGRID); CHKERRQ(ierr);
	PetscFunctionReturn(0);	
}

/** MatDestroy_SeqSG : Destroys the struct grid representation
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatDestroy_SeqSG"
PetscErrorCode MatDestroy_SeqSG(Mat A)
{
	Mat_SeqSG * a = (Mat_SeqSG *)A->data;
	PetscErrorCode ierr;

	PetscFunctionBegin;
#if defined(PETSC_USE_LOG)
	PetscLogObjectState((PetscObject)A,"X-size= %D, Y-size=%D, Z-size=%D, NZ=%D",a->m,a->n,a->p,a->nz*a->stpoints);
#endif
	ierr = PetscFree(a->a);CHKERRQ(ierr);
	ierr = PetscFree(a->xt);CHKERRQ(ierr);
	ierr = PetscFree(a->idx);CHKERRQ(ierr);
	ierr = PetscFree(a->idy);CHKERRQ(ierr);
	ierr = PetscFree(a->idz);CHKERRQ(ierr);
	ierr = PetscFree(a);CHKERRQ(ierr);

	ierr = PetscObjectChangeTypeName((PetscObject)A, 0);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/** MatMult_SeqSG : Performs the Matrix - Vector multiplication y= mat*x on the struct grid representation
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqSG"
PetscErrorCode MatMult_SeqSG(Mat mat, Vec x, Vec y)
{
	PetscErrorCode ierr;
	PetscInt size, i;
	Mat_SeqSG * a = (Mat_SeqSG *) mat->data;
	PetscScalar * v = a->a, *xx,*yy;
	
	PetscFunctionBegin;
	ierr = VecSet(y,0.0); CHKERRQ(ierr);
	ierr = VecGetArray(x, &xx); CHKERRQ(ierr);
	ierr = VecGetArray(y, &yy); CHKERRQ(ierr);
	//openmp version
	//ierr = SG_MatMultOpenmp(v,xx,yy,a->xt,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->dof,a->stpoints); CHKERRQ(ierr);
//......matmultintrinsic version 2
	size =  (a->stpoints/(2*a->dof-1))*a->m*a->n*a->p - (2 * (1+a->m));
	if(a->p > 1)
	{
		ierr = SG_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->dof,a->stpoints,3); 
		size -= 2*a->m*a->n;
	}
	else
		ierr = SG_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->dof,a->stpoints,2); 
	CHKERRQ(ierr);


//......matmultintrinsic
//	ierr = SG_MatMult(v,xx,yy,a->xt,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->dof,a->stpoints); CHKERRQ(ierr);

// version 2 below
//        ierr = SG_MatMult(v,a->xt,yy,xx,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->dof,a->stpoints); CHKERRQ(ierr);
       	
// version 1 below
//        ierr = SG_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->dof,a->stpoints); CHKERRQ(ierr);
	ierr = VecRestoreArray(x,&xx); CHKERRQ(ierr);
	ierr = VecRestoreArray(y,&yy); CHKERRQ(ierr);
	ierr = PetscLogFlops(2*size*a->dof*a->dof); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}




/** SG_MatMult : (version 1) Performs the matrix vector multiplication on linearized matrix representation coeff and the vector x
 Added by Deepan */
/*#undef __FUNCT__
#define __FUNCT__ "SG_MatMult"
PetscErrorCode SG_MatMult(PetscScalar * coeff, PetscScalar * x, PetscScalar * y,PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos )
{
	PetscInt i,j,k,l,xdisp,ydisp,zdisp;
	PetscInt lda1 = m*n*p*dof;
	PetscInt lda2 = m*n*dof;
	PetscInt lda3 = m*dof;

	PetscInt vbeg, vend,xval, offset;
	for(l=0;l<nos;l++)
	{
		xdisp = idx[l]; ydisp = idy[l] ; zdisp = idz[l]; offset = l*lda1;
                //printf("offset: %d\n",offset);
	 	xval = xdisp + ydisp*lda3 + zdisp*lda2;
// vbeg and vend specifies the starting and ending indices of the y vector that gets updated.
//Added by Deepan 
			vbeg = 0;
			vend = m*dof*n*p;
		if(xval > 0) 
		{	
			vend -= xval;
		} 	
		else
		{
			vbeg -= xval;
			xval = 0;
		}
		for(k=vbeg;(k+4)<vend;k+=4,xval+=4)
		{
			
		        y[k]   += (coeff[offset+k] * x[(xval)]);
			y[k+1] += (coeff[offset+k+1] * x[(xval+1)]);
			y[k+2] += (coeff[offset+k+2] * x[(xval+2)]);
			y[k+3] += (coeff[offset+k+3] * x[(xval+3)]);
		}	
		for(;k<vend;k++,xval++)
			y[k] += (coeff[offset+k] * x[(xval)]);
	}
	PetscFunctionReturn(0);
}
*/





/*
#undef __FUNCT__
#define __FUNCT__ "SG_MatMult"
PetscErrorCode SG_MatMult(PetscScalar * coeff, PetscScalar * xi, PetscScalar * x,PetscScalar * y, PetscInt * idx, PetscInt * idy, PetscInt * idz, PetscInt m, PetscInt n, PetscInt p,PetscInt dof, PetscInt nos )
{
	PetscInt i,j,k,l,xdisp,ydisp,zdisp;
	PetscInt lda1 = m*n*p*dof;
	PetscInt lda2 = m*n*dof;
	PetscInt lda3 = m*dof;
	__m256d yv, xv, coeffv,xv1,coeffv1, xv2, coeffv2, xv3, coeffv3;
		
	PetscInt xval[nos], offset[nos];
	memcpy(x+lda2,xi,sizeof(PetscScalar)*lda1);
	for(l=0;l<nos;l++)
	{
		xdisp = idx[l]; ydisp = idy[l] ; zdisp = idz[l]; offset[l] = l*lda1;
	 	xval[l] = xdisp + ydisp*lda3 + zdisp*lda2;
	}
	for(k=0;(k+4)<lda1;k+=4)
	{
		yv = _mm256_loadu_pd((PetscScalar *)(y+k));
		for(l=0;(l+4)<nos;l+=4)
		{
			xv = _mm256_loadu_pd((PetscScalar *)(x+lda2+xval[l]+k));
			xv1 = _mm256_loadu_pd((PetscScalar *)(x+lda2+xval[l+1]+k));
			xv2 = _mm256_loadu_pd((PetscScalar *)(x+lda2+xval[l+2]+k));
			xv3 = _mm256_loadu_pd((PetscScalar *)(x+lda2+xval[l+3]+k));
			coeffv = _mm256_loadu_pd((PetscScalar *)(coeff+offset[l]+k));
			coeffv1 = _mm256_loadu_pd((PetscScalar *)(coeff+offset[l+1]+k));
			coeffv2 = _mm256_loadu_pd((PetscScalar *)(coeff+offset[l+2]+k));
			coeffv3 = _mm256_loadu_pd((PetscScalar *)(coeff+offset[l+3]+k));
			yv = _mm256_add_pd(yv,_mm256_mul_pd(coeffv,xv));
			yv = _mm256_add_pd(yv,_mm256_mul_pd(coeffv1,xv1));
			yv = _mm256_add_pd(yv,_mm256_mul_pd(coeffv2,xv2));
			yv = _mm256_add_pd(yv,_mm256_mul_pd(coeffv3,xv3));
			//y[k] += (coeff[offset+k] * x[(xval)]);
			//y[k+1] += (coeff[offset+k+1] * x[(xval+1)]);
			//y[k+2] += (coeff[offset+k+2] * x[(xval+2)]);
			//y[k+3] += (coeff[offset+k+3] * x[(xval+3)]);
			

}
		for(;l<nos;l++)
		{
			xv = _mm256_loadu_pd((PetscScalar *)(x+lda2+xval[l]+k));
			coeffv = _mm256_loadu_pd((PetscScalar *)(coeff+offset[l]+k));
			yv = _mm256_add_pd(yv,_mm256_mul_pd(coeffv,xv));
		}
		_mm256_storeu_pd((PetscScalar *)(y+k),yv);	
	}
	for(;k<lda1;k++)
		for(l=0;l<nos;l++)
	{
			y[k] += (coeff[offset[l]+k] * x[lda2+(xval[l]+k)]);
		
	}
	PetscFunctionReturn(0);
}

*/


/** MatSetValuesBlocked_SeqSG : Sets the values in the matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_SeqSG"

PetscErrorCode MatSetValuesBlocked_SeqSG(Mat A, PetscInt nrow,const PetscInt irow[], PetscInt ncol,const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
	PetscErrorCode ierr;
	Mat_SeqSG *mat = (Mat_SeqSG*)A->data;
	PetscInt * idm, *idn, i, j, bs = mat->dof; 
	PetscFunctionBegin;
	ierr = PetscMalloc(sizeof(PetscInt)* nrow * bs,&idm);CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)* ncol * bs,&idn);CHKERRQ(ierr);
	for(i=0;i<nrow;i++)
		for(j=0;j<bs;j++)
			idm[i*bs+j]= irow[i]*bs+j;
	for(i=0;i<ncol;i++)
		for(j=0;j<bs;j++)
			idn[i*bs+j]= icol[i]*bs+j;
	MatSetValues_SeqSG(A,nrow*bs,idm,ncol*bs,idn,y,is);
	ierr = PetscFree(idm);CHKERRQ(ierr);
	ierr = PetscFree(idn);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}	

/** MatSetValuesBlocked_SeqSG : Sets the values in the matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetValues_SeqSG"

PetscErrorCode MatSetValues_SeqSG(Mat A, PetscInt nrow,const PetscInt irow[], PetscInt ncol ,const PetscInt icol[],const PetscScalar y[], InsertMode is)
{//(*AIJ,1,&i,nz,cols,vals,INSERT_VALUES);
  // printf("..............call to MatSetValues_SeqSG()\n");
	PetscErrorCode ierr;
	Mat_SeqSG * mat = (Mat_SeqSG *) A->data;
	PetscInt * idx, * idy, * idz;
	PetscInt i,j,count = 0, m,n,p, offset,dis,xdis,ydis,zdis, cdis, k ,stp,dof,rshift,cshift;
	PetscFunctionBegin;	
	
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&idx);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&idy);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&idz);CHKERRQ(ierr);
	
	m = mat->m;
	n = mat->n;
	p = mat->p;
	stp = mat->stpoints;
	dis = mat->dis;
	dof = mat->dof;
	
        //printf("In MatSetValues_SeqSG\n");
        //printf("m=%d, n=%d, p=%d, stp=%d, dis=%d, dof=%d\n",m,n,p,stp,dis,dof);
//	ierr = PetscIntView(m,PETSC_VIWER_STDOUT_WORLD);CHKERRQ(ierr);
//        ierr = PetscIntView(n,PETSC_VIWER_STDOUT_WORLD);CHKERRQ(ierr);
//        ierr = PetscIntView(p,PETSC_VIWER_STDOUT_WORLD);CHKERRQ(ierr);
//        ierr = PetscIntView(stp,PETSC_VIWER_STDOUT_WORLD);CHKERRQ(ierr);
//        ierr = PetscIntView(dis,PETSC_VIWER_STDOUT_WORLD);CHKERRQ(ierr);
//        ierr = PetscIntView(dof,PETSC_VIWER_STDOUT_WORLD);CHKERRQ(ierr);
	fflush(stdout);
	for(i=0;i< nrow ; i++)
	{
		rshift = irow[i]%dof;
		for(j=0;j<ncol;j++)
		{
			cshift = icol[j]%dof;	
			cdis = (icol[j] - irow[i]) + rshift - cshift;
			zdis = (cdis)/(m*n*dof);
			cdis = (cdis)%(m*n*dof);
			ydis = (cdis)/(m*dof);
			cdis = (cdis)%(m*dof);
	//		xdis = (cdis%dof);
			xdis = cdis + cshift - rshift;
			for(k=0;k<stp;k++)
				if(mat->idx[k] == xdis &&  mat->idy[k] == ydis &&  mat->idz[k] == zdis)	
				{
					offset = k; //corresponding band
					break;
				}
			idx[count] = irow[i] % (m*dof);
			idy[count] = (irow[i]/(m*dof))%n;
			//printf("irow[%d]: %d, m: %d, n: %d, dof: %d\n",i,irow[i],m,n,dof);
			idz[count] = (offset*p) + ((irow[i]/(m*dof*n)) %p) ;
			count ++;
		}	
	}
	ierr = SetValues_SeqSG(mat,nrow*ncol,idx,idy,idz,y,is); CHKERRQ(ierr);
	ierr = PetscFree(idx);CHKERRQ(ierr);
	ierr = PetscFree(idy);CHKERRQ(ierr);
	ierr = PetscFree(idz);CHKERRQ(ierr);

	
	//Set the flag that indicates that matrix has changed on the CPU side.
        //This flag is used while copying the matrix to GPU.
        //Added by Chekuri S. Choudary 
	if (A->valid_GPU_matrix != PETSC_CUSP_UNALLOCATED)
	    A->valid_GPU_matrix = PETSC_CUSP_CPU;

	PetscFunctionReturn(0);
}

/** MatSetValuesBlocked_SeqSG : Sets the values in the matrix with the 3d indices supplied
Added by Deepan */
PetscErrorCode SetValues_SeqSG(Mat_SeqSG *  mat, PetscInt n , const PetscInt idx[], const PetscInt idy[],const PetscInt idz[],const  PetscScalar data[], InsertMode is)
{
	PetscInt i, mx = mat->m, ny= mat->n, pz = mat->p, dof = mat->dof;
	PetscInt lda1 = mx*dof*ny, lda2 = mx*dof;
	if(is == ADD_VALUES)
	{
		for(i=0;i<n;i++)
			mat->a[lda1*idz[i]+lda2*idy[i]+idx[i]] += data[i];
	}
	else
	{
		for(i=0;i<n;i++)
			mat->a[lda1*idz[i]+lda2*idy[i]+idx[i]] = data[i];
	}


	PetscFunctionReturn(0);
}


/** MatSetStencil_SeqSG : Sets the stencil neighbor points
 * sets various stencil related info such as neighbor points displacements and dimension of the grid
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetStencil_SeqSG"

PetscErrorCode MatSetStencil_SeqSG(Mat A, PetscInt dim,const PetscInt dims[],const PetscInt starts[], PetscInt dof)
{
	Mat_SeqSG * mat = (Mat_SeqSG *)A->data;
	PetscInt i,cnt=0;
	PetscErrorCode ierr;

	mat->dof = dof;
	// number of stencil displacements per stencil neighbor is given by 2*dof-1
	mat->stpoints = (2*dim+1)*(2*dof-1);
	// neighbors are considered to be one step away
	mat->dis = 1;
	mat->m=1;
	mat->n=1;
	mat->p=1;
	if(dim>0)
	mat->m=dims[0];
	if(dim>1)
	mat->n=dims[1];
	if(dim>2)
	mat->p=dims[2];
	mat->nz = mat->dof * mat->m * mat->n * mat->p;
	//printf("m=%d, n=%d,p=%d\n",mat->m,mat->n,mat->p);

	ierr = PetscMalloc(sizeof(PetscInt)*mat->stpoints,&(mat->idx));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*mat->stpoints,&(mat->idy));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*mat->stpoints,&(mat->idz));CHKERRQ(ierr);
/*Do not change the order of the stencil. MatMult depends on the order of these stencils.*/	
	mat->idx[cnt] = 0; mat->idy[cnt]=0; mat->idz[cnt++]= 0;
	for(i=1;i<dof;i++)
	{
		mat->idx[cnt] = i; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
	}
	if(dim>0)
	{	
		mat->idx[cnt] = dof; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = dof+i; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
		}
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = dof-i; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
		}
	}	
	if(dim>1)
	{
		mat->idx[cnt] = 0; mat->idy[cnt] = 1; mat->idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = i; mat->idy[cnt] = 1; mat->idz[cnt++] = 0;
		}
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = -i; mat->idy[cnt] = 1; mat->idz[cnt++] = 0;
		}
	}	
	if(dim>2)
	{
		mat->idx[cnt] = 0; mat->idy[cnt] = 0; mat->idz[cnt++] = 1;
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = i; mat->idy[cnt] = 0; mat->idz[cnt++] = 1;
		}
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = -i; mat->idy[cnt] = 0; mat->idz[cnt++] = 1;
		}
	}	

	for(i=1;i<dof;i++)
	{
		mat->idx[cnt] = -i; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
	}
	if(dim>0)
	{	
		mat->idx[cnt] = -dof; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = -dof+i; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
		}
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = -dof-i; mat->idy[cnt] = 0; mat->idz[cnt++] = 0;
		}
	}
	if(dim>1)
	{
		mat->idx[cnt] = 0; mat->idy[cnt] = -1; mat->idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = i; mat->idy[cnt] = -1; mat->idz[cnt++] = 0;
		}
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = -i; mat->idy[cnt] = -1; mat->idz[cnt++] = 0;
		}
	}
	if(dim>2)
	{
		mat->idx[cnt] = 0; mat->idy[cnt] = 0; mat->idz[cnt++] = -1;
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = i; mat->idy[cnt] = 0; mat->idz[cnt++] = -1;
		}
		for(i=1;i<dof;i++)
		{
			mat->idx[cnt] = -i; mat->idy[cnt] = 0; mat->idz[cnt++] = -1;
		}
	}
 	PetscFunctionReturn(0);	


}

/** MatSetUpPreallocation_SeqSG : Allocates space for coefficient matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetUpPreallocation_SeqSG"

PetscErrorCode MatSetUpPreallocation_SeqSG(Mat mat)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	Mat_SeqSG * a = (Mat_SeqSG *)mat->data;
	PetscFunctionBegin;
	ierr = PetscMalloc(sizeof(PetscScalar)*a->nz*a->stpoints,&(a->a));CHKERRQ(ierr);
	memset(a->a, 0,sizeof(PetscScalar)*a->nz*a->stpoints);
	mat->preallocated = PETSC_TRUE;
	PetscFunctionReturn(0);
}

/** MatZeroEntries_SeqSG : Sets the values in the matrix to be zero
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_SeqSG"

PetscErrorCode MatZeroEntries_SeqSG(Mat A)
{
	Mat_SeqSG * a = (Mat_SeqSG *)A->data;
	PetscFunctionBegin;
	memset(a->a,0,sizeof(PetscScalar)*a->nz*a->stpoints);
	if (A->valid_GPU_matrix != PETSC_CUSP_UNALLOCATED)
	    A->valid_GPU_matrix = PETSC_CUSP_CPU;

	PetscFunctionReturn(0);
}

/** MatGetDiagonal_SeqSG : Returns the diagonal elements of the matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_SeqSG"

PetscErrorCode MatGetDiagonal_SeqSG(Mat A, Vec v)
{
	PetscErrorCode ierr;
	Mat_SeqSG * a = (Mat_SeqSG *) A->data;
	PetscScalar *x;
	PetscInt i,j,n,dof = a->dof,offset;
	PetscFunctionBegin;
	ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
  	if (n != (a->nz)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Nonconforming matrix and vector");
	ierr = VecGetArray(v,&x);CHKERRQ(ierr);
/** Coefficients corresponding to the stencil 0,0,0(a[0] to a[n]) is returned
Added by Deepan */
	for(i=0;i<n;i++){
          x[i] = a->a[i];
          //printf("x[%d]: %f\n",i,x[i]);
        }
	for(i=0;i<n;i++)
	 	x[i] = a->a[i];
	ierr = VecRestoreArray(v,&x);CHKERRQ(ierr);
 	PetscFunctionReturn(0);
}


/** MatDiagonalScale
    Added by Justin Holewinski */
#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_SeqSG"

PetscErrorCode MatDiagonalScale_SeqSG(Mat A, Vec ll, Vec rr) {
  PetscFunctionBegin;
  fprintf(stderr, "WARNING: MatDiagonalScale_SeqSG is a placeholder!\n");
  PetscFunctionReturn(0);
}


/** MatGetRow_SeqSG : Returns the element corresponding to a single row in the matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatGetRow_SeqSG"

PetscErrorCode MatGetRow_SeqSG(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v)
{

  //printf("..........MatGetRow_SeqSG() called\n");
	Mat_SeqSG * a = (Mat_SeqSG *) A->data;
	PetscInt j;
        // printf("grA\n");
        PetscErrorCode ierr;
	PetscFunctionBegin;

	if (row < 0 || row >= A->rmap->n) 
          SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Row %D out of range",row);
        // printf("grB\n");
	*nz = a->stpoints*a->dof;
        // printf("grC\n");

        if(idx!=PETSC_NULL){
          //  printf("in idx PetscMalloc start.\n");
          ierr = PetscMalloc(sizeof(PetscInt)*a->stpoints*a->dof,idx);CHKERRQ(ierr);
          for(j=0;j<a->stpoints;j++)(*idx)[j] = j;

	}
        // printf("grD\n");
        if(v!=PETSC_NULL){
          //  printf("in v PetscMalloc start.\n");
           ierr = PetscMalloc(sizeof(PetscScalar)*a->stpoints*a->dof,v);CHKERRQ(ierr);
           for(j=0;j<a->stpoints;j++)(*v)[j] = a->a[row+ j*a->nz];
         }
        //  printf("grE\n");
	PetscFunctionReturn(0);
}

/** MatRestoreRow_SeqSG : Used in conjunction with MatGetRow to clear temporary data
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_SeqSG"

PetscErrorCode MatRestoreRow_SeqSG(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v)
{
  //printf("..........MatRestoreRow_SeqSG() called\n");
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(*idx);CHKERRQ(ierr);
  ierr = PetscFree(*v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** MatGetRowMaxAbs_SeqSG : Gets the absolute maximum of the elements in a particular row
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_SeqSG"

PetscErrorCode MatGetRowMaxAbs_SeqSG(Mat A, Vec v, PetscInt idx[])
{
	Mat_SeqSG * a = (Mat_SeqSG *) A->data;
	PetscScalar *x;
	PetscReal atmp;
	PetscInt i,j,n;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	if (A->factortype) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Not for factored matrix");
	ierr = VecSet(v,0.0);CHKERRQ(ierr);
 	ierr = VecGetArray(v,&x);CHKERRQ(ierr);
  	ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
	if (n != A->rmap->n) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Nonconforming matrix and vector");
	for(i=0;i<a->nz;i++)
	{
		x[i] = 0.0;
		for(j=0;j<a->stpoints;j++)
		{
                  atmp = PetscAbsScalar(a->a[i+ j*a->nz]);
			if(atmp > x[i]) x[i] = atmp;
		}
	}
	ierr = VecRestoreArray(v,&x);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "MatConvert_SGtoAIJ"
/*
  MatConvert_SGtoAIJ - Converts from any inputSeqstructgrid format to seqaij.
  Daniel Lowell*/
PetscErrorCode MatConvert_SGtoAIJ(Mat A, Mat *AIJ){
  //printf(".................MatConvert_SGtoAIJ() called\n");
  PetscFunctionBegin;
  Mat_SeqSG * a = (Mat_SeqSG *) A->data;
  PetscInt i,j;
  PetscInt m = A->rmap->n,n = A->cmap->n,nz = a->nz;
  PetscScalar *vals;
  PetscInt *cols;
  PetscErrorCode ierr;
  PetscInt *nnz;
  ierr = PetscMalloc(m*sizeof(PetscInt),&nnz); CHKERRQ(ierr);
  for(i=0;i<m;i++)nnz[i]=a->stpoints*a->dof;
  printf("m: %d, n: n: %d\n",m,n);
  ierr = MatCreateSeqAIJ(MPI_COMM_WORLD,m,n,PETSC_NULL,PETSC_NULL,AIJ);CHKERRQ(ierr);
  Mat_SeqAIJ     *aij = (Mat_SeqAIJ*)(*AIJ)->data;
  PetscInt       *ai = aij->i;

  for(i=0;i<m;i++){
      nnz[i]=0;
      ierr = MatGetRow_SeqSG(A,i,&nnz[i],&cols,&vals); CHKERRQ(ierr);
      ierr = MatSetValues(*AIJ,1,&i,nnz[i],cols,vals,INSERT_VALUES);
      ierr = MatRestoreRow_SeqSG(A,i,&nnz[i],&cols,&vals);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(*AIJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*AIJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=PetscFree(nnz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "MatView_SeqSG"
PetscErrorCode MatView_SeqSG(Mat A,PetscViewer viewer){
	//printf(".................MatView_SeqSG() called\n");
	Mat_SeqSG * a = (Mat_SeqSG *) A->data;
	PetscFunctionBegin;
        PetscInt i,j;
        PetscInt m = A->rmap->n,n = A->cmap->n,nz = a->nz;
        PetscScalar *vals;
        PetscInt *cols;
	PetscErrorCode ierr;
	const MatType mtype = MATSEQAIJ;
	Mat AIJ;

        ierr = MatConvert_SGtoAIJ(A,&AIJ); CHKERRQ(ierr);
        ierr = MatView(AIJ,viewer);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


