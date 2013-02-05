#define PETSCMAT_DLL

#include "../src/mat/impls/blockstructgrid/seq/matblockstructgrid.h"
#include "../src/mat/impls/blockstructgrid/seq/commonfunctions.h"
#include "petscblaslapack.h"
#include "petscbt.h"
#include "petscmat.h"
#include <string.h>
#include "../src/mat/impls/aij/seq/aij.h"  

static struct _MatOps MatOps_Values = {
/*0*/ MatSetValues_SeqBSG,0,0,MatMult_SeqBSG,0,
/*5*/0,0,0,0,0,
/*10*/0,0,0,0,0,
/*15*/0,0,0,MatDiagonalScale_SeqBSG,0,
/*20*/0,0,0,MatZeroEntries_SeqBSG,0,
/*25*/0,0,0,0,MatSetUpPreallocation_SeqBSG,
/*30*/0,0,0,0,0,
/*35*/0,0,0,0,0,
/*40*/0,0,0,0,0,
/*45*/0,0,0,0,0,
/*50*/0,0,0,0,0,
/*55*/0,0,0,MatSetValuesBlocked_SeqBSG,MatGetSubMatrix_SeqBSG,
/*60*/MatDestroy_SeqBSG,MatView_SeqBSG,0,0,0,
/*65*/0,0,MatSetValues_SeqBSG,0,0,
/*70*/0,0,0,0,0,
/*75*/0,0,0,0,0,
/*80*/0,0,0,0,0,
/*85*/0,0,MatSetValuesBlocked_SeqBSG,0,0,
/*90*/0,0,0,0,0,
/*95*/0,0,0,0,0,
/*100*/0,0,0,0,0,
/*105*/0,0,0,0,0,
/*110*/0,0,0,0,0,
/*115*/MatCreate_SeqBSG,0,0,0,0,
/*120*/0,0,0,0,0,
/*125*/0,0,0,0,0,
/*130*/0,0,0,0,0,
/*135*/0,0,0,0,0,
/*140*/0,MatSetStencil_SeqBSG,MatSetGrid_SeqBSG
};


/** MatSetGrid_SeqBSG : Sets the 3d physical grid information*/
#undef __FUNCT__ 
#define __FUNCT__ "MatSetGrid_SeqBSG"

PetscErrorCode MatSetGrid_SeqBSG(Mat B, PetscInt m, PetscInt n, PetscInt p)
{
	Mat_SeqBSG * b = (Mat_SeqBSG*) B->data;
	PetscFunctionBegin;

	if(m <= 0 || n <= 0 || p <= 0 ) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Grid Dimension should be atleast 1");

	b->m = m;
	b->n = n;
	b->p = p;
	b->nz = m*n*p;
	PetscFunctionReturn(0);
}

/** MatCreate_SeqBSG : Creates the struct grid representation */ 
#undef __FUNCT__
#define __FUNCT__ "MatCreate_SeqBSG"

PetscErrorCode MatCreate_SeqBSG(Mat B)
{
	Mat_SeqBSG * b;
	PetscErrorCode ierr;
	PetscMPIInt size;
	PetscFunctionBegin;

	ierr = MPI_Comm_size(((PetscObject)B)->comm, &size); CHKERRQ(ierr);
	if (size > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Comm must be size 1");

	ierr = PetscMalloc(sizeof(Mat_SeqBSG),&b);CHKERRQ(ierr);CHKERRQ(ierr);
	B->data = (void *)b;
	memcpy(B->ops,&MatOps_Values,sizeof(struct _MatOps));
	B->same_nonzero= PETSC_FALSE;
	B->spptr = 0;
	b->a = PETSC_NULL;
	b->diag = PETSC_NULL;
	b->coeff = 0;
	b->m = 0;
	b->n = 0;
	b->p = 0;
	b->dof = 0;
	b->nz = 0;
	b->ioff = PETSC_NULL;
	b->stpoffset = PETSC_NULL;
	b->stencil_stride = 1;
	b->stencil_rbeg = 0;
	b->stencil_rend = 0;
	b->stencil_cbeg = 0;
	b->stencil_cend = 0;
	b->sub_matrix = PETSC_FALSE;
	b->rstart = PETSC_NULL;
	b->lbeg = PETSC_NULL;
	b->lend = PETSC_NULL;

	b->block_arrangement = PETSC_NULL;

	b->bs = 0;
	ierr = PetscObjectChangeTypeName((PetscObject)B, MATBLOCKSTRUCTGRID); CHKERRQ(ierr);
	PetscFunctionReturn(0);	
}

/** MatDestroy_SeqBSG : Destroys the struct grid representation*/
#undef __FUNCT__
#define __FUNCT__ "MatDestroy_SeqBSG"
PetscErrorCode MatDestroy_SeqBSG(Mat A)
{
	Mat_SeqBSG * a = (Mat_SeqBSG *)A->data;
	PetscErrorCode ierr;

	PetscFunctionBegin;
#if defined(PETSC_USE_LOG)
	PetscLogObjectState((PetscObject)A,"X-size= %D, Y-size=%D, Z-size=%D, NZ=%D",a->m,a->n,a->p,a->nz*a->stpoints);
#endif
	ierr = PetscFree(a->stpoffset);CHKERRQ(ierr);
	ierr = PetscFree(a->coeff);CHKERRQ(ierr);
	ierr = PetscFree(a->a);CHKERRQ(ierr);
	ierr = PetscFree(a->diag);CHKERRQ(ierr);
	ierr = PetscFree(a->ioff);CHKERRQ(ierr);
#ifdef _VEC4
	ierr = PetscFree(a->block_arrangement);CHKERRQ(ierr);
#endif

	if(a->sub_matrix)
	{
		ierr = PetscFree(a->rstart);CHKERRQ(ierr);
		ierr = PetscFree(a->lbeg);CHKERRQ(ierr);
		ierr = PetscFree(a->lend);CHKERRQ(ierr);
	}
	ierr = PetscFree(a);CHKERRQ(ierr);

	ierr = PetscObjectChangeTypeName((PetscObject)A, 0);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/** MatMult_SeqBSG : Performs the Matrix - Vector multiplication y= mat*x on the struct grid representation*/
#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqBSG"
PetscErrorCode MatMult_SeqBSG(Mat mat, Vec x, Vec y)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * a = (Mat_SeqBSG *) mat->data;
	PetscScalar ** v = a->coeff,*yy;
	const PetscScalar *xx;
	
	PetscFunctionBegin;
 
#if !(defined _VEC2 || defined _VEC4)
	

#endif
	ierr = VecSet(y,0.0); CHKERRQ(ierr);
	ierr = VecGetArrayRead(x, &xx); CHKERRQ(ierr);
	ierr = VecGetArray(y, &yy); CHKERRQ(ierr);
	if(a->sub_matrix)
	{
		ierr = a->multfunc(v, xx, yy, a->ioff, a->m, a->n, a->p, a->dof, a->stpoints, 3, a->bs, &a->stpoffset[0], a->nregion, a->lbeg, a->lend, a->rstart);
	}else{
		ierr = a->multfunc(v, xx, yy, a->ioff, a->m, a->n, a->p, a->dof, a->stpoints, 3, a->bs, &a->stpoffset[0], a->nregion, a->lbeg, a->lend, a->rstart);
	}
	CHKERRQ(ierr);

	ierr = VecRestoreArrayRead(x,&xx); CHKERRQ(ierr);
	ierr = VecRestoreArray(y,&yy); CHKERRQ(ierr);
	ierr = PetscLogFlops(2*a->bs*a->nz*a->stpoints); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix*/
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_SeqBSG"

PetscErrorCode MatSetValuesBlocked_SeqBSG(Mat A, PetscInt nrow,const PetscInt irow[], PetscInt ncol,const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;	

	ierr = 	MatGetSetValues_SeqBSG(A,  nrow, irow, ncol, icol, &y[0], is, PETSC_TRUE, PETSC_TRUE); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}	

/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix*/
#undef __FUNCT__
#define __FUNCT__ "MatSetValues_SeqBSG"

PetscErrorCode MatSetValues_SeqBSG(Mat A, PetscInt nrow,const PetscInt irow[], PetscInt ncol ,const PetscInt icol[],const PetscScalar y[], InsertMode is)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;	

	ierr = 	MatGetSetValues_SeqBSG(A,  nrow, irow, ncol, icol, &y[0], is,PETSC_TRUE, PETSC_FALSE); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode MatGetSetValues_SeqBSG(Mat A, PetscInt nrow,const PetscInt irow[], PetscInt ncol ,const PetscInt icol[],const PetscScalar *y, InsertMode is, PetscBool SetVal, PetscBool Blocked)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * mat = (Mat_SeqBSG *) A->data;
	PetscInt *ipos, *ioffsets, *icols;
	PetscInt i,j,count = 0, m,n,p, offset,dis,xdis,ydis,zdis, cdis, k ,stp,dof,rshift,cshift, rowval, rbeg, cbeg, nstr, divRatio;

	
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&ipos);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&ioffsets);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&icols);CHKERRQ(ierr);
	
	m = mat->m;
	n = mat->n;
	p = mat->p;
	stp = mat->stpoints;
	dis = mat->dis;
	dof = mat->dof;
	rbeg = mat->stencil_rbeg;
	cbeg = mat->stencil_cbeg;
	nstr = mat->stencil_stride;
	if(Blocked)
		divRatio = 1;
	else
		divRatio = dof;
	for(i=0;i< nrow ; i++)
	{
		rshift = (irow[i]*nstr)/divRatio + rbeg;
		for(j=0;j<ncol;j++)
		{
			cshift = (icol[j]*nstr)/divRatio + cbeg;	
			cdis = cshift - rshift;
			for(k=0;k<stp;k++)
				if(mat->ioff[k] == cdis)	
				{
					offset = k; //corresponding band
					break;
				}
			if(k<stp)
			{
				ipos[count] = irow[i];
				icols[count] = icol[j];
				ioffsets[count] = k;
				count ++;
			}
		}	
	}
	if(SetVal)
	{

#if  (defined _VEC2 || defined _VEC4)

			if(Blocked)
			{
				ierr = SetValuesBlocked_Matrix_SeqBSG(mat,nrow*ncol,ipos,ioffsets,y,is); CHKERRQ(ierr);
			}
			else
		        {
				ierr = SetValues_Matrix_SeqBSG(mat,nrow*ncol,ipos,icols,ioffsets,y,is); CHKERRQ(ierr);
			}
#endif

	}
	else
	{
          ierr = GetValues_Matrix_SeqBSG(mat,nrow*ncol,ipos,ioffsets,(PetscScalar*)y); CHKERRQ(ierr);
	}
	ierr = PetscFree(ioffsets);CHKERRQ(ierr);
	ierr = PetscFree(ipos);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode GetValues_Matrix_SeqBSG(Mat_SeqBSG *  mat, PetscInt n , const PetscInt ipos[], const PetscInt ioffsets[], PetscScalar *data)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,c,endpoint,k1,mx = mat->m, ny= mat->n, dof = mat->dof, bs = mat->bs;
	PetscInt pos, l3threshold = WORKINGSETSIZE/bs;
	PetscInt *start, *lbeg, *lend, count, finalpos, nregion, ioff;
	l3threshold = WORKINGSETSIZE/mat->bs;
	start = mat->rstart;
	lbeg = mat->lbeg; 
	lend = mat->lend; 
	nregion = mat->nregion;	
	for(i=0;i<n;i++)
	{
		pos = ipos[i];
		count = 0;
		for(k=0;k<nregion;k++)
		{
			if(start[k] <= pos && pos < start[k+1])
			{
				finalpos = pos - start[k] ;
				ioff = count + ioffsets[i];// for submatrix - lbeg[k];
				if(dof % 2 == 0)
				{
						for(j=0;j<dof;j++)
							for(k1=0;k1<dof;k1++)
								data[i*bs+j*dof+k1] =	mat->coeff[ioff][(finalpos)*bs+j*dof+k1] ;
				}
				else
				{
						for(j=0;j<dof;j++)
						{
							for(k1=0;k1<dof-1;k1++)
								data[i*bs+j*dof+k1] =	mat->coeff[ioff][(finalpos)*bs+j*(dof-1)+k1];

								data[i*bs+(j+1)*dof-1] =	mat->coeff[ioff][(finalpos)*bs+dof*(dof-1)+j];
						}
				}

			}
			count += mat->stpoints;
			//count += (lend[k] - lbeg[k]);
		}
	}
	PetscFunctionReturn(0);
}


/** MatSetStencil_SeqBSG : Sets the stencil neighbor points
 * sets various stencil related info such as neighbor points displacements and dimension of the grid */
#undef __FUNCT__
#define __FUNCT__ "MatSetStencil_SeqBSG"

PetscErrorCode MatSetStencil_SeqBSG(Mat A, PetscInt dim,const PetscInt dims[],const PetscInt starts[], PetscInt dof)
{
	Mat_SeqBSG * mat = (Mat_SeqBSG *)A->data;
	PetscInt mnos = dim;
	PetscErrorCode ierr;

	mat->dof = dof;
	mat->bs = dof*dof;
	// number of stencil displacements per stencil neighbor is given by 2*dof-1
	mat->stpoints = (2*dim+1);
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
	mat->nz =  mat->m * mat->n * mat->p;

	ierr = PetscMalloc(sizeof(PetscInt)*mat->stpoints,&(mat->ioff));CHKERRQ(ierr);
	
	/*Do not change the order of the stencil. MatMult depends on the order of these stencils.*/	
	if(dim>2)
	{
		mat->ioff[mnos-3] = -(mat->m*mat->n); 
	}	
	if(dim>1)
	{
		mat->ioff[mnos-2] = -(mat->m);
	}	
	if(dim>0)
	{	
		mat->ioff[mnos-1] = -(1);
	}	
	mat->ioff[mnos] = 0;
	if(dim > 0)
	{
		mat->ioff[mnos+1] = 1; 
	}
	if(dim > 1)
	{
		mat->ioff[mnos+2] = (mat->m);
	}
	if(dim > 2)
	{
		mat->ioff[mnos+3] = (mat->m*mat->n);
	}


#if  (defined _VEC2 || defined _VEC4)

       	if(mat->dof == 1)
		mat->multfunc = &BSG_MatMult_1;
	else if(mat->dof == 2)
		mat->multfunc = &BSG_MatMult_2;
	else if (mat->dof == 3)
		mat->multfunc = &BSG_MatMult_3;
	else if (mat->dof == 4)
		mat->multfunc = &BSG_MatMult_4;
	else if (mat->dof == 5)
		mat->multfunc = &BSG_MatMult_5;
	else if (mat->dof == 6)
		mat->multfunc = &BSG_MatMult_6;
#endif


#ifdef _VEC4
	else if (mat->dof == 7)
	        mat->multfunc = &BSG_MatMult_7;
#endif


#ifdef _VEC2
   
        else if (mat->dof == 8)
		mat->multfunc = &BSG_MatMult_8;
	else if (mat->dof == 10)
		mat->multfunc = &BSG_MatMult_10;
	else if (mat->dof == 12)
		mat->multfunc = &BSG_MatMult_12;
	else if (mat->dof == 14)
		mat->multfunc = &BSG_MatMult_14;
	else if (mat->dof == 16)
		mat->multfunc = &BSG_MatMult_16;
	else if (mat->dof == 18)
		mat->multfunc = &BSG_MatMult_18;
	else if (mat->dof == 20)
		mat->multfunc = &BSG_MatMult_20;
	else if (mat->dof == 22)
		mat->multfunc = &BSG_MatMult_22;
	else if (mat->dof == 24)
		mat->multfunc = &BSG_MatMult_24;
	else if (mat->dof == 26)
		mat->multfunc = &BSG_MatMult_26;
	else if (mat->dof == 28)
		mat->multfunc = &BSG_MatMult_28;
	else if (mat->dof == 30)
		mat->multfunc = &BSG_MatMult_30;
#endif

#if  (defined _VEC2 || defined _VEC4)

	else if (mat->dof%2 == 0)
		mat->multfunc = &BSG_MatMult_Neven;
	else
		mat->multfunc = &BSG_MatMult_Nodd;
#endif


 	PetscFunctionReturn(0);	
}

/** MatSetUpPreallocation_SeqBSG : Allocates space for coefficient matrix */
#undef __FUNCT__
#define __FUNCT__ "MatSetUpPreallocation_SeqBSG"
PetscErrorCode MatSetUpPreallocation_SeqBSG(Mat mat)
{
	Mat_SeqBSG * a = (Mat_SeqBSG *)mat->data;
	PetscErrorCode ierr;
	PetscFunctionBegin;
	if(a->sub_matrix)
	{
		ierr = MatSetUpPreallocation_SubMatrix_SeqBSG(mat); CHKERRQ(ierr);
	}
	else
	{
		ierr = MatSetUpPreallocation_Matrix_SeqBSG(mat); CHKERRQ(ierr);
	}
	mat->preallocated = PETSC_TRUE;
	PetscFunctionReturn(0);
}


PetscErrorCode MatSetUpPreallocation_SubMatrix_SeqBSG(Mat mat)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * a = (Mat_SeqBSG *)mat->data;
	PetscInt i,j,k,strtval, dof = a->dof;
	a->nz = (a->stencil_rend-a->stencil_rbeg)/a->stencil_stride;
	ierr = PetscMalloc(sizeof(PetscScalar)*a->nz*a->bs*a->stpoints,&(a->a));CHKERRQ(ierr);
	memset(a->a, 0,sizeof(PetscScalar)*a->nz*a->bs*a->stpoints);
	ierr = PetscMalloc(sizeof(PetscScalar)*a->nz*a->bs,&(a->diag));CHKERRQ(ierr);
	memset(a->diag, 0,sizeof(PetscScalar)*a->nz*a->bs);

#ifdef _VEC4
	ierr = PetscMalloc(sizeof(PetscInt)*a->bs,&(a->block_arrangement));CHKERRQ(ierr);
	memset(a->block_arrangement, 0,sizeof(PetscInt)*a->bs);
	strtval = 0;
	for(i = 0; (i+4) <= dof ; i+=4)
	{
		for(j=0; (j+4) <= dof; j+=4,strtval += 16)
		{
			for(k = 0; k < 2 ; k++ ){
				a->block_arrangement [i*dof+k+j] = strtval+k;
				a->block_arrangement [i*dof+k+j+2] = strtval+4+k;
				a->block_arrangement [(i+1)*dof+k+j] = strtval+8+k;
				a->block_arrangement [(i+1)*dof+k+j+2] = strtval+12+k;
			}
			for(k = 2; k < 4 ; k++ ){
				a->block_arrangement [(i+2)*dof+k+j] = strtval+k;
				a->block_arrangement [(i+2)*dof+k+j-2] = strtval+4+k;
				a->block_arrangement [(i+3)*dof+k+j] = strtval+8+k;
				a->block_arrangement [(i+3)*dof+k+j-2] = strtval+12+k;
			}
		}
		for(; (j+2) <= dof; j+=2,strtval += 8)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [(i+2)*dof+j] = strtval+2;
				a->block_arrangement [(i+2)*dof+j+1] = strtval+3;
				a->block_arrangement [(i+1)*dof+j] = strtval+4;
				a->block_arrangement [(i+1)*dof+j+1] = strtval+5;
				a->block_arrangement [(i+3)*dof+j] = strtval+6;
				a->block_arrangement [(i+3)*dof+j+1] = strtval+7;
			
		}
		for(; j < dof; j++, strtval += 4)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [(i+1)*dof+j] = strtval+1;
				a->block_arrangement [(i+2)*dof+j] = strtval+2;
				a->block_arrangement [(i+3)*dof+j] = strtval+3;
		}
	}
	for(; (i+2) <= dof ; i+=2)
	{
		for(j=0; (j+4) <= dof; j+=4,strtval += 8)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [i*dof+j+2] = strtval+2;
				a->block_arrangement [i*dof+j+3] = strtval+3;
				a->block_arrangement [(i+1)*dof+j] = strtval+4;
				a->block_arrangement [(i+1)*dof+j+1] = strtval+5;
				a->block_arrangement [(i+1)*dof+j+2] = strtval+6;
				a->block_arrangement [(i+1)*dof+j+3] = strtval+7;
		}
		for(; (j+2) <= dof; j+=2,strtval += 4)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [(i+1)*dof+j] = strtval+2;
				a->block_arrangement [(i+1)*dof+j+1] = strtval+3;
		}
		for(; j < dof; j++, strtval += 2)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [(i+1)*dof+j] = strtval+1;
		}
	}
	for(; i < dof ; i++)
	{
		for(j=0; (j+4) <= dof; j+=4,strtval += 4)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [i*dof+j+2] = strtval+2;
				a->block_arrangement [i*dof+j+3] = strtval+3;
		}
		for(; (j+2) <= dof; j+=2,strtval += 2)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
		}
		for(; j < dof; j++, strtval += 1)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
		}
	}

#endif

	ierr = MatSetUpRegion_SeqBSG(mat); CHKERRQ(ierr);
	ierr = PetscLayoutSetBlockSize(mat->rmap,a->dof);CHKERRQ(ierr);
	ierr = PetscLayoutSetBlockSize(mat->cmap,a->dof);CHKERRQ(ierr);
	ierr = PetscLayoutSetUp(mat->rmap);CHKERRQ(ierr);
	ierr = PetscLayoutSetUp(mat->cmap);CHKERRQ(ierr);
	mat->preallocated = PETSC_TRUE;
	PetscFunctionReturn(0);
}



PetscErrorCode MatSetUpSubRegion_SeqBSG(Mat A)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * a = (Mat_SeqBSG*) A->data;
	PetscInt nregion = a->nregion;
	PetscInt * lbeg = a->lbeg;
	PetscInt * lend = a->lend;
	PetscInt * start = a->rstart;
	PetscInt i,c,j, l3threshold, endpoint;
	PetscInt count = 0;
	PetscInt rcount = 0;	
	l3threshold = WORKINGSETSIZE/a->bs;
	
	for(i=0;i<nregion;i++)
	{
		if(a->sub_matrix){
			for(c = start[i]; c < start[i+1]; c+= l3threshold)
			{
				count += (lend[i] - lbeg[i]);
				rcount ++;
			}
		}
		else{
			for(c = start[i]; c < start[i+1]; c+= l3threshold)
			{
				count += a->stpoints;//(lend[i] - lbeg[i]);
				rcount ++;
			}
		}
	}
	ierr = PetscMalloc(sizeof(PetscScalar*)*(count+1),&(a->coeff));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(rcount+1),&(a->stpoffset));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(rcount+1),&(a->rstart));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(rcount),&(a->lbeg));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(rcount),&(a->lend));CHKERRQ(ierr);

	a->nregion = rcount;
	
	a->coeff[0] = a->a;
	count = 0;
	rcount = 0;
	for(i=0;i<nregion;i++)
	{
		if(a->sub_matrix){
			for(c = start[i]; c < start[i+1]; c+= l3threshold)
			{
				endpoint = (c+l3threshold) < start[i+1] ? (c+l3threshold) : start[i+1]; 
				for(j=1;j<=(lend[i] - lbeg[i]);j++) // counter starting from 1 is very important
				{
					a->coeff[count+j] = a->coeff[count]+(j*(endpoint-c)*a->bs);
				}
				a->stpoffset[rcount] = count;
				count += (lend[i] - lbeg[i]);
				a->rstart[rcount] = c;
				a->lbeg[rcount] = lbeg[i];
				a->lend[rcount] = lend[i];
				rcount++;
			}
		}
		else
		{
			for(c = start[i]; c < start[i+1]; c+= l3threshold)
			{
				endpoint = (c+l3threshold) < start[i+1] ? (c+l3threshold) : start[i+1]; 
				//for(j=1;j<=(lend[i] - lbeg[i]);j++)
				for(j=1;j<= a->stpoints;j++) 
				{
					a->coeff[count+j] = a->coeff[count]+(j*(endpoint-c)*a->bs);
				}
				a->stpoffset[rcount] = count;
				//count += (lend[i] - lbeg[i]);
				count += a->stpoints;//(lend[i] - lbeg[i]);
				a->rstart[rcount] = c;
				a->lbeg[rcount] = lbeg[i];
				a->lend[rcount] = lend[i];
				rcount++;
			}
		}
	}
	a->stpoffset[rcount] = count;
	a->rstart[a->nregion] = start[nregion];

		ierr = PetscFree(start);CHKERRQ(ierr);
		ierr = PetscFree(lbeg);CHKERRQ(ierr);
		ierr = PetscFree(lend);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode MatSetUpRegion_SeqBSG(Mat A)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * a = (Mat_SeqBSG*) A->data;
	PetscInt irowbeg = a->stencil_rbeg;
	PetscInt irowend = a->stencil_rend;
	PetscInt icolbeg = a->stencil_cbeg;
	PetscInt icolend = a->stencil_cend;
	PetscInt nstride = a->stencil_stride;
 	PetscInt c,i,j,next_val, maxstpoints = a->stpoints;
	PetscInt *lbeg, *lend, *lstart, region_count, *ioffsets;
	PetscInt ts, te , lb, le;
	const PetscInt lda3 = a->m;
	const PetscInt lda2 = lda3 * a->n;
	
	ierr = PetscMalloc2(maxstpoints, PetscInt, &lbeg, maxstpoints, PetscInt, &lend);CHKERRQ(ierr);
	ierr = PetscMalloc (sizeof(PetscInt) * (2*maxstpoints), &lstart);CHKERRQ(ierr);
	ierr = PetscMalloc (sizeof(PetscInt) * (maxstpoints), &ioffsets);CHKERRQ(ierr);
	
	for(i=0;i<maxstpoints;i++){
		ioffsets[i] =  (a->ioff[i]);
		ts = irowbeg + ioffsets[i];
		te = irowend + ioffsets[i];
		if((icolbeg <= ts && ts < icolend)||(icolbeg <= te && te < icolend)){
			lbeg[i] = irowbeg;
			lend[i] = irowend;
			if(ts < icolbeg)
			{
				 lbeg[i] += (icolbeg - ts);
			}
			if(te > icolend)
			{
				 lend[i] += (icolend - te);
			}
		}
		else{
			lbeg[i] = -1;
			lend[i] = -1;
		}
		
	}
	region_count = 0;
	i=maxstpoints-1; j=maxstpoints-1;

	while(i >= 0 && j >= 0)
	{
		if(lbeg[i] < lend[j])
		{
			next_val = lbeg[i--];
			next_val += ((next_val-irowbeg)%nstride);
		}
		else{
			next_val = lend[j--];
			if(next_val != irowend)
                       		next_val -= ((next_val-irowbeg)%nstride);
		}
		if((region_count !=0 && lstart[region_count-1] == next_val) || next_val < 0 )
			continue;
		lstart[region_count++] = next_val;
	}	

	while( i >= 0)
	{
		next_val = lbeg[i--];
		next_val += ((next_val-irowbeg)%nstride);
		if((region_count !=0 && lstart[region_count-1] == next_val) || next_val < 0 )
			continue;
		lstart[region_count++] = next_val;
	}
	while( j >= 0)
	{
		next_val = lend[j--];
		if(next_val != irowend)
                	next_val -= ((next_val-irowbeg)%nstride);
		if(region_count !=0 && lstart[region_count-1] == next_val)
			continue;
		lstart[region_count++] = next_val;
	}
 	a->nregion = region_count-1; //Very important

	ierr = PetscMalloc(sizeof(PetscInt)*(region_count+1),&(a->rstart));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(region_count),&(a->lbeg));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(region_count),&(a->lend));CHKERRQ(ierr);

	for(c=0; c<a->nregion; c++){
		lb = maxstpoints+1; le = -1;	
		for(i= 0; i< maxstpoints; i++)
		{
			ts = lstart[c]+ioffsets[i];
			if(icolbeg <= ts && ts < icolend){
				lb = i;
				break;
			}
		}
		for(; i< maxstpoints; i++){
			ts = lstart[c]+ioffsets[i];
			if(icolbeg <= ts && ts < icolend){
				le = i+1;
			}
		}
		a->lbeg[c] = lb;
		a->lend[c] = le;
		a->rstart[c] = (lstart[c]-irowbeg );
		if(irowbeg != 0) a->rstart[c] += 1;
		a->rstart[c] /= nstride;
	}
	a->rstart[a->nregion] = (lstart[a->nregion] - irowbeg );
		if(irowbeg != 0) a->rstart[a->nregion] += 1;
		a->rstart[a->nregion] /= nstride;
	
	ierr = MatSetUpSubRegion_SeqBSG(A); CHKERRQ(ierr);
	PetscFree(ioffsets);	
	PetscFree(lstart);
	PetscFree2(lbeg, lend);
	PetscFunctionReturn(0);  
}



PetscErrorCode MatSetUpPreallocation_Matrix_SeqBSG(Mat mat)
{
	PetscErrorCode ierr;
	PetscInt i,c,j,k,strtval, l3threshold, endpoint;
	PetscInt start[8], count = 0;
	Mat_SeqBSG * a = (Mat_SeqBSG *)mat->data;
	PetscInt dof = a->dof;
	ierr = PetscMalloc(sizeof(PetscScalar)*a->nz*a->bs*a->stpoints,&(a->a));CHKERRQ(ierr);
	memset(a->a, 0,sizeof(PetscScalar)*a->nz*a->bs*a->stpoints);
	ierr = PetscMalloc(sizeof(PetscScalar)*a->nz*a->bs,&(a->diag));CHKERRQ(ierr);
	memset(a->diag, 0,sizeof(PetscScalar)*a->nz*a->bs);


#ifdef _VEC4
	ierr = PetscMalloc(sizeof(PetscInt)*a->bs,&(a->block_arrangement));CHKERRQ(ierr);
	memset(a->block_arrangement, 0,sizeof(PetscInt)*a->bs);
	strtval = 0;
	for(i = 0; (i+4) <= dof ; i+=4)
	{
		for(j=0; (j+4) <= dof; j+=4,strtval += 16)
		{
			for(k = 0; k < 2 ; k++ ){
				a->block_arrangement [i*dof+k+j] = strtval+k;
				a->block_arrangement [i*dof+k+j+2] = strtval+4+k;
				a->block_arrangement [(i+1)*dof+k+j] = strtval+8+k;
				a->block_arrangement [(i+1)*dof+k+j+2] = strtval+12+k;
			}
			for(k = 2; k < 4 ; k++ ){
				a->block_arrangement [(i+2)*dof+k+j] = strtval+k;
				a->block_arrangement [(i+2)*dof+k+j-2] = strtval+4+k;
				a->block_arrangement [(i+3)*dof+k+j] = strtval+8+k;
				a->block_arrangement [(i+3)*dof+k+j-2] = strtval+12+k;
			}
		}
		for(; (j+2) <= dof; j+=2,strtval += 8)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [(i+2)*dof+j] = strtval+2;
				a->block_arrangement [(i+2)*dof+j+1] = strtval+3;
				a->block_arrangement [(i+1)*dof+j] = strtval+4;
				a->block_arrangement [(i+1)*dof+j+1] = strtval+5;
				a->block_arrangement [(i+3)*dof+j] = strtval+6;
				a->block_arrangement [(i+3)*dof+j+1] = strtval+7;
			
		}
		for(; j < dof; j++, strtval += 4)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [(i+1)*dof+j] = strtval+1;
				a->block_arrangement [(i+2)*dof+j] = strtval+2;
				a->block_arrangement [(i+3)*dof+j] = strtval+3;
		}
	}
	for(; (i+2) <= dof ; i+=2)
	{
		for(j=0; (j+4) <= dof; j+=4,strtval += 8)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [i*dof+j+2] = strtval+2;
				a->block_arrangement [i*dof+j+3] = strtval+3;
				a->block_arrangement [(i+1)*dof+j] = strtval+4;
				a->block_arrangement [(i+1)*dof+j+1] = strtval+5;
				a->block_arrangement [(i+1)*dof+j+2] = strtval+6;
				a->block_arrangement [(i+1)*dof+j+3] = strtval+7;
		}
		for(; (j+2) <= dof; j+=2,strtval += 4)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [(i+1)*dof+j] = strtval+2;
				a->block_arrangement [(i+1)*dof+j+1] = strtval+3;
		}
		for(; j < dof; j++, strtval += 2)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [(i+1)*dof+j] = strtval+1;
		}
	}
	for(; i < dof ; i++)
	{
		for(j=0; (j+4) <= dof; j+=4,strtval += 4)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
				a->block_arrangement [i*dof+j+2] = strtval+2;
				a->block_arrangement [i*dof+j+3] = strtval+3;
		}
		for(; (j+2) <= dof; j+=2,strtval += 2)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
				a->block_arrangement [i*dof+j+1] = strtval+1;
		}
		for(; j < dof; j++, strtval += 1)
		{
				a->block_arrangement [i*dof+j] = strtval+0;
		}
	}

#endif

	l3threshold = WORKINGSETSIZE/a->bs;
	a->stencil_stride = 1;
	a->stencil_rbeg = 0;
	a->stencil_cbeg = 0;
	a->stencil_rend = a->nz;
	a->stencil_cend = a->nz;
	ierr = MatSetUpRegion_SeqBSG(mat); CHKERRQ(ierr);

	ierr = PetscLayoutSetBlockSize(mat->rmap,a->dof);CHKERRQ(ierr);
	ierr = PetscLayoutSetBlockSize(mat->cmap,a->dof);CHKERRQ(ierr);
	ierr = PetscLayoutSetUp(mat->rmap);CHKERRQ(ierr);
	ierr = PetscLayoutSetUp(mat->cmap);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}



PetscErrorCode MatGetSubMatrix_SeqBSG_Private(Mat A,IS isrow,IS iscol,MatReuse scall,Mat *B)
{
	Mat_SeqBSG    *a = (Mat_SeqBSG*)A->data,*c;
	PetscErrorCode ierr;
	PetscInt      i,k,kstart,kend;
	PetscInt       row,xrow,*ioff;
	const PetscInt *irow,*icol;
	PetscInt       nrows,ncols,bs=a->bs;
	Mat            C;
	PetscBool      sorted;
	PetscInt	rrange, crange, nstride=0;
	PetscInt	*dims, *starts;
	PetscScalar	*val;


	ierr = ISSorted(isrow,&sorted);CHKERRQ(ierr);
	if (!sorted) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"IS is not sorted");
	ierr = ISSorted(iscol,&sorted);CHKERRQ(ierr);
	if (!sorted) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"IS is not sorted");

	ierr = ISGetIndices(isrow,&irow);CHKERRQ(ierr);
	ierr = ISGetIndices(iscol,&icol);CHKERRQ(ierr);
	ierr = ISGetLocalSize(isrow,&nrows);CHKERRQ(ierr);
	ierr = ISGetLocalSize(iscol,&ncols);CHKERRQ(ierr);

	rrange = irow[nrows-1] - irow[0];
	crange = icol[ncols-1] - icol[0];
	if(nrows > 1)	nstride = irow[1] - irow[0];


	/* Create and fill new matrix */
	if (scall == MAT_REUSE_MATRIX) {
	c = (Mat_SeqBSG *)((*B)->data);

	if ((*B)->rmap->bs!=bs || !c->sub_matrix || c->nz != (rrange) || c->nz != (crange)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Submatrix wrong size");
	c->stencil_stride = nstride;
	c->stencil_rbeg = irow[0];
	c->stencil_cbeg = icol[0];
	c->stencil_rend = irow[nrows-1];
	c->stencil_cend = icol[ncols-1];
	ierr = MatSetUpRegion_SeqBSG(*B); CHKERRQ(ierr);
	C = *B;
	} else { 
	 PetscInt dim = 3;
	ierr = PetscMalloc(dim*sizeof(PetscInt),&dims);CHKERRQ(ierr);
	ierr = PetscMalloc(dim*sizeof(PetscInt),&starts);CHKERRQ(ierr);
	dims[0] = a->m; dims[1] = a->n; dims[2] = a->p;
	ierr = MatCreate(((PetscObject)A)->comm,&C);CHKERRQ(ierr);
	ierr = MatSetSizes(C,nrows*bs,ncols*bs,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
	ierr = MatSetType(C,((PetscObject)A)->type_name);CHKERRQ(ierr);
	ierr = MatSetStencil(C,dim, dims, starts, a->dof);
	c = (Mat_SeqBSG *)((C)->data);
	c->stencil_stride = nstride;
	c->stencil_rbeg = irow[0];
	c->stencil_cbeg = icol[0];
	c->stencil_rend = irow[nrows-1];
	c->stencil_cend = icol[ncols-1];
	c->sub_matrix = PETSC_TRUE;
	ierr = MatSetUpPreallocation_SeqBSG(C); CHKERRQ(ierr);
	ierr = PetscFree(dims); CHKERRQ(ierr);
	ierr = PetscFree(starts); CHKERRQ(ierr);
	}

	/* copy*/
	c = (Mat_SeqBSG *)(C->data);
	ierr = PetscMalloc(c->bs*sizeof(PetscScalar),&val);CHKERRQ(ierr);
	for (i=0; i<nrows; i++) {
		row    = irow[i];
		for (k=0; k<c->stpoints; k++) {
			kstart = row + a->ioff[k];
                        kend = (kstart-c->stencil_cbeg)%c->stencil_stride;
			if(kend == 0 && kstart >=0 )
                      	{
//				*mat_j++ = tcol - 1;
//				ierr     = PetscMemcpy(mat_a,a->a+k*bs2,bs2*sizeof(MatScalar));CHKERRQ(ierr);
//				mat_a   += bs2;
//				(*mat_ilen)++;
//				GetOriginalMatrixBlock
				ierr = MatGetSetValues_SeqBSG(A, 1, &row, 1 , &kstart, val, INSERT_VALUES, PETSC_FALSE, PETSC_TRUE); CHKERRQ(ierr);
				xrow = (row-c->stencil_rbeg)/c->stencil_stride;
                                kstart = (kstart-c->stencil_cbeg)/c->stencil_stride;
				MatSetValues_SeqBSG(C, 1, &xrow, 1 , &kstart, val,INSERT_VALUES);	
			}
		}
	}
	ierr = PetscFree(val); CHKERRQ(ierr);
			    
			  /* Free work space */
	ierr = ISRestoreIndices(iscol,&icol);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
			  
	ierr = ISRestoreIndices(isrow,&irow);CHKERRQ(ierr);
	*B = C;
	PetscFunctionReturn(0);
}





/*** MatGetSubMatrix_SeqBSG : REturns sub matrix
 * Should be strides
 * Stride should be same value
 * Range of row and cols should be same
 * */

#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrix_SeqBSG"
PetscErrorCode MatGetSubMatrix_SeqBSG(Mat A,IS isrow,IS iscol,MatReuse scall,Mat *B)
{
	Mat_SeqBSG    *a = (Mat_SeqBSG*)A->data;
	IS             is1,is2;
	PetscErrorCode ierr;
	PetscInt       *vary,*iary,nrows,ncols,i,bs=A->rmap->bs,count;
	PetscInt 	stride_count=0, rstride=0, cstride=0, rrange=0 , crange= 0;
	const PetscInt *irow,*icol;

	PetscFunctionBegin;
	ierr = ISGetIndices(isrow,&irow);CHKERRQ(ierr);
	ierr = ISGetIndices(iscol,&icol);CHKERRQ(ierr);
	ierr = ISGetLocalSize(isrow,&nrows);CHKERRQ(ierr);
	ierr = ISGetLocalSize(iscol,&ncols);CHKERRQ(ierr);
	  
	  /* Verify if the indices corespond to each element in a block 
	   and form the IS with compressed IS */
	ierr = PetscMalloc2(a->nz,PetscInt,&vary,a->nz,PetscInt,&iary);CHKERRQ(ierr);
	ierr = PetscMemzero(vary,a->nz*sizeof(PetscInt));CHKERRQ(ierr);
	for (i=0; i<nrows; i++) vary[irow[i]/bs]++;
	count = 0;
	for (i=0; i<a->nz; i++) {
	  if (vary[i]!=0 && vary[i]!=bs) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Index set does not match blocks");
	  if (vary[i]==bs) iary[count++] = i;
	}
	
	if(count > 1){
		rstride = iary[1] - iary[0];
		for(i=1; i<count-1; i++)
		{
			stride_count = iary[i+1] - iary[i];
			if(stride_count != rstride)	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Not constant strides");
		}
		rrange = iary[count-1] - iary[0];
	}
	
	ierr = ISCreateGeneral(PETSC_COMM_SELF,count,iary,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
	
	ierr = PetscMemzero(vary,(a->nz)*sizeof(PetscInt));CHKERRQ(ierr);
	for (i=0; i<ncols; i++) vary[icol[i]/bs]++;
	count = 0;
	for (i=0; i<a->nz; i++) {
	  if (vary[i]!=0 && vary[i]!=bs) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Internal error in PETSc");
	  if (vary[i]==bs) iary[count++] = i;
	}
	
	if(count > 1){
		cstride = iary[1] - iary[0];
		for(i=1; i<count-1; i++)
		{
			stride_count = iary[i+1] - iary[i];
			if(stride_count != cstride)	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Not constant strides");
		}
		crange = iary[count-1] - iary[0];
	}
	
	if(rstride > 0 && cstride > 0 && rstride != cstride)	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Not constant strides");
	
	if(rrange != crange)	SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Not square matrix");

	ierr = ISCreateGeneral(PETSC_COMM_SELF,count,iary,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
	ierr = ISRestoreIndices(isrow,&irow);CHKERRQ(ierr);
	ierr = ISRestoreIndices(iscol,&icol);CHKERRQ(ierr);
	ierr = PetscFree2(vary,iary);CHKERRQ(ierr);

	ierr = MatGetSubMatrix_SeqBSG_Private(A,is1,is2,scall,B);CHKERRQ(ierr);
	ierr = ISDestroy(&is1);CHKERRQ(ierr);
	ierr = ISDestroy(&is2);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}



/*** MatGetSubMatrices_SeqBSG : REturns sub matrices
 * Should be strides
 * Stride should be same value
 * Range of row and cols should be same
 * */

#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrices_SeqBSG"
PetscErrorCode MatGetSubMatrices_SeqBSG(Mat A,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *B[])
{
	PetscErrorCode ierr;
	PetscInt       i;

	PetscFunctionBegin;
	if (scall == MAT_INITIAL_MATRIX) {
		ierr = PetscMalloc((n+1)*sizeof(Mat),B);CHKERRQ(ierr);
	}

	for (i=0; i<n; i++) {
		ierr = MatGetSubMatrix_SeqBSG(A,irow[i],icol[i],scall,&(*B)[i]);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}


/** MatZeroEntries_SeqBSG : Sets the values in the matrix to be zero */
#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_SeqBSG"

PetscErrorCode MatZeroEntries_SeqBSG(Mat A)
{
	Mat_SeqBSG * a = (Mat_SeqBSG *)A->data;
	PetscFunctionBegin;
	memset(a->a,0,sizeof(PetscScalar)*a->nz*a->bs*a->stpoints);
	memset(a->diag, 0,sizeof(PetscScalar)*a->nz*a->bs);
#ifdef PETSC_HAVE_CUSP
	if (A->valid_GPU_matrix != PETSC_CUSP_UNALLOCATED)
	    A->valid_GPU_matrix = PETSC_CUSP_CPU;
#endif
	PetscFunctionReturn(0);
}



/** MatDiagonalScale_SeqBSG : Scales matrix based on left and right vectors */
#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_SeqBSG"

PetscErrorCode MatDiagonalScale_SeqBSG (Mat A, Vec ll, Vec rr){
	Mat_SeqBSG *a = (Mat_SeqBSG *) A->data;
	const PetscScalar *l,*r, *li, *ri[a->stpoints];
	PetscInt dof = a->dof, bs = a->bs, m = a->nz*dof;
	PetscInt vm, stp, i, c, j, bc1, bc2, endpoint;
	PetscInt l3threshold = WORKINGSETSIZE/bs;
	PetscInt lda3 = a->m;
	PetscInt lda2 = lda3*a->n;
	PetscInt lda1 = lda2*a->p;

	PetscErrorCode ierr;

	PetscFunctionBegin;
	if(ll){
		ierr = VecGetArrayRead(ll,&l); CHKERRQ(ierr);
		ierr = VecGetLocalSize(ll,&vm); CHKERRQ(ierr);
		if (vm != m) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Left scaling vector wrong length");
		if(dof%2 == 0){
			for(i=0;i<a->nregion;i++){
				for(c = a->rstart[i]; c < a->rstart[i+1]; c++)
				{
						li = l+c*dof;
						for(stp = a->lbeg[i]; stp < a->lend[i]; stp++){
							for(bc1 = 0; bc1 < dof; bc1++){
								for(bc2 = 0; bc2 < dof; bc2++){
									a->coeff[i*a->stpoints+stp][(c-a->rstart[i])*bs+bc1*dof+bc2] *= li[bc1];	
								}
							}			
						} 
				}
			}
		}
		else{
			for(i=0;i<a->nregion;i++){
				for(c = a->rstart[i]; c < a->rstart[i+1]; c++)
				{
						li = l+c*dof;
						for(stp = a->lbeg[i]; stp < a->lend[i]; stp++){
							for(bc1 = 0; bc1 < dof; bc1++){
								for(bc2 = 0; bc2 < dof-1; bc2++){
									a->coeff[i*a->stpoints+stp][(c-a->rstart[i])*bs+bc1*(dof-1)+bc2] *= li[bc1];	
								}
							}			
							for(bc1 = 0; bc1 < dof; bc1++){
								a->coeff[i*a->stpoints+stp][(c-a->rstart[i])*bs+(dof-1)*dof+bc1] *= li[bc1];	
							}			
						} 
				}
			}
		}
		ierr = VecRestoreArrayRead(ll,&l); CHKERRQ(ierr);
    		ierr = PetscLogFlops(a->nz*a->stpoints*a->bs);CHKERRQ(ierr);
		
	}

	if(rr){
		ierr = VecGetArrayRead(rr,&l); CHKERRQ(ierr);
		ierr = VecGetLocalSize(rr,&vm); CHKERRQ(ierr);
		if (vm != m) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Left scaling vector wrong length");
		for(i=0;i<a->stpoints;i++){
			ri[i] = l + dof*( a->ioff[i]);
		}
		if(dof%2 == 0){
			for(i=0;i<a->nregion; i++){
				for(c = a->rstart[i] ; c< a->rstart[i+1] ; c++){
					for(stp = a->lbeg[i];stp<a->lend[i]; stp++){
						for(bc1 = 0; bc1 < dof; bc1++){
							for(bc2 = 0; bc2 < dof; bc2++){
								a->coeff[i*a->stpoints+stp][(c-a->rstart[i])*bs+bc1*dof+bc2]*= ri[stp][c*dof+bc2];
							}
						}
					} 
				}
			}
		}
		else{
			for(i=0;i<a->nregion; i++){
				for(c = a->rstart[i] ; c< a->rstart[i+1] ; c++){
					for(stp = a->lbeg[i];stp<a->lend[i]; stp++){
						for(bc1 = 0; bc1 < dof; bc1++){
							for(bc2 = 0; bc2 < dof-1; bc2++){
								a->coeff[i*a->stpoints+stp][(c-a->rstart[i])*bs+bc1*(dof-1)+bc2]*= ri[stp][c*dof+bc2];
							}
							a->coeff[i*a->stpoints+stp][(c-a->rstart[i])*bs+(dof-1)*dof+bc1]*= ri[stp][c*dof+dof-1];
						}
					} 
				}
			}
		}	
		ierr = VecRestoreArrayRead(rr,&l); CHKERRQ(ierr);
    		ierr = PetscLogFlops(a->nz*a->stpoints*a->bs);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatView_SeqBSG"
PetscErrorCode MatView_SeqBSG(Mat A,PetscViewer viewer)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * a = (Mat_SeqBSG *)A->data;
	PetscInt stpoints = a->stpoints, bs= a->bs;
	PetscInt stcount, icount, jcount;
	PetscScalar ** data = a->coeff;
	PetscInt start[8], count;
	PetscFunctionBegin;

	ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
	ierr = PetscObjectPrintClassNamePrefixType((PetscObject)A,viewer,"Matrix Object");CHKERRQ(ierr);
	
/*	start[0] = 0; 
	start[1] = 1; 
	start[2] = a->m; 
	start[3] = a->m*a->n; 
	start[4] = a->nz - (a->m*a->n); 
	start[5] = a->nz - a->m; 
	start[6] = a->nz - 1; 
	start[7] = a->nz; 
*/
	for(count = 0; count < a->nregion; count ++)
	{
		for(icount = a->rstart[count]; icount < a->rstart[count+1] ; icount++)
		{
			ierr = PetscViewerASCIIPrintf(viewer,"\n%D : [ ",icount);CHKERRQ(ierr);
			//for(stcount = a->lbeg[count] ; stcount < a->lend[count]; stcount++){
			for(stcount = 0 ; stcount < a->stpoints; stcount++){

				ierr = PetscViewerASCIIPrintf(viewer,"\t %D,( ",stcount);CHKERRQ(ierr);
				for(jcount = 0; jcount < bs ; jcount++)
				{
					ierr = PetscViewerASCIIPrintf(viewer,"%G  ",a->coeff[count*a->stpoints + stcount][icount-a->rstart[count]+jcount]);CHKERRQ(ierr);
				}
				ierr = PetscViewerASCIIPrintf(viewer,") \n");CHKERRQ(ierr);
			}
			ierr = PetscViewerASCIIPrintf(viewer," ]\n");CHKERRQ(ierr);
		}
	} 
	

	ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}






#ifdef _VEC2
/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix with the 3d indices supplied */
PetscErrorCode SetValues_Matrix_SeqBSG(Mat_SeqBSG *  mat, PetscInt n , const PetscInt ipos[],const PetscInt icols[], const PetscInt ioffsets[], const  PetscScalar data[], InsertMode is)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,c,endpoint,k1,mx = mat->m, ny= mat->n, dof = mat->dof, bs = mat->bs;
	PetscInt pos, l3threshold = WORKINGSETSIZE/bs;
	PetscInt rdis, cdis;
	PetscInt *start, *lbeg, *lend, finalpos, nregion, ioff;
	//PetscInt count
	l3threshold = WORKINGSETSIZE/mat->bs;
		start = mat->rstart;
		lbeg = mat->lbeg; 
		lend = mat->lend; 
		nregion = mat->nregion;

	for(i=0;i<n;i++)
	{
		pos = ipos[i]/dof;
		rdis = ipos[i] %dof;
		cdis = icols[i] %dof;
		//count = 0;
		if(ioffsets[i] == 3){
			//diag
			if(is == ADD_VALUES)
			{
					mat->diag[pos*bs+rdis*dof+cdis ] += data[i];	
			} else {
					mat->diag[pos*bs+rdis*dof+cdis ] = data[i];	
			}
		}
		for(k=0;k<nregion;k++)
		{
			if(start[k] <= pos && pos < start[k+1])
			{
				finalpos = pos - start[k];
				//ioff = count + ioffsets[i];// for submatrix - lbeg[k];
				ioff = k*mat->stpoints + ioffsets[i];// for submatrix - lbeg[k];
				if(dof % 2 == 0)
				{
					if(is == ADD_VALUES)
					{
							mat->coeff[ioff][(finalpos)*bs+rdis*dof+cdis] += data[i];
					}
					else
					{
							mat->coeff[ioff][(finalpos)*bs+rdis*dof+cdis] = data[i];
							
					}
				}
				else
				{
					if(is == ADD_VALUES)
					{
						if(cdis != dof-1)
							mat->coeff[ioff][(finalpos)*bs+rdis*(dof-1)+cdis] += data[i];
						else
							mat->coeff[ioff][(finalpos)*bs+dof*(dof-1)+rdis] += data[i];
					}
					else
					{
						if(cdis != dof-1)
							mat->coeff[ioff][(finalpos)*bs+rdis*(dof-1)+cdis] = data[i];
						else
							mat->coeff[ioff][(finalpos)*bs+dof*(dof-1)+rdis] = data[i];
					}
				}
				//count += mat->stpoints;
				//	count += (lend[k] - lbeg[k]);
			}
		}
	}
	PetscFunctionReturn(0);
}


/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix with the 3d indices supplied */
PetscErrorCode SetValuesBlocked_Matrix_SeqBSG(Mat_SeqBSG *  mat, PetscInt n , const PetscInt ipos[], const PetscInt ioffsets[], const  PetscScalar data[], InsertMode is)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,c,endpoint,k1,mx = mat->m, ny= mat->n, dof = mat->dof, bs = mat->bs;
	PetscInt pos, l3threshold = WORKINGSETSIZE/bs;
	PetscInt *start, *lbeg, *lend, finalpos, nregion, ioff;
	//PetscInt count
	l3threshold = WORKINGSETSIZE/mat->bs;
		start = mat->rstart;
		lbeg = mat->lbeg; 
		lend = mat->lend; 
		nregion = mat->nregion;

	for(i=0;i<n;i++)
	{
		pos = ipos[i];
		//count = 0;
		if(ioffsets[i] == 3){
			//diag
			if(is == ADD_VALUES)
			{
				for(j=0;j<bs;j++){
					mat->diag[pos*bs+j] += data[i*bs+j];	
				}
			} else {
				for(j=0;j<bs;j++){
					mat->diag[pos*bs+j] = data[i*bs+j];	
				}
			}
		}
		for(k=0;k<nregion;k++)
		{
			if(start[k] <= pos && pos < start[k+1])
			{
				finalpos = pos - start[k];
				//ioff = count + ioffsets[i];// for submatrix - lbeg[k];
				ioff = k*mat->stpoints + ioffsets[i];// for submatrix - lbeg[k];
				if(dof % 2 == 0)
				{
					if(is == ADD_VALUES)
					{
						for(j=0;j<dof;j++)
							for(k1=0;k1<dof;k1++)
								mat->coeff[ioff][(finalpos)*bs+j*dof+k1] += data[i*bs+j*dof+k1];
					}
					else
					{
						for(j=0;j<dof;j++)
							for(k1=0;k1<dof;k1++)
								mat->coeff[ioff][(finalpos)*bs+j*dof+k1] = data[i*bs+j*dof+k1];
							
					}
				}
				else
				{
					if(is == ADD_VALUES)
					{
						for(j=0;j<dof;j++)
						{
							for(k1=0;k1<dof-1;k1++)
								mat->coeff[ioff][(finalpos)*bs+j*(dof-1)+k1] += data[i*bs+j*dof+k1];
		
								mat->coeff[ioff][(finalpos)*bs+dof*(dof-1)+j] += data[i*bs+(j+1)*dof-1];
						}
					}
					else
					{
						for(j=0;j<dof;j++)
						{
							for(k1=0;k1<dof-1;k1++)
								mat->coeff[ioff][(finalpos)*bs+j*(dof-1)+k1] = data[i*bs+j*dof+k1];

								mat->coeff[ioff][(finalpos)*bs+dof*(dof-1)+j] = data[i*bs+(j+1)*dof-1];
						}
					}
				}
				//count += mat->stpoints;
				//	count += (lend[k] - lbeg[k]);
			}
		}
	}
	PetscFunctionReturn(0);
}
#endif



#ifdef _VEC4
/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix with the 3d indices supplied
 * Note: Not preferred for AVX */

PetscErrorCode SetValues_Matrix_SeqBSG(Mat_SeqBSG *  mat, PetscInt n , const PetscInt ipos[],const PetscInt icols[], const PetscInt ioffsets[], const  PetscScalar data[], InsertMode is)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,c,endpoint,k1,mx = mat->m, ny= mat->n, dof = mat->dof, bs = mat->bs;
	PetscInt pos, l3threshold = WORKINGSETSIZE/bs;
	PetscInt rdis, cdis, rdish, cdish;
	PetscInt *start, *lbeg, *lend, finalpos, nregion, ioff;
	PetscInt strtval;
	//PetscInt count
	l3threshold = WORKINGSETSIZE/mat->bs;
		start = mat->rstart;
		lbeg = mat->lbeg; 
		lend = mat->lend; 
		nregion = mat->nregion;

	for(i=0;i<n;i++)
	{
		pos = ipos[i]/dof;
		rdis = ipos[i] %dof;
		cdis = icols[i] %dof;
		//count = 0;
		if(ioffsets[i] == 3){
			//diag
			if(is == ADD_VALUES)
			{
					mat->diag[pos*bs+rdis*dof+cdis ] += data[i];	
			} else {
					mat->diag[pos*bs+rdis*dof+cdis ] = data[i];	
			}
		}
		for(k=0;k<nregion;k++)
		{
			if(start[k] <= pos && pos < start[k+1])
			{
				finalpos = pos - start[k];
				//ioff = count + ioffsets[i];// for submatrix - lbeg[k];
				ioff = k*mat->stpoints + ioffsets[i];// for submatrix - lbeg[k];
                                if (is == ADD_VALUES) {
                                  mat->coeff[ioff][(finalpos)*bs+mat->block_arrangement[rdis*dof+cdis]] += data[i];
                                } else {
                                  mat->coeff[ioff][(finalpos)*bs+mat->block_arrangement[rdis*dof+cdis]] = data[i];
                                }
                                  
			/*	if(rdis +4 <= dof){
					if(cdis +4 <= dof){
						strtval = rdis*4*dof + cdis*4*4; 
						if(is == ADD_VALUES)
						{
							if(cdish < 2){
								if(rdish == 0){
									mat->coeff[ioff][(finalpos)*bs+strtval+cdis+cdish] += data[i];
								}
								else if(rdish == 1){
									mat->coeff[ioff][(finalpos)*bs+strtval+8+cdis+cdish] += data[i];
								}
								else if(rdish == 2){
									mat->coeff[ioff][(finalpos)*bs+strtval+4+cdis+cdish+2] += data[i];
								}
								else{
									mat->coeff[ioff][(finalpos)*bs+strtval+12+cdis+cdish+2] += data[i];
								}
							}
							else{
								if(rdish == 0){
									mat->coeff[ioff][(finalpos)*bs+strtval+4+cdis+cdish-2] += data[i];
								}
								else if(rdish == 1){
									mat->coeff[ioff][(finalpos)*bs+strtval+12+cdis+cdish-2] += data[i];
								}
								else if(rdish == 2){
									mat->coeff[ioff][(finalpos)*bs+strtval+cdis+cdish] += data[i];
								}
								else{
									mat->coeff[ioff][(finalpos)*bs+strtval+8+cdis+cdish] += data[i];
								}
							}	
						}
						else
						{
							if(cdish < 2){
								if(rdish == 0){
									mat->coeff[ioff][(finalpos)*bs+strtval+cdis+cdish] = data[i];
								}
								else if(rdish == 1){
									mat->coeff[ioff][(finalpos)*bs+strtval+8+cdis+cdish] = data[i];
								}
								else if(rdish == 2){
									mat->coeff[ioff][(finalpos)*bs+strtval+4+cdis+cdish+2] = data[i];
								}
								else{
									mat->coeff[ioff][(finalpos)*bs+strtval+12+cdis+cdish+2] = data[i];
								}
							}
							else{
								if(rdish == 0){
									mat->coeff[ioff][(finalpos)*bs+strtval+4+cdis+cdish-2] = data[i];
								}
								else if(rdish == 1){
									mat->coeff[ioff][(finalpos)*bs+strtval+12+cdis+cdish-2] = data[i];
								}
								else if(rdish == 2){
									mat->coeff[ioff][(finalpos)*bs+strtval+cdis+cdish] = data[i];
								}
								else{
									mat->coeff[ioff][(finalpos)*bs+strtval+8+cdis+cdish] = data[i];
								}
							}	
						}
					}
				}*/
				//count += mat->stpoints;
				//	count += (lend[k] - lbeg[k]);
			}
		}
	}
	PetscFunctionReturn(0);
}

/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix with the 3d indices supplied
Added by Deepan */
PetscErrorCode SetValuesBlocked_Matrix_SeqBSG(Mat_SeqBSG *  mat, PetscInt n , const PetscInt ipos[], const PetscInt ioffsets[], const  PetscScalar data[], InsertMode is)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,c,endpoint,k1,k2,mx = mat->m, ny= mat->n, dof = mat->dof, bs = mat->bs;
	PetscInt pos, l3threshold = WORKINGSETSIZE/bs;
	PetscInt *start, *lbeg, *lend, finalpos, nregion, ioff, strtval;
	//PetscInt count
	l3threshold = WORKINGSETSIZE/mat->bs;
		start = mat->rstart;
		lbeg = mat->lbeg; 
		lend = mat->lend; 
		nregion = mat->nregion;

	for(i=0;i<n;i++)
	{
		pos = ipos[i];
		//count = 0;
		if(ioffsets[i] == 3){
			//diag
			if(is == ADD_VALUES)
			{
				for(j=0;j<bs;j++){
					mat->diag[pos*bs+j] += data[i*bs+j];	
				}
			} else {
				for(j=0;j<bs;j++){
					mat->diag[pos*bs+j] = data[i*bs+j];	
				}
			}
		}
		for(k=0;k<nregion;k++)
		{
			if(start[k] <= pos && pos < start[k+1])
			{
				finalpos = pos - start[k];
				//ioff = count + ioffsets[i];// for submatrix - lbeg[k];
				ioff = k*mat->stpoints + ioffsets[i];// for submatrix - lbeg[k];
				if(is == ADD_VALUES)
				{
					/*for(j=0;(j+4)<=dof;j+=4){
						for(k1=0;(k1+4)<=dof;k1+=4){
							strtval = j*4*dof + k1*4*4; 
							for(k2=0; k2<2; k2++){
								mat->coeff[ioff][(finalpos)*bs+strtval+k2] += data[i*bs+j*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+4+k2] += data[i*bs+j*dof+k1+k2+2];
								mat->coeff[ioff][(finalpos)*bs+strtval+8+k2] += data[i*bs+(j+1)*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+12+k2] += data[i*bs+(j+1)*dof+k1+k2+2];
							}
							for(k2=2; k2<4; k2++){
								mat->coeff[ioff][(finalpos)*bs+strtval+k2] += data[i*bs+(j+2)*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+4+k2] += data[i*bs+(j+2)*dof+k1+k2-2];
								mat->coeff[ioff][(finalpos)*bs+strtval+8+k2] += data[i*bs+(j+3)*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+12+k2] += data[i*bs+(j+3)*dof+k1+k2-2];
							}
						}
					}*/
					for(j = 0; j< bs ; j++)
						mat->coeff[ioff][(finalpos)*bs+mat->block_arrangement[j]] += data[i*bs+j];
						
				}
				else
				{
					/*for(j=0;(j+4)<=dof;j+=4){
						for(k1=0;(k1+4)<=dof;k1+=4){
							strtval = j*4*dof + k1*4*4; 
							for(k2=0; k2<2; k2++){
								mat->coeff[ioff][(finalpos)*bs+strtval+k2] = data[i*bs+j*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+4+k2] = data[i*bs+j*dof+k1+k2+2];
								mat->coeff[ioff][(finalpos)*bs+strtval+8+k2] = data[i*bs+(j+1)*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+12+k2] = data[i*bs+(j+1)*dof+k1+k2+2];
							}
							for(k2=2; k2<4; k2++){
								mat->coeff[ioff][(finalpos)*bs+strtval+k2] = data[i*bs+(j+2)*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+4+k2] = data[i*bs+(j+2)*dof+k1+k2-2];
								mat->coeff[ioff][(finalpos)*bs+strtval+8+k2] = data[i*bs+(j+3)*dof+k1+k2];
								mat->coeff[ioff][(finalpos)*bs+strtval+12+k2] = data[i*bs+(j+3)*dof+k1+k2-2];
							}
						}
					}*/
					for(j = 0; j< bs ; j++)
						mat->coeff[ioff][(finalpos)*bs+mat->block_arrangement[j]] = data[i*bs+j];
				}
				//count += mat->stpoints;
				//	count += (lend[k] - lbeg[k]);
			}
		}
	}
	PetscFunctionReturn(0);
}
#endif



/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix with the 3d indices supplied
Added by Deepan *//*
PetscErrorCode SetValues_Matrix_SeqBSG(Mat_SeqBSG *  mat, PetscInt n , const PetscInt idx[], const PetscInt idy[],const PetscInt idz[], const PetscInt ioffsets[], const  PetscScalar data[], InsertMode is)
{
	PetscInt i,j,k,k1,c,mx = mat->m, ny= mat->n, dof = mat->dof, bs = mat->bs;
	PetscInt lda1 = mx*ny, lda2 = mx;
	PetscInt pos, l3threshold = WORKINGSETSIZE/bs;
	PetscInt start[8], count, endpoint, finalpos;
	start[0] = 0; 
	start[1] = 1; 
	start[2] = mat->m; 
	start[3] = mat->m*mat->n; 
	start[4] = mat->nz - (mat->m*mat->n); 
	start[5] = mat->nz - mat->m; 
	start[6] = mat->nz - 1; 
	start[7] = mat->nz; 
	l3threshold = WORKINGSETSIZE/mat->bs;

	for(i=0;i<n;i++)
	{
		pos = (lda1*idz[i]+lda2*idy[i]+idx[i]);
		count = 0;
		for(k=0;k<7;k++)
		{
			for(c = start[k]; c < start[k+1]; c+= l3threshold,count += mat->stpoints )
			{
				endpoint = (c+l3threshold) < start[k+1] ? (c+l3threshold) : start[k+1]; 
				if(c <= pos && pos < endpoint)
				{
					finalpos = pos - c;
					if(dof % 2 == 0)
					{
						if(is == ADD_VALUES)
						{
							for(j=0;j<dof;j++)
								for(k1=0;k1<dof;k1++)
									mat->coeff[ioffsets[i]+count][(finalpos)*bs+j*dof+k1] += data[i*bs+j*dof+k1];
						}
						else
						{
							for(j=0;j<dof;j++)
								for(k1=0;k1<dof;k1++)
									mat->coeff[ioffsets[i]+count][(finalpos)*bs+j*dof+k1] = data[i*bs+j*dof+k1];
						}
					}
					else
					{
						if(is == ADD_VALUES)
						{
							for(j=0;j<dof;j++)
							{
								for(k1=0;k1<dof-1;k1++)
									mat->coeff[ioffsets[i]+count][(finalpos)*bs+j*(dof-1)+k1] += data[i*bs+j*dof+k1];
			
									mat->coeff[ioffsets[i]+count][(finalpos)*bs+dof*(dof-1)+j] += data[i*bs+(j+1)*dof-1];
							}
						}
						else
						{
							for(j=0;j<dof;j++)
							{
								for(k1=0;k1<dof-1;k1++)
									mat->coeff[ioffsets[i]+count][(finalpos)*bs+j*(dof-1)+k1] = data[i*bs+j*dof+k1];
			
									mat->coeff[ioffsets[i]+count][(finalpos)*bs+dof*(dof-1)+j] = data[i*bs+(j+1)*dof-1];
							}
						}
					}
	 
				}
				
			}
		}
	}
	PetscFunctionReturn(0);
}
*/
