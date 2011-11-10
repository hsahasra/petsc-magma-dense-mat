#define PETSCMAT_DLL


#include "../src/mat/impls/blockstructgrid/matblockstructgrid.h"
#include "petscblaslapack.h"
#include "petscbt.h"
#include "petscmat.h"
#include "../src/mat/impls/aij/seq/aij.h"  


//#include <immintrin.h>


static struct _MatOps MatOps_Values = {
/*0*/ MatSetValues_SeqBSG,0,0,MatMult_SeqBSG,0,
/*5*/0,0,0,0,0,
/*10*/0,0,0,0,0,
/*15*/0,0,0,0,0,
/*20*/0,0,0,MatZeroEntries_SeqBSG,0,
/*25*/0,0,0,0,MatSetUpPreallocation_SeqBSG,
/*30*/0,0,0,0,0,
/*35*/0,0,0,0,0,
/*40*/0,0,0,0,0,
/*45*/0,0,0,0,0,
/*50*/0,0,0,0,0,
/*55*/0,0,0,MatSetValuesBlocked_SeqBSG,0,
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


/** MatSetGrid_SeqBSG : Sets the 3d physical grid information
Added by Deepan */
#undef __FUNCT__ 
#define __FUNCT__ "MatSetGrid_SeqBSG"

PetscErrorCode MatSetGrid_SeqBSG(Mat B, PetscInt m, PetscInt n, PetscInt p)
{
	Mat_SeqBSG * b = (Mat_SeqBSG*) B->data;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	if(m <= 0 || n <= 0 || p <= 0 ) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Grid Dimension should be atleast 1");

	b->m = m;
	b->n = n;
	b->p = p;
	b->nz = m*n*p;
	PetscFunctionReturn(0);
}

/** MatCreate_SeqBSG : Creates the struct grid representation 
Added by Deepan */
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
	b->a = 0;
	b->m = 0;
	b->n = 0;
	b->p = 0;
	b->dof = 0;
	b->nz = 0;
	b->idx = PETSC_NULL;
	b->idy = PETSC_NULL;
	b->idz = PETSC_NULL;

	b->tnz = 0;
	b->bs = 0;
	ierr = PetscObjectChangeTypeName((PetscObject)B, MATBLOCKSTRUCTGRID); CHKERRQ(ierr);
	PetscFunctionReturn(0);	
}

/** MatDestroy_SeqBSG : Destroys the struct grid representation
Added by Deepan */
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
	ierr = PetscFree(a->a);CHKERRQ(ierr);
	ierr = PetscFree(a->idx);CHKERRQ(ierr);
	ierr = PetscFree(a->idy);CHKERRQ(ierr);
	ierr = PetscFree(a->idz);CHKERRQ(ierr);

	ierr = PetscFree(a);CHKERRQ(ierr);

	ierr = PetscObjectChangeTypeName((PetscObject)A, 0);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/** MatMult_SeqBSG : Performs the Matrix - Vector multiplication y= mat*x on the struct grid representation
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqBSG"
PetscErrorCode MatMult_SeqBSG(Mat mat, Vec x, Vec y)
{
	PetscErrorCode ierr;
	PetscInt size, i;
	Mat_SeqBSG * a = (Mat_SeqBSG *) mat->data;
	PetscScalar * v = a->a, *xx,*yy;
	
	PetscFunctionBegin;
	ierr = VecSet(y,0.0); CHKERRQ(ierr);
	ierr = VecGetArray(x, &xx); CHKERRQ(ierr);
	ierr = VecGetArray(y, &yy); CHKERRQ(ierr);
	ierr = BSG_MatMult(v, xx, yy, a->idx, a->idy, a->idz, a->m, a->n, a->p, a->dof, a->stpoints, 3, a->bs);
	CHKERRQ(ierr);

	ierr = VecRestoreArray(x,&xx); CHKERRQ(ierr);
	ierr = VecRestoreArray(y,&yy); CHKERRQ(ierr);
	ierr = PetscLogFlops(2*a->nz*a->stpoints); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_SeqBSG"

PetscErrorCode MatSetValuesBlocked_SeqBSG(Mat A, PetscInt nrow,const PetscInt irow[], PetscInt ncol,const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
	PetscErrorCode ierr;
	ierr = MatSetValues_SeqBSG(A,nrow,irow,ncol,icol,y,is); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}	

/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetValues_SeqBSG"

PetscErrorCode MatSetValues_SeqBSG(Mat A, PetscInt nrow,const PetscInt irow[], PetscInt ncol ,const PetscInt icol[],const PetscScalar y[], InsertMode is)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * mat = (Mat_SeqBSG *) A->data;
	PetscInt * idx, * idy, * idz, *ioffsets;
	PetscInt i,j,count = 0, m,n,p, offset,dis,xdis,ydis,zdis, cdis, k ,stp,dof,rshift,cshift, rowval;
	PetscFunctionBegin;	
	
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&idx);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&idy);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&idz);CHKERRQ(ierr);
	ierr = PetscMalloc(nrow*ncol*sizeof(PetscInt),&ioffsets);CHKERRQ(ierr);
	
	m = mat->m;
	n = mat->n;
	p = mat->p;
	stp = mat->stpoints;
	dis = mat->dis;
	dof = mat->dof;
	
	for(i=0;i< nrow ; i++)
	{
		rshift = irow[i];
		for(j=0;j<ncol;j++)
		{
			cshift = icol[j];	
			cdis = cshift - rshift;
			zdis = (cdis)/(m*n);
			cdis = (cdis)%(m*n);
			ydis = (cdis)/(m);
			cdis = (cdis)%(m);
			xdis = cdis;
			for(k=0;k<stp;k++)
				if(mat->idx[k] == xdis &&  mat->idy[k] == ydis &&  mat->idz[k] == zdis)	
				{
					offset = k; //corresponding band
					break;
				}
			if(k<stp)
			{
				rowval = rshift;
				ioffsets[count] = k;
				idx[count] = rowval % (m);
				idy[count] = (rowval/(m))%n;
				idz[count] = ((rowval/(m*n)) %p) ;
				count ++;
			}
				//printf("\n x:%d y:%d ix:%d iy:%d iz:%d k:%d val:%d",irow[i],icol[j],idx[count],idy[count],idz[count],ioffsets[count],(stp*mat->bs*(mat->m*mat->n*idz[count]+mat->m*idy[count]+idx[count]))+(ioffsets[count]*mat->bs));
		}	
	}
	ierr = SetValues_SeqBSG(mat,nrow*ncol,idx,idy,idz,ioffsets,y,is); CHKERRQ(ierr);
	ierr = PetscFree(ioffsets);CHKERRQ(ierr);
	ierr = PetscFree(idx);CHKERRQ(ierr);
	ierr = PetscFree(idy);CHKERRQ(ierr);
	ierr = PetscFree(idz);CHKERRQ(ierr);

	if (A->valid_GPU_matrix != PETSC_CUSP_UNALLOCATED)
	    A->valid_GPU_matrix = PETSC_CUSP_CPU;
	PetscFunctionReturn(0);
}

/** MatSetValuesBlocked_SeqBSG : Sets the values in the matrix with the 3d indices supplied
Added by Deepan */
PetscErrorCode SetValues_SeqBSG(Mat_SeqBSG *  mat, PetscInt n , const PetscInt idx[], const PetscInt idy[],const PetscInt idz[], const PetscInt ioffsets[], const  PetscScalar data[], InsertMode is)
{
	PetscInt i,j,mx = mat->m, ny= mat->n, pz = mat->p, dof = mat->dof, bs = mat->bs , nos = mat->stpoints;
	PetscInt lda1 = mx*ny, lda2 = mx, pos;
	if(is == ADD_VALUES)
	{
		for(i=0;i<n;i++)
		{
			pos = (nos*bs*(lda1*idz[i]+lda2*idy[i]+idx[i]))+(ioffsets[i]*bs);
			for(j=0;j<bs;j++)
				mat->a[pos + j] += data[i*bs+j];
		}
	}
	else
	{
		for(i=0;i<n;i++)
		{
			pos = (nos*bs*(lda1*idz[i]+lda2*idy[i]+idx[i]))+(ioffsets[i]*bs);
			for(j=0;j<bs;j++)
			{
				//printf("\n%d  %f", pos+j, data[i*bs+j]);
				mat->a[pos + j] = data[i*bs+j];
			}
		}
	}


	PetscFunctionReturn(0);
}


/** MatSetStencil_SeqBSG : Sets the stencil neighbor points
 * sets various stencil related info such as neighbor points displacements and dimension of the grid
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetStencil_SeqBSG"

PetscErrorCode MatSetStencil_SeqBSG(Mat A, PetscInt dim,const PetscInt dims[],const PetscInt starts[], PetscInt dof)
{
	Mat_SeqBSG * mat = (Mat_SeqBSG *)A->data;
	PetscInt i,cnt=0, mnos = dim;
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
	mat->tnz = mat->stpoints*mat->nz;
	//printf("m=%d, n=%d,p=%d\n",mat->m,mat->n,mat->p);

	ierr = PetscMalloc(sizeof(PetscInt)*mat->stpoints,&(mat->idx));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*mat->stpoints,&(mat->idy));CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*mat->stpoints,&(mat->idz));CHKERRQ(ierr);
/*Do not change the order of the stencil. MatMult depends on the order of these stencils.*/	
	mat->idx[mnos] = 0; mat->idy[mnos]=0; mat->idz[mnos]= 0;
	if(dim > 0)
	{
		mat->idx[mnos+1] = 1; mat->idy[mnos+1] = 0; mat->idz[mnos+1] = 0; 
	}
	if(dim > 1)
	{
		mat->idx[mnos+2] = 0; mat->idy[mnos+2] = 1; mat->idz[mnos+2] = 0;
	}
	if(dim > 2)
	{
		mat->idx[mnos+3] = 0; mat->idy[mnos+3] = 0; mat->idz[mnos+3] = 1;
	}
	if(dim>0)
	{	
		mat->idx[mnos-1] = -1; mat->idy[mnos-1] = 0; mat->idz[mnos-1] = 0;
	}	
	if(dim>1)
	{
		mat->idx[mnos-2] = 0; mat->idy[mnos-2] = -1; mat->idz[mnos-2] = 0;
	}	
	if(dim>2)
	{
		mat->idx[mnos-3] = 0; mat->idy[mnos-3] = 0; mat->idz[mnos-3] = -1; 
	}	
 	PetscFunctionReturn(0);	
}

/** MatSetUpPreallocation_SeqBSG : Allocates space for coefficient matrix
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatSetUpPreallocation_SeqBSG"

PetscErrorCode MatSetUpPreallocation_SeqBSG(Mat mat)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	Mat_SeqBSG * a = (Mat_SeqBSG *)mat->data;
	PetscFunctionBegin;
	ierr = PetscMalloc(sizeof(PetscScalar)*a->tnz*a->bs,&(a->a));CHKERRQ(ierr);
	memset(a->a, 0,sizeof(PetscScalar)*a->tnz*a->bs);
	mat->preallocated = PETSC_TRUE;
	PetscFunctionReturn(0);
}

/** MatZeroEntries_SeqBSG : Sets the values in the matrix to be zero
Added by Deepan */
#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_SeqBSG"

PetscErrorCode MatZeroEntries_SeqBSG(Mat A)
{
	Mat_SeqBSG * a = (Mat_SeqBSG *)A->data;
	PetscFunctionBegin;
	memset(a->a,0,sizeof(PetscScalar)*a->tnz*a->bs);
	if (A->valid_GPU_matrix != PETSC_CUSP_UNALLOCATED)
	    A->valid_GPU_matrix = PETSC_CUSP_CPU;

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatView_SeqBSG"
PetscErrorCode MatView_SeqBSG(Mat A,PetscViewer viewer)
{
	PetscErrorCode ierr;
	Mat_SeqBSG * a = (Mat_SeqBSG *)A->data;
	PetscInt dof = a->dof, m = a->m , n = a->n , p = a->p, stpoints = a->stpoints, tnz = a->tnz, bs= a->bs;
	PetscInt stcount, icount, jcount, kcount;
	PetscScalar * data = a->a;
	PetscFunctionBegin;

	ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
	ierr = PetscObjectPrintClassNamePrefixType((PetscObject)A,viewer);CHKERRQ(ierr);

/*	for(stcount = 0; stcount < stpoints-1; stcount++)
	{
		ierr = PetscViewerASCIIPrintf(viewer,"\n\ndiagonal %D:\n",stcount);CHKERRQ(ierr);
		for(icount = coeffOffset[stcount]*bs; icount < coeffOffset[stcount+1]*bs; icount+=bs)
		{
			ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
			for(jcount = 0; jcount < bs; jcount +=dof)
			{
				ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
				for (kcount = 0 ; kcount < dof; kcount++)
				{
					ierr = PetscViewerASCIIPrintf(viewer," (%G) ",data[icount+jcount+kcount]);CHKERRQ(ierr);
				}
			}
		}
	}

	ierr = PetscViewerASCIIPrintf(viewer,"\n\ndiagonal %D:\n",stpoints-1);CHKERRQ(ierr);
	for(icount = coeffOffset[stpoints-1]*bs; icount < tnz*bs; icount+=bs)
	{
		ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
		for(jcount = 0; jcount < bs; jcount +=dof)
		{
			ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
			for (kcount = 0 ; kcount < dof; kcount++)
			{
				ierr = PetscViewerASCIIPrintf(viewer," (%G) ",data[icount+jcount+kcount]);CHKERRQ(ierr);
			}
		}
	}
*/
	for(stcount = 0; stcount < tnz*bs ; stcount++)
	{
		ierr = PetscViewerASCIIPrintf(viewer,"\n%D : (%G) ",stcount , data[stcount]);CHKERRQ(ierr);
	}
	ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
	ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

