#include <../src/mat/impls/sggpu/mpi/mpisggpu.h>

// Matrix function table
static struct _MatOps MatOps_Values = {
/*0*/ MatSetValues_MPISGGPU,MatGetRow_MPISGGPU,MatRestoreRow_MPISGGPU,MatMult_MPISGGPU,0,
/*5*/0,0,0,0,0,
/*10*/0,0,0,0,0,
/*15*/0,0,MatGetDiagonal_MPISGGPU,MatDiagonalScale_MPISGGPU,0,
/*20*/MatAssemblyBegin_MPISGGPU,MatAssemblyEnd_MPISGGPU,0,MatZeroEntries_MPISGGPU,0,
/*25*/0,0,0,0,MatSetUp_MPISGGPU,
/*30*/0,0,0,0,0,
/*35*/0,0,0,0,0,
/*40*/0,0,0,0,0,
/*45*/0,0,0,0,0,
/*50*/0,0,MatGetColumnIJ_MPISGGPU,0,MatFDColoringCreate_MPISGGPU,
/*55*/0,0,0,MatSetValuesBlocked_MPISGGPU,0,
/*60*/MatDestroy_MPISGGPU,MatView_MPISGGPU,0,0,0,
/*65*/0,0,MatSetValues_MPISGGPU,0,MatGetRowMaxAbs_MPISGGPU,
/*70*/0,0,0,0,0,
/*75*/MatFDColoringApply_MPISGGPU,0,0,0,0,
/*80*/0,0,0,0,0,
/*85*/0,0,MatSetValuesBlocked_MPISGGPU,0,0,
/*90*/0,0,0,0,0,
/*95*/0,0,0,0,0,
/*100*/0,0,0,0,0,
/*105*/0,0,0,0,0,
/*110*/0,0,0,0,0,
/*115*/MatCreate_MPISGGPU,0,0,0,0,
/*120*/0,0,0,0,0,
/*125*/0,0,0,0,0,
/*130*/0,0,0,0,0,
/*135*/0,0,0,0,0,
/*140*/0,0,
/*142*/MatSetGrid_MPISGGPU
};



EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatCreate_MPISGGPU"
PetscErrorCode MatCreate_MPISGGPU(Mat A)
{
  Mat_SeqSGGPU * mat;
  PetscErrorCode ierr;
  PetscMPIInt size;

  PetscFunctionBegin;
  SGTrace;

//  ierr = MPI_Comm_size(((PetscObject)A)->comm, &size); CHKERRQ(ierr);
//  if (size > 1)
//    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Comm must be size 1");

  // Create internal matrix structure
  ierr = PetscMalloc(sizeof(Mat_SeqSGGPU), &mat); CHKERRQ(ierr);
  memset(mat, 0, sizeof(Mat_SeqSGGPU));
  mat->diag_starts = new std::map<int, int>();
  mat->diagonals = new std::vector<int>();

  checkCudaError(cudaStreamCreate(&mat->stream));

  // Fill out PETSc matrix structure
  A->data = mat;
  memcpy(A->ops, &MatOps_Values, sizeof(struct _MatOps));
  A->same_nonzero= PETSC_FALSE;
  A->spptr = 0;

  // Set object type
  ierr = PetscObjectChangeTypeName((PetscObject)A, MATMPISGGPU); CHKERRQ(ierr);

  ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,
        "MatMPISGGPUSetPreallocation_C","MatMPISGGPUSetPreallocation_MPIDIA",
        MatMPISGGPUSetPreallocation_MPISGGPU);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
EXTERN_C_END




#undef __FUNCT__
#define __FUNCT__ "MatDestroy_MPISGGPU"
PetscErrorCode MatDestroy_MPISGGPU(Mat A) {
  MatDestroy_SeqSGGPU(A);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetGrid_MPISGGPU"
PetscErrorCode MatSetGrid_MPISGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p) 
{
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)B->data;

  PetscFunctionBegin;
  SGTrace;

  mat->m = m;
  mat->n = n;
  mat->p = p;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatMult_MPISGGPU"
PetscErrorCode MatMult_MPISGGPU(Mat mat, Vec x, Vec y) {

  MatView_MPISGGPU(mat,PETSC_VIEWER_STDOUT_WORLD);	

  MatMult_SeqSGGPU(mat, x, y); 

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_MPISGGPU"
PetscErrorCode MatSetValuesBlocked_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is) {
  MatSetValuesBlocked_SeqSGGPU(A, nrow, irow, ncol, icol, y, is);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValues_MPISGGPU"
PetscErrorCode MatSetValues_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is) {

  int i, j;
  PetscErrorCode ierr;
  PetscBool resizegpu = PETSC_FALSE;
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;

  PetscInt row, col;	
  PetscInt rank;
  PetscInt numprocs; 	

  PetscInt buf[8192],*bufr=0,*bufc=0,*irowm,*icolm;

  PetscFunctionBegin;
  SGTrace;
  
  PetscInt rstart = A->rmap->rstart, rend = A->rmap->rend;
  PetscInt cstart = A->cmap->rstart, cend = A->cmap->rend;  
 
    if ((nrow+ncol) <= (PetscInt)(sizeof(buf)/sizeof(PetscInt))) 
	{
	      irowm = buf; icolm = buf+nrow;
        } 
    else 
	{
	      ierr = PetscMalloc2(nrow,PetscInt,&bufr,ncol,PetscInt,&bufc);CHKERRQ(ierr);
	      irowm = bufr; icolm = bufc;
	}
   
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
   MPI_Comm_size(PETSC_COMM_WORLD, &numprocs);
   
   ierr = ISLocalToGlobalMappingApply(A->rmap->mapping,nrow,irow,irowm);CHKERRQ(ierr);
   ierr = ISLocalToGlobalMappingApply(A->cmap->mapping,ncol,icol,icolm);CHKERRQ(ierr);


  // Handle each element
  for (i = 0; i < nrow; i++) {
	if (irowm[i] < 0) continue;
	        row = irowm[i]; 

    for (j = 0; j < ncol; j++) {
  	
	if(irowm[i] >= rstart && irowm[i] < rend) 
	{
		if (icolm[j] >= cstart && icolm[j] < cend) 
		{
	        	col = icolm[j]; //-cstart;
		} 
		else if (icolm[j] < 0) 
			continue;
		else
			col = icolm[j];


      	      // Compute the diagonal and offset into the diagonal storage
	      // for the element
	      //int row = irow[i];
	      //int col = icol[j];
      
	      int diff = col - row;
	      int left = row % mat->dof;
	      int diag = int(floor((double)(diff + left) / mat->dof));
	      int col_offset = col % mat->dof;
	      int num_elems = mat->m * mat->n * mat->p * mat->dof;
	      int offset = col_offset * num_elems + row;

#if _TRACE
      printf("- row: %d  col: %d  val: %lf  diag: %d  offset: %d\n", row, col, y[i*ncol+j], diag, offset);
#endif

      std::map<int, int> &diag_starts = *(mat->diag_starts);
      std::map<int, int>::iterator I = diag_starts.find(diag);
      int diag_offset = 0;
      if (I == diag_starts.end()) {
        printf("WARNING: malloc() in MatSetValues\n");
        resizegpu = PETSC_TRUE;
        // The diagonal does not yet exist, so add a new diagonal
        int num_diags = diag_starts.size() + 1;
        int size = num_diags * mat->m * mat->n * mat->p * mat->dof * mat->dof;
        PetscScalar *newData;
        ierr = PetscMalloc(size * sizeof(PetscScalar), &newData); CHKERRQ(ierr);
        memset(newData, 0, size * sizeof(PetscScalar));
        size -= mat->m * mat->n * mat->p * mat->dof * mat->dof;
        if (num_diags > 1) {
          // This is not the first diagonal, so copy
#if _TRACE
          printf("- Memcpy of %d elements\n", size);
#endif
          memcpy(newData, mat->hostData, size * sizeof(PetscScalar));
        }
        PetscFree(mat->hostData);
        mat->hostData = newData;
        diag_offset = size;
        diag_starts[diag] = diag_offset;
        mat->diagonals->push_back(diag);
      } else {
        // The diagonal already exists, so get the base offset
        diag_offset = I->second;
      }

      diag_offset += offset;

      if (is == INSERT_VALUES)
        mat->hostData[diag_offset] = y[i * ncol + j];
      else
        mat->hostData[diag_offset] += y[i * ncol + j];

      mat->non_zeros++;
	}
    }
  }

  if (resizegpu) {
    int size,mat_size;
    // Create GPU buffer
    if (mat->deviceData) {
      cudaFree(mat->deviceData);
    }
    size = mat->diag_starts->size() * mat->m * mat->n * mat->p * mat->dof * mat->dof;
    checkCudaError(cudaMalloc(&mat->deviceData, sizeof(PetscScalar) * size));


    mat_size = mat->m * mat->n * mat->p * mat->dof;

    if (mat->deviceX) {
      cudaFree(mat->deviceX);
    }
    if (mat->deviceY) {
      cudaFree(mat->deviceY);
    }
    if (mat->deviceDiags) {
      cudaFree(mat->deviceDiags);
    }
    // We know the expected size of x, y, so go ahead and allocate them now
    checkCudaError(cudaMalloc(&mat->deviceX, mat_size * sizeof(PetscScalar)));
    checkCudaError(cudaMalloc(&mat->deviceY, mat_size * sizeof(PetscScalar)));

    // We also know how many diagonals we have, and their indices
    checkCudaError(cudaMalloc(&mat->deviceDiags, sizeof(int) * mat->diagonals->size()));
  }

  PetscFunctionReturn(0);

}


//#undef __FUNCT__
//#define __FUNCT__ "MatSetStencil_MPISGGPU"
//PetscErrorCode MatSetStencil_MPISGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof) {
//  MatSetStencil_SeqSGGPU(A, dim, dims, starts, dof);
//  PetscFunctionReturn(0);
//}


#undef __FUNCT__
#define __FUNCT__ "MatSetUp_MPISGGPU"
PetscErrorCode MatSetUp_MPISGGPU(Mat mat) {
  MatSetUp_SeqSGGPU(mat);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_MPISGGPU"
PetscErrorCode MatZeroEntries_MPISGGPU(Mat A) {
  MatZeroEntries_SeqSGGPU(A);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_MPISGGPU"
PetscErrorCode MatGetDiagonal_MPISGGPU(Mat A, Vec v) {
  MatGetDiagonal_SeqSGGPU(A, v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_MPISGGPU"
PetscErrorCode MatDiagonalScale_MPISGGPU(Mat A, Vec ll, Vec rr) {
  MatDiagonalScale_SeqSGGPU(A, ll, rr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRow_MPISGGPU"
PetscErrorCode MatGetRow_MPISGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v) {
  MatGetRow_SeqSGGPU(A, row, nz, idx , v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_MPISGGPU"
PetscErrorCode MatRestoreRow_MPISGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v) {
  MatRestoreRow_SeqSGGPU(A, row, nz, idx, v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_MPISGGPU"
PetscErrorCode MatGetRowMaxAbs_MPISGGPU(Mat A, Vec v, PetscInt idx[]) {
  MatGetRowMaxAbs_SeqSGGPU(A, v, idx);
  PetscFunctionReturn(0);
}


void DisplayLocalMatrix(Mat A);

#undef __FUNCT__
#define __FUNCT__ "MatView_MPISGGPU"
PetscErrorCode MatView_MPISGGPU(Mat A, PetscViewer viewer) 
{

  PetscInt rank, i, numprocs;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);	
	
  for(i = 0; i < numprocs;++i)
	{
	 if(i == rank)
	    DisplayLocalMatrix(A);
         MPI_Barrier(PETSC_COMM_WORLD);   
    	}
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyBegin_MPISGGPU"
PetscErrorCode MatAssemblyBegin_MPISGGPU(Mat A, MatAssemblyType type) {
  MatAssemblyBegin_SeqSGGPU(A, type);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyEnd_MPISGGPU"
PetscErrorCode MatAssemblyEnd_MPISGGPU(Mat A, MatAssemblyType type) {
  MatAssemblyEnd_SeqSGGPU(A, type);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatFDColoringApply_MPISGGPU"
PetscErrorCode  MatFDColoringApply_MPISGGPU(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx) {
  MatFDColoringApply_SeqSGGPU(J, coloring, x1, flag, sctx);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatFDColoringCreate_MPISGGPU"
PetscErrorCode MatFDColoringCreate_MPISGGPU(Mat mat,ISColoring iscoloring,MatFDColoring c) {
  MatFDColoringCreate_SeqSGGPU(mat, iscoloring, c);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetColumnIJ_MPISGGPU"
PetscErrorCode MatGetColumnIJ_MPISGGPU(Mat A,PetscInt oshift,PetscBool  symmetric,PetscBool  inodecompressed,PetscInt *nn, const PetscInt *ia[], const PetscInt *ja[],PetscBool  *done) {
  MatGetColumnIJ_SeqSGGPU(A, oshift, symmetric, inodecompressed, nn, ia, ja, done);
  PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatMPISGGPUSetPreallocation"
PetscErrorCode MatMPISGGPUSetPreallocation(Mat A,PetscInt stencil_type, PetscInt dof)
{
  PetscErrorCode ierr;
  Mat_SeqSGGPU *mat = (Mat_SeqSGGPU*)A->data;

  PetscFunctionBegin;

  mat->stencil_type = stencil_type;
  mat->dof = dof;
  if(A->preallocated)PetscFunctionReturn(0);
  PetscValidHeaderSpecific(A,MAT_CLASSID,1);
  
  ierr = PetscTryMethod(A,"MatMPISGGPUSetPreallocation_C",(Mat,PetscInt,const PetscInt []),(A,0,0));CHKERRQ(ierr);
  A->preallocated=PETSC_TRUE;
  PetscFunctionReturn(0);
}
EXTERN_C_END


EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatMPISGGPUSetPreallocation_MPISGGPU"
extern PetscErrorCode MatMPISGGPUSetPreallocation_MPISGGPU(Mat A,PetscInt nz, const PetscInt nnz[])
{
  PetscErrorCode ierr;
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;

  PetscInt dim,diag_size,size,num_diags,i,vecsize;

  ierr = PetscLayoutSetBlockSize(A->rmap,1);CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(A->cmap,1);CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(A->rmap);CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(A->cmap);CHKERRQ(ierr);

  dim = A->stencil.dim;
  if (mat->dof > 1) {
    dim--;
  }

//  mat->m = mat->n = mat->p = 1;
  mat->dim = dim;
//  if (mat->dim > 0) mat->m = A->stencil.dims[dim-1];
//  if (mat->dim > 1) mat->n = A->stencil.dims[dim-2];
//  if (mat->dim > 2) mat->p = A->stencil.dims[dim-3];

  if (mat->stencil_type == 0) {
    /* star stencil */
    num_diags = 2*mat->dim + 1;
  } else {
    /* box stencil */
    num_diags =  1;
    for (i=0;i<mat->dim;i++) num_diags*=3;
  }

  diag_size = mat->m * mat->n * mat->p * mat->dof * mat->dof;
  size = num_diags * diag_size;

  if (mat->m == 0 || mat->n == 0 || mat->p == 0 || mat->dof == 0) {
    SETERRQ(PETSC_COMM_SELF,0,"MatSetPreallocation_SeqSGGPU called without valid m, n, p, and dof!");
  }


  ierr = PetscMalloc(sizeof(PetscInt)*num_diags,&mat->diag_offsets);
  ierr = PetscMalloc(size * sizeof(PetscScalar), &mat->hostData); CHKERRQ(ierr);
  memset(mat->hostData, 0, size * sizeof(PetscScalar));

  (*mat->diag_starts)[0]  = 0 * diag_size;
  (*mat->diagonals).push_back(0);
  (*mat->diag_starts)[1]  = 1 * diag_size;
  (*mat->diagonals).push_back(1);
  (*mat->diag_starts)[-1] = 2 * diag_size;
  (*mat->diagonals).push_back(-1);
  if (mat->stencil_type == 0) {
    if (mat->dim == 2) {
      (*mat->diag_starts)[mat->m] = 3 * diag_size;
      (*mat->diagonals).push_back(mat->m);
      (*mat->diag_starts)[-mat->m] = 4 * diag_size;
      (*mat->diagonals).push_back(-mat->m);
    } else if (mat->dim == 3) {
      (*mat->diag_starts)[mat->m] = 3 * diag_size;
      (*mat->diagonals).push_back(mat->m);
      (*mat->diag_starts)[-mat->m] = 4 * diag_size;
      (*mat->diagonals).push_back(-mat->m);

      (*mat->diag_starts)[mat->m*mat->n] = 5 * diag_size;
      (*mat->diagonals).push_back(mat->m*mat->n);
      (*mat->diag_starts)[-mat->m*mat->n] = 6 * diag_size;
      (*mat->diagonals).push_back(-mat->m*mat->n);
    }
  } else {
    if (mat->dim == 2) {
      (*mat->diag_starts)[mat->n-1] = 3 * diag_size;
      (*mat->diagonals).push_back(mat->m);
      (*mat->diag_starts)[-mat->n-1] = 4 * diag_size;
      (*mat->diagonals).push_back(-mat->m);
      (*mat->diag_starts)[mat->n] = 5 * diag_size;
      (*mat->diagonals).push_back(mat->m);
      (*mat->diag_starts)[-mat->n] = 6 * diag_size;
      (*mat->diagonals).push_back(-mat->m);
      (*mat->diag_starts)[mat->n+1] = 7 * diag_size;
      (*mat->diagonals).push_back(mat->m);
      (*mat->diag_starts)[-mat->n+1] = 8 * diag_size;
      (*mat->diagonals).push_back(-mat->m);
    }
  }
  /*
  printf("Diagonals preallocated:\n");
  for (std::map<int, int>::iterator I = mat->diag_starts->begin(),
         E = mat->diag_starts->end(); I != E; ++I) {
    printf("%4d --> %4d\n",I->first,I->second);
  }
   */
  
  
  // Create GPU buffer
  if (mat->deviceData) {
    cudaFree(mat->deviceData);
  }
  checkCudaError(cudaMalloc(&mat->deviceData, sizeof(PetscScalar) * size));
  checkCudaError(cudaMemset(mat->deviceData,0.0,sizeof(PetscScalar)*size));

  // Copy data to device
  checkCudaError(cudaMemcpy(mat->deviceData, mat->hostData, sizeof(PetscScalar) * size, cudaMemcpyHostToDevice));


  vecsize = mat->m * mat->n * mat->p * mat->dof;

  // We know the expected size of x, y, so go ahead and allocate them now
  checkCudaError(cudaMalloc(&mat->deviceX, vecsize * sizeof(PetscScalar)));
  checkCudaError(cudaMalloc(&mat->deviceY, vecsize * sizeof(PetscScalar)));

  // We also know how many diagonals we have, and their indices
  checkCudaError(cudaMalloc(&mat->deviceDiags, sizeof(int) * mat->diagonals->size()));
  A->preallocated = PETSC_TRUE;
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  

  PetscFunctionReturn(0);
}
EXTERN_C_END



void DisplayLocalMatrix(Mat A)
{

  Mat_SeqSGGPU *a;
  a  = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode ierr;
  PetscInt nrows,ndiag,dof,i,j,iblock,col,index,offset;
  std::map<int, int> &diag_starts = *(a->diag_starts);
  
  nrows = a->m * a->n * a->p * a->dof;
  ndiag = a->diagonals->size();
  dof = a->dof;

  for (std::map<int, int>::iterator I = diag_starts.begin(),
         E = diag_starts.end(); I != E; ++I) 
	{
	fprintf(stdout,"- Diag %d:%d\n", I->first, I->second);
  	}

	fprintf(stdout,"\n");
	fprintf(stdout,"hostData:\n");

  for (i=0;i<nrows;i++) 
	{
	fprintf(stdout,"row %2.2d:",i); 

	    for (j=0;j<ndiag*dof;j++) 
		{
			fprintf(stdout," %4g ",a->hostData[i+j*nrows]);
		}
	fprintf(stdout,"\n");
  	}
	fprintf(stdout,"\n\n");


  for (iblock=0;iblock<nrows/dof;iblock++)  
	{
	    for (i=iblock*dof;i<(iblock+1)*dof;i++) 
		{
		fprintf(stdout,"row %d:",i);
	
		  for (std::map<int, int>::iterator I = a->diag_starts->begin(),
        	     E = a->diag_starts->end(); I != E; ++I) 
		  {
	        	/* Ignore 0 padding */
		        offset = I->first;
		        if (offset + iblock < 0) 
			{
	        	  continue;
        		}
	        	if (offset + iblock >= (nrows/dof)) 
			{
        	  	break;
        		}
        
	        	for (j=0;j<dof;j++) 
			{
	        	  index = i + I->second + j*nrows; // column-major
	        	  col = offset*dof+(iblock*dof) + j;
	                  fprintf(stdout," (%d, %g) ",col,a->hostData[index]);
        		}
      	 	  }
		fprintf(stdout,"\n");
    		}
  	}

}