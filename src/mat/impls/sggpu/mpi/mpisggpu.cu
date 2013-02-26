#include <../src/mat/impls/sggpu/mpi/mpisggpu.h>


// Direct access to seqgpu vector type
#include "../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"

// Interop with CUSP vector
#include "../src/vec/vec/impls/seq/seqcusp/cuspvecimpl.h"

#include "cuPrintf.cu"

#define BLOCKWIDTH_X 128
#define BLOCKWIDTH_Y 1



//===-- CUDA Device Code -------------------------------------------------===//
 
texture<int2, 1> vector_x;
     
static __inline__ __device__ double fetch_doubleMPI(texture<int2, 1> tex, int i)
     {
       int2 v = tex1Dfetch(tex, i);
       return __hiloint2double(v.y, v.x);
     }
     
__global__ void MatMultKernelMPI(PetscScalar * coeff, PetscScalar * y, PetscScalar *x,PetscInt mat_size, PetscInt num_diags, int * diagonals, PetscInt dof, PetscInt vec_size) {
       
int idx = blockDim.x * blockIdx.x * 1 + threadIdx.x * 1;

     if (idx >= mat_size) 
      return;
     
int diag_size = mat_size * dof;
     
PetscScalar yval0 = 0.0;
int idx0 = idx;


//#pragma unroll 4
for (int i = 0; i < num_diags; ++i) 
	{
	    int d = diagonals[i];
    
	    int offset0 = diag_size * i + idx0;
	    int block0 = (idx0/dof + d) * dof;
     
	    //#pragma unroll 12
	    for (int j = 0; j < dof; ++j) 
		{
		      // Get coefficient
		      PetscScalar aval0 = coeff[offset0 + mat_size*j];
		      // Get X value
    	              PetscScalar xval0 = 0;

	              if ((block0 >= 0) && (block0 < vec_size))
			//xval0 = x[block0 + j];
			xval0 = fetch_doubleMPI(vector_x, block0 + j);

		      yval0 += aval0 * xval0;
	
		      //21, 22, 25, 26, 37, 38, 41, 42 are the only non zero entries in the 64x1 vector for ex14 when grid size is 4x4x4 	
		      //The following if statement is for testing and debugging and can be removed. 
		      // if ((idx0 == 21) || (idx0 == 22) || (idx0 == 25) || (idx0 == 26) ||  (idx0 == 37) ||  (idx0 == 38) ||  (idx0 == 41) ||  (idx0 == 42))
		      //     cuPrintf("d:%d \t offset:%d \t block0:%d \t aval0:%lf \t xval0:%lf \t yval0 :%lf\n",d,offset0,block0,aval0,xval0, yval0);
	        }
	}
   
      y[idx0] = yval0;

}
    
//===-- Host Code --------------------------------------------------------===//


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
  Mat_MPISGGPU * mat;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  SGTrace;

  PetscInt rank, size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(((PetscObject)A)->comm, &size); CHKERRQ(ierr);

  // Create internal matrix structure
  ierr = PetscMalloc(sizeof(Mat_MPISGGPU), &mat); CHKERRQ(ierr);
  memset(mat, 0, sizeof(Mat_MPISGGPU));
  ierr = PetscMalloc(sizeof(Mat_SeqSGGPU), &mat->mat_seq); CHKERRQ(ierr);
  memset(mat->mat_seq, 0, sizeof(Mat_SeqSGGPU));
  mat->mat_seq->diag_starts = new std::map<int, int>();
  mat->mat_seq->diagonals = new std::vector<int>();
  mat->rank = rank;   
  mat->size = size;

  checkCudaError(cudaStreamCreate(&mat->mat_seq->stream));

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

  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;  
  PetscErrorCode ierr;

  PetscFunctionBegin;
 
  if (mat_seq->hostData) {
    ierr = PetscFree(mat_seq->hostData); CHKERRQ(ierr);
  }
  if (mat_seq->deviceData) {
    cudaFree(mat_seq->deviceData);
  }
  if (mat_seq->diag_starts) {
    delete mat_seq->diag_starts;
  }
  ierr = PetscFree(mat_seq->diag_offsets); CHKERRQ(ierr);
  if (mat_seq->diagonals) {
    delete mat_seq->diagonals;
  }
  if (mat_seq->deviceX) {
    cudaFree(mat_seq->deviceX);
  }
  if (mat_seq->deviceY) {
    cudaFree(mat_seq->deviceY);
  }
  if (mat_seq->deviceDiags) {
    cudaFree(mat_seq->deviceDiags);
  }
  
  if(mat_seq->ja) { ierr = PetscFree(mat_seq->ja); CHKERRQ(ierr); }
  
  if(mat_seq->ia) { ierr = PetscFree(mat_seq->ia); CHKERRQ(ierr); }
  
  checkCudaError(cudaStreamDestroy(mat_seq->stream));
  if(mat_seq)
	  ierr = PetscFree(mat_seq); CHKERRQ(ierr);
    
  ierr = VecDestroy(&mat->lvec); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&mat->Mvctx); CHKERRQ(ierr);
  ierr = PetscFree(A->data); CHKERRQ(ierr);
  
  ierr = PetscObjectChangeTypeName((PetscObject)A, 0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetGrid_MPISGGPU"
PetscErrorCode MatSetGrid_MPISGGPU(Mat A, PetscInt m, PetscInt n, PetscInt p) 
{
  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;

  PetscFunctionBegin;
  SGTrace;

  mat_seq->m = m;
  mat_seq->n = n;
  mat_seq->p = p;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatMult_MPISGGPU"
PetscErrorCode MatMult_MPISGGPU(Mat A, Vec x, Vec y) {

  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;

//  MatView_MPISGGPU(A,PETSC_VIEWER_STDOUT_WORLD);	

  PetscBool isseqcusp,isseqgpu,ismpicusp,iscusp;
  PetscErrorCode ierr;
  PetscInt mat_size, vec_size;
  CUSPARRAY *xgpu,*ygpu;
  PetscScalar *devX,*devY;

  PetscFunctionBegin;
  SGTrace;

  // Initialize y to zero
  ierr = VecSet(y, 0.0); CHKERRQ(ierr);

  ierr = PetscObjectTypeCompare((PetscObject)x,VECSEQCUSP,&isseqcusp);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)x,VECMPICUSP,&ismpicusp);CHKERRQ(ierr);
  iscusp = (isseqcusp || ismpicusp) ? PETSC_TRUE : PETSC_FALSE;
  ierr = PetscObjectTypeCompare((PetscObject)x,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);

  if (isseqgpu) {
    dim3 block(BLOCKWIDTH_X, BLOCKWIDTH_Y);
    dim3 grid((int)ceil((float)((mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/mat->size)/(float)BLOCKWIDTH_X / 1.0), 1);

//     ierr = VecScatterBegin(mat->Mvctx,x,mat->lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
//    ierr = VecScatterEnd(mat->Mvctx,x,mat->lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

    int shared_size = 0;
    Vec_SeqGPU *vx = (Vec_SeqGPU*) x->data;
    Vec_SeqGPU *vy = (Vec_SeqGPU*) y->data;
    /* Make sure y is also VECSEQGPU */
    ierr = PetscObjectTypeCompare((PetscObject)x,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
    if (!isseqgpu) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and y must be same type");
    }


    /* synch up x */
    if (vx->syncState==VEC_CPU) {
      ierr = VecCopyOverH2D(x,vx->cpuptr);CHKERRQ(ierr);
      vx->syncState=VEC_SYNCHED;
    }
    /* Get device pointer for X */
    devX = vx->devptr;
    devY = vy->devptr;
 

   /* Bind X to device texture */
    mat_size = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/mat->size;
    vec_size = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof);

    checkCudaError(cudaBindTexture(0, vector_x, devX, vec_size * sizeof(PetscScalar)));

    MatMultKernelMPI<<<grid, block, shared_size, mat_seq->stream>>>(mat_seq->deviceData, devY, devX, mat_size, mat_seq->diagonals->size(), mat_seq->deviceDiags, mat_seq->dof, vec_size);

    cudaUnbindTexture(vector_x);
    cudaDeviceSynchronize();

  } else if (iscusp) {
    dim3 block(BLOCKWIDTH_X, BLOCKWIDTH_Y);
    dim3 grid((int)ceil((float)((mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/mat->size)/(float)BLOCKWIDTH_X / 1.0), 1);
 
    int shared_size = 0;
    /* Make sure y is also VECCUSP */
    ierr = PetscObjectTypeCompare((PetscObject)x,VECCUSP,&isseqgpu);CHKERRQ(ierr);
    if (!iscusp) 
    {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and y must be same type");
    }
 
    ierr = VecScatterBegin(mat->Mvctx,x,mat->lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(mat->Mvctx,x,mat->lvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
 
    mat_size = ((mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/mat->size);
    vec_size = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof);
 
    ierr = VecCUSPGetArrayWrite(y, &ygpu);CHKERRQ(ierr);
    ierr = VecCUSPGetArrayRead(mat->lvec, &xgpu);CHKERRQ(ierr);
    devY = thrust::raw_pointer_cast(&(*ygpu)[0]);
    devX = thrust::raw_pointer_cast(&(*xgpu)[0]);
 
    /* Bind X to device texture */
    checkCudaError(cudaBindTexture(0, vector_x, devX, vec_size * sizeof(PetscScalar)));
 
#if _TRACE
    printf("Host diagonals:\n");
    for (int i = 0; i < mat_seq->diagonals->size(); ++i) {
      printf("- %d\n", (*mat_seq->diagonals)[i]);
    }
#endif

    /* Invoke */

#if _TIME
    double start, end;
    start = getclock();
#endif

cudaPrintfInit();
MatMultKernelMPI<<<grid, block, shared_size, mat_seq->stream>>>(mat_seq->deviceData, devY, devX, mat_size, mat_seq->diagonals->size(), mat_seq->deviceDiags, mat_seq->dof,vec_size);
cudaPrintfDisplay(stdout,true);
cudaPrintfEnd();

#if _TIME
    checkCudaError(cudaStreamSynchronize(mat_seq->stream));
    end = getclock();
    double elapsed = end - start;
    double gflops = (2.0 * mat_seq->non_zeros / elapsed / 1e9);

    double nos = ((mat_seq->p == 1 ? 2 : 3) * 2 + 1) * (2*mat_seq->dof - 1);
    double nz = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/numprocs;
    double alt_gflops = (2.0 * nos * nz) / ((end - start)*1024*1024*1024);

#if _CSV_OUT
    fprintf(stderr, "%d,%lf,%lf,\n", (mat_seq->m, mat_seq->n, mat_seq->p, mat_seq->dof)/numprocs, elapsed, gflops);
#endif
    printf("SGGPU Kernel Time:           %lf sec\n", elapsed);
    printf("SGGPU Kernel GFlop/s:        %lf\n", gflops);
    printf("SGGPU Kernel GFlop/s (alt):  %lf\n", alt_gflops);
#endif

    /* Cleanup */
    cudaUnbindTexture(vector_x);

    ierr = VecCUSPRestoreArrayRead(x, &xgpu); CHKERRQ(ierr);
    ierr = VecCUSPRestoreArrayWrite(y, &ygpu); CHKERRQ(ierr);
    ierr = WaitForGPU() ; CHKERRCUSP(ierr);
    cudaDeviceSynchronize();
  } else {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Vec must be seqgpu or cusp type");
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_MPISGGPU"
PetscErrorCode MatSetValuesBlocked_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is) {
  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatSetValuesBlocked_MPISGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValues_MPISGGPU"
PetscErrorCode MatSetValues_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is) {

  int i, j;
  PetscErrorCode ierr;
  PetscBool resizegpu = PETSC_FALSE;
  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;

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
      
	      int diff = col - row;
	      int left = row % mat_seq->dof;
	      int diag = int(floor((double)(diff + left) / mat_seq->dof));
             
              if(mat->rank > 0)
		diag = rstart + diag;

	      int col_offset = col % mat_seq->dof;
	      int num_elems = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/mat->size;
	      int offset = col_offset * num_elems + row - (rank*num_elems);

#if _TRACE
      printf("- row: %d  col: %d  val: %lf  diag: %d  offset: %d\n", row, col, y[i*ncol+j], diag, offset);
#endif

      std::map<int, int> &diag_starts = *(mat_seq->diag_starts);
      std::map<int, int>::iterator I = diag_starts.find(diag);
      int diag_offset = 0;
      if (I == diag_starts.end()) {
        printf("WARNING: malloc() in MatSetValues\n");
        resizegpu = PETSC_TRUE;
        // The diagonal does not yet exist, so add a new diagonal
        int num_diags = diag_starts.size() + 1;
        int size = num_diags * ((mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof * mat_seq->dof)/mat->size);
        PetscScalar *newData;
        ierr = PetscMalloc(size * sizeof(PetscScalar), &newData); CHKERRQ(ierr);
        memset(newData, 0, size * sizeof(PetscScalar));
        size -= ((mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof * mat_seq->dof)/mat->size);
        if (num_diags > 1) {
          // This is not the first diagonal, so copy
#if _TRACE
          printf("- Memcpy of %d elements\n", size);
#endif
          memcpy(newData, mat_seq->hostData, size * sizeof(PetscScalar));
        }
        PetscFree(mat_seq->hostData);
        mat_seq->hostData = newData;
        diag_offset = size;
        diag_starts[diag] = diag_offset;
        mat_seq->diagonals->push_back(diag);
      }
	
	else 
      {
        // The diagonal already exists, so get the base offset
        diag_offset = I->second;
      }

      diag_offset += offset;

      if (is == INSERT_VALUES)
        mat_seq->hostData[diag_offset] = y[i * ncol + j];
      else
        mat_seq->hostData[diag_offset] += y[i * ncol + j];

      mat_seq->non_zeros++;
	}
    }
  }

  if (resizegpu) {
    int size,mat_size;
    // Create GPU buffer
    if (mat_seq->deviceData) {
      cudaFree(mat_seq->deviceData);
    }
    size = mat_seq->diag_starts->size() * ((mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof * mat_seq->dof)/numprocs);
    checkCudaError(cudaMalloc(&mat_seq->deviceData, sizeof(PetscScalar) * size));

    mat_size = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/numprocs;

    if (mat_seq->deviceX) {
      cudaFree(mat_seq->deviceX);
    }
    if (mat_seq->deviceY) {
      cudaFree(mat_seq->deviceY);
    }
    if (mat_seq->deviceDiags) {
      cudaFree(mat_seq->deviceDiags);
    }
    // We know the expected size of x, y, so go ahead and allocate them now
    checkCudaError(cudaMalloc(&mat_seq->deviceX, mat_size * sizeof(PetscScalar)));
    checkCudaError(cudaMalloc(&mat_seq->deviceY, mat_size * sizeof(PetscScalar)));

    // We also know how many diagonals we have, and their indices
    checkCudaError(cudaMalloc(&mat_seq->deviceDiags, sizeof(int) * mat_seq->diagonals->size()));
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
PetscErrorCode MatSetUp_MPISGGPU(Mat A) {

  PetscFunctionBegin;
  SGTrace;

  PetscPrintf(PETSC_COMM_WORLD,"MatSetUp_MPISGGPU() not implemented\n");
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_MPISGGPU"
PetscErrorCode MatZeroEntries_MPISGGPU(Mat A) {

  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;
  PetscInt size;
  PetscFunctionBegin;
  SGTrace;
  
  size = mat_seq->diag_starts->size() * ((mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof * mat_seq->dof)/mat->size);
  memset(mat_seq->hostData, 0, size * sizeof(PetscScalar));
  
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_MPISGGPU"
PetscErrorCode MatGetDiagonal_MPISGGPU(Mat A, Vec v) {

  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatGetDiagonal_MPISGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_MPISGGPU"
PetscErrorCode MatDiagonalScale_MPISGGPU(Mat A, Vec ll, Vec rr) {

  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatDiagonalScale_MPISGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRow_MPISGGPU"
PetscErrorCode MatGetRow_MPISGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v) {

  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatGetRow_MPISGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_MPISGGPU"
PetscErrorCode MatRestoreRow_MPISGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v) {

  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatRestoreRow_MPISGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_MPISGGPU"
PetscErrorCode MatGetRowMaxAbs_MPISGGPU(Mat A, Vec v, PetscInt idx[]) {

  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatGetRowMaxAbs_MPISGGPU not implemented");
}


void DisplayLocalMatrix(Mat_SeqSGGPU *mat_seq);

#undef __FUNCT__
#define __FUNCT__ "MatView_MPISGGPU"
PetscErrorCode MatView_MPISGGPU(Mat A, PetscViewer viewer) 
{
  PetscInt i;
  Mat_MPISGGPU* mat = (Mat_MPISGGPU*)A->data;
  
  for(i = 0; i < mat->size;++i)
	{
	 if(i == mat->rank)
	    DisplayLocalMatrix(mat->mat_seq);
         MPI_Barrier(PETSC_COMM_WORLD);   
    	}
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyBegin_MPISGGPU"
PetscErrorCode MatAssemblyBegin_MPISGGPU(Mat A, MatAssemblyType type) {

  PetscFunctionBegin;
  SGTrace;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyEnd_MPISGGPU"
PetscErrorCode MatAssemblyEnd_MPISGGPU(Mat A, MatAssemblyType type) {

  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;

  PetscInt size;

  PetscFunctionBegin;

#if _TRACE
  printf("[SeqSGGPU] MatAssemblyEnd_SeqSGGPU\n");

  for (std::map<int, int>::iterator I = mat_seq->diag_starts->begin(),
       E = mat_seq->diag_starts->end(); I != E; ++I) {
    printf("- Diag %d:\n", I->first);
    for (int i = 0; i < mat_seq->dof; ++i) {
      for (int j = 0; j < (mat_seq->dof * mat_seq->m * mat_seq->n * mat_seq->p)/numprocs; ++j) {
        int offset = i * ((mat_seq->dof * mat_seq->m * mat_seq->n * mat_seq->p)/numprocs) + j;
        printf(" %lf ", mat_seq->hostData[offset + I->second]);
      }
      printf("\n");
    }
  }
#endif

  size = (mat_seq->diag_starts->size()*mat_seq->m*mat_seq->n*mat_seq->p*mat_seq->dof*mat_seq->dof)/mat->size;

  checkCudaError(cudaMemcpyAsync(mat_seq->deviceDiags, &(*mat_seq->diagonals)[0], sizeof(int) * mat_seq->diagonals->size(), cudaMemcpyHostToDevice, mat_seq->stream));

  checkCudaError(cudaMemcpy(mat_seq->deviceData, mat_seq->hostData, sizeof(PetscScalar) * size, cudaMemcpyHostToDevice));

  cudaDeviceSynchronize();
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "MatFDColoringApply_MPISGGPU"
PetscErrorCode  MatFDColoringApply_MPISGGPU(Mat A,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx) {

  MatFDColoringApply_SeqSGGPU(A, coloring, x1, flag, sctx);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatFDColoringCreate_MPISGGPU"
PetscErrorCode MatFDColoringCreate_MPISGGPU(Mat A,ISColoring iscoloring,MatFDColoring c) {
 
//  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
//  Mat_SeqSGGPU *mat_seq = mat->mat_seq;

  PetscErrorCode ierr;
  PetscInt       i,n,nrows,N,j,k,m,ncols,col;
  const PetscInt *is,*ci,*cj,*rows;
  PetscInt       nis = iscoloring->n,*rowhit,*columnsforrow,l,bs = 1;
  IS             *isa;
  PetscBool      done,flg = PETSC_FALSE;

  PetscFunctionBegin;


	PetscPrintf(PETSC_COMM_WORLD,"MatFDColoringCreate_SeqSGGPU\n");


  ierr = ISColoringGetIS(iscoloring,PETSC_IGNORE,&isa);CHKERRQ(ierr);
  /* this is ugly way to get blocksize but cannot call MatGetBlockSize() because AIJ can have bs > 1 */

  N          = A->cmap->N/bs;
  c->M       = A->rmap->N/bs;  /* set total rows, columns and local rows */
  c->N       = A->cmap->N/bs;
  c->m       = A->rmap->N/bs;
  c->rstart  = 0;

  c->ncolors = nis;
  ierr       = PetscMalloc(nis*sizeof(PetscInt),&c->ncolumns);CHKERRQ(ierr);
  ierr       = PetscMalloc(nis*sizeof(PetscInt*),&c->columns);CHKERRQ(ierr); 
  ierr       = PetscMalloc(nis*sizeof(PetscInt),&c->nrows);CHKERRQ(ierr);
  ierr       = PetscMalloc(nis*sizeof(PetscInt*),&c->rows);CHKERRQ(ierr);
  ierr       = PetscMalloc(nis*sizeof(PetscInt*),&c->columnsforrow);CHKERRQ(ierr);

  ierr = MatGetColumnIJ(A,0,PETSC_FALSE,PETSC_FALSE,&ncols,&ci,&cj,&done);CHKERRQ(ierr);
  if (!done) SETERRQ1(((PetscObject)A)->comm,PETSC_ERR_SUP,"MatGetColumnIJ() not supported for matrix type %s",((PetscObject)A)->type_name);

  /*
     Temporary option to allow for debugging/testing
  */
  ierr = PetscOptionsGetBool(PETSC_NULL,"-matfdcoloring_slow",&flg,PETSC_NULL);CHKERRQ(ierr);

  ierr = PetscMalloc((N+1)*sizeof(PetscInt),&rowhit);CHKERRQ(ierr);
  ierr = PetscMalloc((N+1)*sizeof(PetscInt),&columnsforrow);CHKERRQ(ierr);

  for (i=0; i<nis; i++) {
    ierr = ISGetLocalSize(isa[i],&n);CHKERRQ(ierr);
    ierr = ISGetIndices(isa[i],&is);CHKERRQ(ierr);
    c->ncolumns[i] = n;
    if (n) {
      ierr = PetscMalloc(n*sizeof(PetscInt),&c->columns[i]);CHKERRQ(ierr);
      ierr = PetscMemcpy(c->columns[i],is,n*sizeof(PetscInt));CHKERRQ(ierr);
    } else {
      c->columns[i]  = 0;
    }

    if (!flg) { /* ------------------------------------------------------------------------------*/
      /* fast, crude version requires O(N*N) work */
      ierr = PetscMemzero(rowhit,N*sizeof(PetscInt));CHKERRQ(ierr);
      /* loop over columns*/
      for (j=0; j<n; j++) {
        col  = is[j];
        rows = cj + ci[col]; 
        m    = ci[col+1] - ci[col];
        /* loop over columns marking them in rowhit */
        for (k=0; k<m; k++) {
          rowhit[*rows++] = col + 1;
        }
      }
      /* count the number of hits */
      nrows = 0;
      for (j=0; j<N; j++) {
        if (rowhit[j]) nrows++;
      }
      c->nrows[i] = nrows;
      ierr        = PetscMalloc((nrows+1)*sizeof(PetscInt),&c->rows[i]);CHKERRQ(ierr);
      ierr        = PetscMalloc((nrows+1)*sizeof(PetscInt),&c->columnsforrow[i]);CHKERRQ(ierr);
      nrows       = 0;
      for (j=0; j<N; j++) {
        if (rowhit[j]) {
          c->rows[i][nrows]          = j;
          c->columnsforrow[i][nrows] = rowhit[j] - 1;
          nrows++;
        }
      }
    } else {  /*-------------------------------------------------------------------------------*/
      /* slow version, using rowhit as a linked list */
      PetscInt currentcol,fm,mfm;
      rowhit[N] = N;
      nrows     = 0;
      /* loop over columns */
      for (j=0; j<n; j++) {
        col   = is[j];
        rows  = cj + ci[col]; 
        m     = ci[col+1] - ci[col];
        /* loop over columns marking them in rowhit */
        fm    = N; /* fm points to first entry in linked list */
        for (k=0; k<m; k++) {
          currentcol = *rows++;
	  /* is it already in the list? */
          do {
            mfm  = fm;
            fm   = rowhit[fm];
          } while (fm < currentcol);
          /* not in list so add it */
          if (fm != currentcol) {
            nrows++;
            columnsforrow[currentcol] = col;
            /* next three lines insert new entry into linked list */
            rowhit[mfm]               = currentcol;
            rowhit[currentcol]        = fm;
            fm                        = currentcol; 
            /* fm points to present position in list since we know the columns are sorted */
          } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Detected invalid coloring");
        }
      }
      c->nrows[i] = nrows;
      ierr        = PetscMalloc((nrows+1)*sizeof(PetscInt),&c->rows[i]);CHKERRQ(ierr);
      ierr        = PetscMalloc((nrows+1)*sizeof(PetscInt),&c->columnsforrow[i]);CHKERRQ(ierr);
      /* now store the linked list of rows into c->rows[i] */
      nrows       = 0;
      fm          = rowhit[N];
      do {
        c->rows[i][nrows]            = fm;
        c->columnsforrow[i][nrows++] = columnsforrow[fm];
        fm                           = rowhit[fm];
      } while (fm < N);
    } /* ---------------------------------------------------------------------------------------*/
    ierr = ISRestoreIndices(isa[i],&is);CHKERRQ(ierr);  
  }
  ierr = MatRestoreColumnIJ(A,0,PETSC_FALSE,PETSC_FALSE,&ncols,&ci,&cj,&done);CHKERRQ(ierr);

  ierr = PetscFree(rowhit);CHKERRQ(ierr);
  ierr = PetscFree(columnsforrow);CHKERRQ(ierr);

  /* Optimize by adding the vscale, and scaleforrow[][] fields */
  /*
       see the version for MPIAIJ
  */
  ierr = VecCreateGhost(((PetscObject)A)->comm,A->rmap->n,PETSC_DETERMINE,0,PETSC_NULL,&c->vscale);CHKERRQ(ierr);
  ierr = PetscMalloc(c->ncolors*sizeof(PetscInt*),&c->vscaleforrow);CHKERRQ(ierr);
  for (k=0; k<c->ncolors; k++) { 
    ierr = PetscMalloc((c->nrows[k]+1)*sizeof(PetscInt),&c->vscaleforrow[k]);CHKERRQ(ierr);
    for (l=0; l<c->nrows[k]; l++) {
      col = c->columnsforrow[k][l];
      c->vscaleforrow[k][l] = col;
    }
  }
  ierr = ISColoringRestoreIS(iscoloring,&isa);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetColumnIJ_MPISGGPU"
PetscErrorCode MatGetColumnIJ_MPISGGPU(Mat A,PetscInt oshift,PetscBool  symmetric,PetscBool  inodecompressed,PetscInt *nn, const PetscInt *ia[], const PetscInt *ja[],PetscBool  *done) {

  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *a = mat->mat_seq;  

  PetscErrorCode ierr;
  PetscInt       n = A->cmap->n;
  PetscInt       ndiag = a->diagonals->size();
  PetscInt       nrows = a->m*a->n*a->p*a->dof;
  PetscInt       nz=a->dof*ndiag*nrows;
  PetscInt       iblock,i,j,col,index,colblock,offset;

  PetscFunctionBegin;  

  *nn = nrows;

  if (!ia) PetscFunctionReturn(0);
  if (a->ja) {
    ierr = PetscFree(a->ja); CHKERRQ(ierr);
  }
  if (a->ia) {
    ierr = PetscFree(a->ia); CHKERRQ(ierr);
  }
  ierr = PetscMalloc((n+1)*sizeof(PetscInt),&a->ia);CHKERRQ(ierr);
  ierr = PetscMalloc((nz+1)*sizeof(PetscInt),&a->ja);CHKERRQ(ierr);

  /* Assuming symmetric nonzero structure */
  index=0;
  for (iblock=0;iblock<nrows/a->dof;iblock++) {
    for (i=iblock*a->dof;i<(iblock+1)*a->dof;i++) {
      a->ia[i] = index;
      for (std::map<int, int>::iterator I = a->diag_starts->begin(),
             E = a->diag_starts->end(); I != E; ++I) {
        offset = I->first;
        colblock = offset + iblock;
        /* Ignore 0 padding */
        if (colblock < 0) {
          continue;
        }
        if (colblock >= (nrows/a->dof)) {
          break;
        }
        /* skip some blocks for nonperiodic da */
        if (a->stencil_type==0 &&  a->dim==2 && 
            ((colblock - iblock == 1 && !(colblock % a->n)) ||
             (iblock - colblock == 1 && !(iblock % a->n)))) {
          continue;

        }
      
        for (j=0;j<a->dof;j++) {
          col = (colblock*a->dof)  + j;
          a->ja[index++] = col;
	}

      }
    }
  }
  a->ia[nrows] = index;
  *ia = a->ia;
  *ja = a->ja;

  PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatMPISGGPUSetPreallocation"
PetscErrorCode MatMPISGGPUSetPreallocation(Mat A,PetscInt stencil_type, PetscInt dof)
{
  PetscErrorCode ierr;
  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;  

  PetscFunctionBegin;

  mat_seq->stencil_type = stencil_type;
  mat_seq->dof = dof;
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
  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;

  PetscInt dim,diag_size,size,num_diags,i,vecsize;

  ierr = PetscLayoutSetBlockSize(A->rmap,1);CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(A->cmap,1);CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(A->rmap);CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(A->cmap);CHKERRQ(ierr);

  dim = A->stencil.dim;
  if (mat_seq->dof > 1) {
    dim--;
  }

  PetscInt rstart = A->rmap->rstart;
  mat_seq->dim = dim;

  if (mat_seq->stencil_type == 0) {
    /* star stencil */
    num_diags = 2*mat_seq->dim + 1;
  } else {
    /* box stencil */
    num_diags =  1;
    for (i=0;i<mat_seq->dim;i++) num_diags*=3;
  }

  diag_size = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof * mat_seq->dof)/mat->size;
  size = num_diags * diag_size;

  if (mat_seq->m == 0 || mat_seq->n == 0 || mat_seq->p == 0 || mat_seq->dof == 0) {
    SETERRQ(PETSC_COMM_SELF,0,"MatSetPreallocation_SeqSGGPU called without valid m, n, p, and dof!");
  }

  ierr = PetscMalloc(sizeof(PetscInt)*num_diags,&mat_seq->diag_offsets);
  ierr = PetscMalloc(size * sizeof(PetscScalar), &mat_seq->hostData); CHKERRQ(ierr);
  memset(mat_seq->hostData, 0, size * sizeof(PetscScalar));

  (*mat_seq->diag_starts)[rstart + 0]  = 0 * diag_size;
  (*mat_seq->diagonals).push_back(rstart + 0);
  (*mat_seq->diag_starts)[rstart + 1]  = 1 * diag_size;
  (*mat_seq->diagonals).push_back(rstart + 1);
  (*mat_seq->diag_starts)[rstart - 1] = 2 * diag_size;
  (*mat_seq->diagonals).push_back(rstart - 1);
  if (mat_seq->stencil_type == 0) {
    if (mat_seq->dim == 2) 
    {
      (*mat_seq->diag_starts)[rstart + mat_seq->m] = 3 * diag_size;
      (*mat_seq->diagonals).push_back(rstart + mat_seq->m);
      (*mat_seq->diag_starts)[rstart - mat_seq->m] = 4 * diag_size;
      (*mat_seq->diagonals).push_back(rstart - mat_seq->m);

    } 
    else if (mat_seq->dim == 3) 
    {
      (*mat_seq->diag_starts)[rstart + mat_seq->m] = 3 * diag_size;
      (*mat_seq->diagonals).push_back(rstart + mat_seq->m);
      (*mat_seq->diag_starts)[rstart - mat_seq->m] = 4 * diag_size;
      (*mat_seq->diagonals).push_back(rstart - mat_seq->m);

      (*mat_seq->diag_starts)[rstart + mat_seq->m*mat_seq->n] = 5 * diag_size;
      (*mat_seq->diagonals).push_back(rstart + mat_seq->m*mat_seq->n);
      (*mat_seq->diag_starts)[rstart - mat_seq->m*mat_seq->n] = 6 * diag_size;
      (*mat_seq->diagonals).push_back(rstart - mat_seq->m*mat_seq->n);
    }
  } else {
    if (mat_seq->dim == 2) {
      (*mat_seq->diag_starts)[rstart + mat_seq->n-1] = 3 * diag_size;
      (*mat_seq->diagonals).push_back(rstart + mat_seq->m);
      (*mat_seq->diag_starts)[rstart - mat_seq->n-1] = 4 * diag_size;
      (*mat_seq->diagonals).push_back(rstart - mat_seq->m);
      (*mat_seq->diag_starts)[rstart + mat_seq->n] = 5 * diag_size;
      (*mat_seq->diagonals).push_back(rstart + mat_seq->m);
      (*mat_seq->diag_starts)[rstart - mat_seq->n] = 6 * diag_size;
      (*mat_seq->diagonals).push_back(rstart - mat_seq->m);
      (*mat_seq->diag_starts)[rstart + mat_seq->n+1] = 7 * diag_size;
      (*mat_seq->diagonals).push_back(rstart + mat_seq->m);
      (*mat_seq->diag_starts)[rstart - mat_seq->n+1] = 8 * diag_size;
      (*mat_seq->diagonals).push_back(rstart - mat_seq->m);
    }
  }
  /*
  printf("Diagonals preallocated:\n");
  for (std::map<int, int>::iterator I = mat_seq->diag_starts->begin(),
         E = mat_seq->diag_starts->end(); I != E; ++I) {
    printf("%4d --> %4d\n",I->first,I->second);
  }
   */
  
  
  // Create GPU buffer
  if (mat_seq->deviceData) {
    cudaFree(mat_seq->deviceData);
  }
  checkCudaError(cudaMalloc(&mat_seq->deviceData, sizeof(PetscScalar) * size));
  checkCudaError(cudaMemset(mat_seq->deviceData,0,sizeof(PetscScalar)*size));

  // Copy data to device
  checkCudaError(cudaMemcpy(mat_seq->deviceData, mat_seq->hostData, sizeof(PetscScalar) * size, cudaMemcpyHostToDevice));

  vecsize = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/mat->size;

  // We know the expected size of x, y, so go ahead and allocate them now
  checkCudaError(cudaMalloc(&mat_seq->deviceX, vecsize * sizeof(PetscScalar)));
  checkCudaError(cudaMalloc(&mat_seq->deviceY, vecsize * sizeof(PetscScalar)));

  // We also know how many diagonals we have, and their indices
  checkCudaError(cudaMalloc(&mat_seq->deviceDiags, sizeof(int) * mat_seq->diagonals->size()));
  A->preallocated = PETSC_TRUE;
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  MatSetUpMultiply_MPISGGPU(A);

  PetscFunctionReturn(0);
}
EXTERN_C_END



void DisplayLocalMatrix(Mat_SeqSGGPU *mat_seq)
{

  PetscInt nrows,ndiag,dof,i,j,iblock,col,index,offset;
  std::map<int, int> &diag_starts = *(mat_seq->diag_starts);

  PetscInt numprocs, rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);

  nrows = (mat_seq->m * mat_seq->n * mat_seq->p * mat_seq->dof)/numprocs;
  ndiag = mat_seq->diagonals->size();
  dof = mat_seq->dof;

  for (std::map<int, int>::iterator I = diag_starts.begin(),
         E = diag_starts.end(); I != E; ++I) 
	{
		fprintf(stdout,"- Diag %d:%d\n", I->first, I->second);
  	}

	fprintf(stdout,"\n");
	fprintf(stdout,"hostData:\n");

  for (i=0;i<nrows;i++) 
	{
	fprintf(stdout,"row %2.2d:",rank*nrows + i); 

	    for (j=0;j<ndiag*dof;j++) 
		{
			fprintf(stdout," %4g ",mat_seq->hostData[i+j*nrows]);
		}
	fprintf(stdout,"\n");
  	}
	fprintf(stdout,"\n\n");


  for (iblock=0;iblock<(nrows/dof);iblock++)  
	{
	    for (i=iblock*dof;i<(iblock+1)*dof;i++) 
		{
		fprintf(stdout,"row %d:",rank*nrows + i);
	
		  for (std::map<int, int>::iterator I = mat_seq->diag_starts->begin(),
        	     E = mat_seq->diag_starts->end(); I != E; ++I) 
		  {
	        	/* Ignore 0 padding */
		        offset = I->first;

		        if (offset + iblock + (rank*nrows) < 0) 
			{
	        	  continue;
        		}

	        	if (offset + iblock + (rank*nrows) >= ((nrows*numprocs)/dof)) 
			{
        	  	break;
        		}
	        
	        	for (j=0;j<dof;j++) 
			{
	        	  index = i + I->second + j*nrows; // column-major
	        	  col = offset*dof+((iblock+(rank*nrows))*dof) + j;
	                  fprintf(stdout," (%d, %g) ",col,mat_seq->hostData[index]);
        		}
      	 	  }
		fprintf(stdout,"\n");
    		}
  	}

}
