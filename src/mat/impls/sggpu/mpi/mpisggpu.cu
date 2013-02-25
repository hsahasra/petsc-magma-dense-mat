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
     
__global__ void MatMultKernelMPI(PetscScalar * coeff, PetscScalar * y, PetscScalar *x,PetscInt mat_size, PetscInt num_diags, int * diagonals, PetscInt dof) {
       
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
	      ///PetscScalar xval0 = fetch_doubleMPI(vector_x, block0 + j);
	      PetscScalar xval0 = x[block0 + j];

	      yval0 += aval0 * xval0;

	      //21, 22, 25, 26, 37, 38, 41, 42 are the only non zero entries in the 64x1 vector for ex14 when grid size is 4x4x4 	
	      //The following if statement is for testing and debugging and can be removed. 
	      //if ((idx0 == 21) || (idx0 == 22) || (idx0 == 25) || (idx0 == 26) ||  (idx0 == 37) ||  (idx0 == 38) ||  (idx0 == 41) ||  (idx0 == 42))
	      //  cuPrintf("d:%d \t offset:%d \t block0:%d \t aval0:%lf \t xval0:%lf \t yval0 :%lf\n",d,offset0,block0,aval0,xval0, yval0);
        }
}
   
      y[idx0] = yval0;

//      cuPrintf("y[%d]:%g\n",idx0,y[idx0]);

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
  Mat_SeqSGGPU * mat;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  SGTrace;

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


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

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  mat->m = m;
  mat->n = n;
  mat->p = p;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatMult_MPISGGPU"
PetscErrorCode MatMult_MPISGGPU(Mat A, Vec x, Vec y) {

  Mat_SeqSGGPU *mat = (Mat_SeqSGGPU*)A->data;

//  MatView_MPISGGPU(A,PETSC_VIEWER_STDOUT_WORLD);	

  PetscInt rank;
  PetscInt numprocs;	
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);	

  PetscBool isseqcusp,isseqgpu,ismpicusp,iscusp;
  PetscErrorCode ierr;
  PetscInt mat_size;
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
    dim3 grid((int)ceil((float)((mat->m * mat->n * mat->p * mat->dof)/numprocs)/(float)BLOCKWIDTH_X / 1.0), 1);

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
    mat_size = (mat->m * mat->n * mat->p * mat->dof)/numprocs;

    checkCudaError(cudaBindTexture(0, vector_x, devX, mat_size * sizeof(PetscScalar)));

    MatMultKernelMPI<<<grid, block, shared_size, mat->stream>>>(mat->deviceData, devY, devX, mat_size, mat->diagonals->size(), mat->deviceDiags, mat->dof);

    cudaUnbindTexture(vector_x);
    cudaDeviceSynchronize();


  } else if (iscusp) {
    dim3 block(BLOCKWIDTH_X, BLOCKWIDTH_Y);
    dim3 grid((int)ceil((float)((mat->m * mat->n * mat->p * mat->dof)/numprocs)/(float)BLOCKWIDTH_X / 1.0), 1);

    int shared_size = 0;
    /* Make sure y is also VECCUSP */
    ierr = PetscObjectTypeCompare((PetscObject)x,VECCUSP,&isseqgpu);CHKERRQ(ierr);
    if (!iscusp) 
    {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and y must be same type");
    }
	
    mat_size = ((mat->m * mat->n * mat->p * mat->dof)/numprocs);

   // Vec xx;
   // ierr = VecCreateSeq(PETSC_COMM_SELF,(mat->m*mat->n*mat->p*mat->dof),&xx);
   // VecCopy(x,xx);	 
   // VecView(x,PETSC_VIEWER_STDOUT_WORLD);
   // VecView(xx,PETSC_VIEWER_STDOUT_WORLD);

    ierr = VecCUSPGetArrayWrite(y, &ygpu); CHKERRQ(ierr);
    ierr = VecCUSPGetArrayRead(x, &xgpu); CHKERRQ(ierr);
    devY = thrust::raw_pointer_cast(&(*ygpu)[0]);
    devX = thrust::raw_pointer_cast(&(*xgpu)[0]);

    /* Bind X to device texture */
    checkCudaError(cudaBindTexture(0, vector_x, devX, mat_size * sizeof(PetscScalar)));

#if _TRACE
    printf("Host diagonals:\n");
    for (int i = 0; i < mat->diagonals->size(); ++i) {
      printf("- %d\n", (*mat->diagonals)[i]);
    }
#endif

    /* Invoke */

#if _TIME
    double start, end;
    start = getclock();
#endif

cudaPrintfInit();
MatMultKernelMPI<<<grid, block, shared_size, mat->stream>>>(mat->deviceData, devY, devX, mat_size, mat->diagonals->size(), mat->deviceDiags, mat->dof);
cudaPrintfDisplay(stdout,true);
cudaPrintfEnd();

#if _TIME
    checkCudaError(cudaStreamSynchronize(mat->stream));
    end = getclock();
    double elapsed = end - start;
    double gflops = (2.0 * mat->non_zeros / elapsed / 1e9);

    double nos = ((mat->p == 1 ? 2 : 3) * 2 + 1) * (2*mat->dof - 1);
    double nz = (mat->m * mat->n * mat->p * mat->dof)/numprocs;
    double alt_gflops = (2.0 * nos * nz) / ((end - start)*1024*1024*1024);

#if _CSV_OUT
    fprintf(stderr, "%d,%d,%d,%d,%lf,%lf,\n", (mat->m, mat->n, mat->p, mat->dof)/numprocs, elapsed, gflops);
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

//	VecView(x,PETSC_VIEWER_STDOUT_WORLD);
//	VecView(y,PETSC_VIEWER_STDOUT_WORLD);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_MPISGGPU"
PetscErrorCode MatSetValuesBlocked_MPISGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is) {
  MatSetValuesBlocked_SeqSGGPU(A, nrow, irow, ncol, icol, y, is);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
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
	      int num_elems = (mat->m * mat->n * mat->p * mat->dof)/numprocs;
	      int offset = col_offset * num_elems + row - (rank*num_elems);

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
        int size = num_diags * ((mat->m * mat->n * mat->p * mat->dof * mat->dof)/numprocs);
        PetscScalar *newData;
        ierr = PetscMalloc(size * sizeof(PetscScalar), &newData); CHKERRQ(ierr);
        memset(newData, 0, size * sizeof(PetscScalar));
        size -= ((mat->m * mat->n * mat->p * mat->dof * mat->dof)/numprocs);
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
      }
	
	else 
      {
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
    size = mat->diag_starts->size() * ((mat->m * mat->n * mat->p * mat->dof * mat->dof)/numprocs);
    checkCudaError(cudaMalloc(&mat->deviceData, sizeof(PetscScalar) * size));

    mat_size = (mat->m * mat->n * mat->p * mat->dof)/numprocs;

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

	PetscInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  MatSetUp_SeqSGGPU(mat);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_MPISGGPU"
PetscErrorCode MatZeroEntries_MPISGGPU(Mat A) {

  Mat_SeqSGGPU *mat = (Mat_SeqSGGPU*)A->data;
  PetscInt size;
  PetscInt rank, numprocs;
  PetscFunctionBegin;
  SGTrace;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);

  size = mat->diag_starts->size() * ((mat->m * mat->n * mat->p * mat->dof * mat->dof)/numprocs);
  memset(mat->hostData, 0, size * sizeof(PetscScalar));
  
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_MPISGGPU"
PetscErrorCode MatGetDiagonal_MPISGGPU(Mat A, Vec v) {

	PetscInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  MatGetDiagonal_SeqSGGPU(A, v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_MPISGGPU"
PetscErrorCode MatDiagonalScale_MPISGGPU(Mat A, Vec ll, Vec rr) {

	PetscInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  MatDiagonalScale_SeqSGGPU(A, ll, rr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRow_MPISGGPU"
PetscErrorCode MatGetRow_MPISGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v) {

	PetscInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  MatGetRow_SeqSGGPU(A, row, nz, idx , v);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_MPISGGPU"
PetscErrorCode MatRestoreRow_MPISGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v) {
  MatRestoreRow_SeqSGGPU(A, row, nz, idx, v);

	PetscInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_MPISGGPU"
PetscErrorCode MatGetRowMaxAbs_MPISGGPU(Mat A, Vec v, PetscInt idx[]) {

	PetscInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;
  PetscInt size;
  PetscInt numprocs;
  PetscFunctionBegin;

  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);	

#if _TRACE
  printf("[SeqSGGPU] MatAssemblyEnd_SeqSGGPU\n");

  for (std::map<int, int>::iterator I = mat->diag_starts->begin(),
       E = mat->diag_starts->end(); I != E; ++I) {
    printf("- Diag %d:\n", I->first);
    for (int i = 0; i < mat->dof; ++i) {
      for (int j = 0; j < (mat->dof * mat->m * mat->n * mat->p)/numprocs; ++j) {
        int offset = i * ((mat->dof * mat->m * mat->n * mat->p)/numprocs) + j;
        printf(" %lf ", mat->hostData[offset + I->second]);
      }
      printf("\n");
    }
  }
#endif

  size = (mat->diag_starts->size()*mat->m*mat->n*mat->p*mat->dof*mat->dof)/numprocs;

  checkCudaError(cudaMemcpyAsync(mat->deviceDiags, &(*mat->diagonals)[0], sizeof(int) * mat->diagonals->size(), cudaMemcpyHostToDevice, mat->stream));

  checkCudaError(cudaMemcpy(mat->deviceData, mat->hostData, sizeof(PetscScalar) * size, cudaMemcpyHostToDevice));


  cudaDeviceSynchronize();
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
  PetscInt rank;
  PetscInt numprocs;

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);

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

  diag_size = (mat->m * mat->n * mat->p * mat->dof * mat->dof)/numprocs;
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
    if (mat->dim == 2) 
    {
      (*mat->diag_starts)[mat->m] = 3 * diag_size;
      (*mat->diagonals).push_back(mat->m);
      (*mat->diag_starts)[-mat->m] = 4 * diag_size;
      (*mat->diagonals).push_back(-mat->m);

    } 
    else if (mat->dim == 3) 
    {

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

  vecsize = (mat->m * mat->n * mat->p * mat->dof)/numprocs;

  // We know the expected size of x, y, so go ahead and allocate them now
  checkCudaError(cudaMalloc(&mat->deviceX, vecsize * sizeof(PetscScalar)));
  checkCudaError(cudaMalloc(&mat->deviceY, vecsize * sizeof(PetscScalar)));

  // We also know how many diagonals we have, and their indices
  checkCudaError(cudaMalloc(&mat->deviceDiags, sizeof(int) * mat->diagonals->size()));
  A->preallocated = PETSC_TRUE;
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
//  MatSetUpMultiply_MPISGGPU(A);

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

  PetscInt numprocs, rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);

  nrows = (a->m * a->n * a->p * a->dof)/numprocs;
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
	fprintf(stdout,"row %2.2d:",rank*nrows + i); 

	    for (j=0;j<ndiag*dof;j++) 
		{
			fprintf(stdout," %4g ",a->hostData[i+j*nrows]);
		}
	fprintf(stdout,"\n");
  	}
	fprintf(stdout,"\n\n");


  for (iblock=0;iblock<(nrows/dof);iblock++)  
	{
	    for (i=iblock*dof;i<(iblock+1)*dof;i++) 
		{
		fprintf(stdout,"row %d:",rank*nrows + i);
	
		  for (std::map<int, int>::iterator I = a->diag_starts->begin(),
        	     E = a->diag_starts->end(); I != E; ++I) 
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
	                  fprintf(stdout," (%d, %g) ",col,a->hostData[index]);
        		}
      	 	  }
		fprintf(stdout,"\n");
    		}
  	}

}