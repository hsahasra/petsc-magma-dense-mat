/// SGGPU Matrix Type

#define PETSCMAT_DLL

// Debugging flags
#define _TRACE 0
#define _TIME 0
#define _CSV_OUT 0

#if _TRACE
#define SGTrace printf("[SeqSGGPU] %s\n",__FUNCT__);
#else
#define SGTrace
#endif

#include "petsc-private/matimpl.h"

// Direct access to seqgpu vector type
#include "../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"

// Interop with CUSP vector
#include "../src/vec/vec/impls/seq/seqcusp/cuspvecimpl.h"

#include <stdio.h>
#include <cuda.h>
#include <sys/time.h>

// C++ library headers
#include <map>

#include "sggpu.h"

// Hard-coded block size
#define BLOCKWIDTH_X 128
#define BLOCKWIDTH_Y 1


// ----------------------------------------------------------
// helper function for error checking
// pops the CUDA error stack and exits on nonzero error code
// written by: dlowell ANL-MCS
// ----------------------------------------------------------
EXTERN_C_BEGIN
void checkCudaError(cudaError_t err) {
  if(cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s.\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}
EXTERN_C_END
// ----------------------------------------------------------



//------------------------------------------------------
// general timer function using unix system call
// dlowell ANL-MCS
//------------------------------------------------------
double getclock() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}




//===-- CUDA Device Code -------------------------------------------------===//
 
texture<int2, 1> vector_x;
     
static __inline__ __device__ double fetch_double(texture<int2, 1> tex, int i)
     {
       int2 v = tex1Dfetch(tex, i);
       return __hiloint2double(v.y, v.x);
     }
     
__global__ void MatMultKernel(PetscScalar * coeff, PetscScalar * y, PetscInt mat_size, PetscInt num_diags, int * diagonals, PetscInt dof) {
       
int idx = blockDim.x * blockIdx.x * 1 + threadIdx.x * 1;
     
if (idx >= mat_size)
      return;
     
int diag_size = mat_size * dof;
     
PetscScalar yval0 = 0.0;
int idx0 = idx;
     
//#pragma unroll 4
for (int i = 0; i < num_diags; ++i) {
    int d = diagonals[i];
     
    int offset0 = diag_size * i + idx0;
    int block0 = (idx0 / dof + d) * dof;
     
    //#pragma unroll 12
    for (int j = 0; j < dof; ++j) {
      // Get coefficient
      PetscScalar aval0 = coeff[offset0 + mat_size*j];
      // Get X value
      PetscScalar xval0 = fetch_double(vector_x, block0 + j);
     
      yval0 += aval0 * xval0;
           }
         }
   
      y[idx0] = yval0;
      }
    
//===-- Host Code --------------------------------------------------------===//




// Matrix function table
static struct _MatOps MatOps_Values = {
/*0*/ MatSetValues_SeqSGGPU,MatGetRow_SeqSGGPU,MatRestoreRow_SeqSGGPU,MatMult_SeqSGGPU,0,
/*5*/0,0,0,0,0,
/*10*/0,0,0,0,0,
/*15*/0,0,MatGetDiagonal_SeqSGGPU,MatDiagonalScale_SeqSGGPU,0,
/*20*/MatAssemblyBegin_SeqSGGPU,MatAssemblyEnd_SeqSGGPU,0,MatZeroEntries_SeqSGGPU,0,
/*25*/0,0,0,0,MatSetUp_SeqSGGPU,
/*30*/0,0,0,0,0,
/*35*/0,0,0,0,0,
/*40*/0,0,0,0,0,
/*45*/0,0,0,0,0,
/*50*/0,0,MatGetColumnIJ_SeqSGGPU,0,MatFDColoringCreate_SeqSGGPU,
/*55*/0,0,0,MatSetValuesBlocked_SeqSGGPU,0,
/*60*/MatDestroy_SeqSGGPU,MatView_SeqSGGPU,0,0,0,
/*65*/0,0,MatSetValues_SeqSGGPU,0,MatGetRowMaxAbs_SeqSGGPU,
/*70*/0,0,0,0,0,
/*75*/MatFDColoringApply_SeqSGGPU,0,0,0,0,
/*80*/0,0,0,0,0,
/*85*/0,0,MatSetValuesBlocked_SeqSGGPU,0,0,
/*90*/0,0,0,0,0,
/*95*/0,0,0,0,0,
/*100*/0,0,0,0,0,
/*105*/0,0,0,0,0,
/*110*/0,0,0,0,0,
/*115*/MatCreate_SeqSGGPU,0,0,0,0,
/*120*/0,0,0,0,0,
/*125*/0,0,0,0,0,
/*130*/0,0,0,0,0,
/*135*/0,0,0,0,0,
/*140*/0,0,
/*142*/MatSetGrid_SeqSGGPU
};




EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatCreate_SeqSGGPU"
PetscErrorCode MatCreate_SeqSGGPU(Mat A)
{
  Mat_SeqSGGPU * mat;
  PetscErrorCode ierr;
  PetscMPIInt size;

  PetscFunctionBegin;
  SGTrace;

	PetscPrintf(PETSC_COMM_WORLD,"MatCreate_SeqSGGPU\n");
	

  ierr = MPI_Comm_size(((PetscObject)A)->comm, &size); CHKERRQ(ierr);
  if (size > 1)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Comm must be size 1");

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
  ierr = PetscObjectChangeTypeName((PetscObject)A, MATSEQSGGPU); CHKERRQ(ierr);

  ierr = PetscObjectComposeFunctionDynamic((PetscObject)A,
        "MatSeqSGGPUSetPreallocation_C","MatSeqSGGPUSetPreallocation_SeqDIA",
        MatSeqSGGPUSetPreallocation_SeqSGGPU);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
EXTERN_C_END


#undef __FUNCT__
#define __FUNCT__ "MatDestroy_SeqSGGPU"
PetscErrorCode MatDestroy_SeqSGGPU(Mat A)
{
  Mat_SeqSGGPU * mat;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  SGTrace;

  mat = (Mat_SeqSGGPU*)A->data;

  if (mat->hostData) {
    ierr = PetscFree(mat->hostData); CHKERRQ(ierr);
  }
  if (mat->deviceData) {
    cudaFree(mat->deviceData);
  }
  if (mat->diag_starts) {
    delete mat->diag_starts;
  }
  ierr = PetscFree(mat->diag_offsets); CHKERRQ(ierr);
  if (mat->diagonals) {
    delete mat->diagonals;
  }
  if (mat->deviceX) {
    cudaFree(mat->deviceX);
  }
  if (mat->deviceY) {
    cudaFree(mat->deviceY);
  }
  if (mat->deviceDiags) {
    cudaFree(mat->deviceDiags);
  }
  if(mat->ja)       { ierr = PetscFree(mat->ja); CHKERRQ(ierr);       }
  if(mat->ia)       { ierr = PetscFree(mat->ia); CHKERRQ(ierr);       }
  checkCudaError(cudaStreamDestroy(mat->stream));
  PetscFree(mat); CHKERRQ(ierr);

  ierr = PetscObjectChangeTypeName((PetscObject)A, 0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetGrid_SeqSGGPU"
PetscErrorCode MatSetGrid_SeqSGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p)
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
#define __FUNCT__ "MatMult_SeqSGGPU"
PetscErrorCode MatMult_SeqSGGPU(Mat A, Vec x, Vec y)
{
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;
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
    dim3 grid((int)ceil((float)(mat->m * mat->n * mat->p * mat->dof)/(float)BLOCKWIDTH_X / 1.0), 1);

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
    mat_size = mat->m * mat->n * mat->p * mat->dof;

    checkCudaError(cudaBindTexture(0, vector_x, devX, mat_size * sizeof(PetscScalar)));
    MatMultKernel<<<grid, block, shared_size, mat->stream>>>(mat->deviceData, devY, mat_size, mat->diagonals->size(), mat->deviceDiags, mat->dof);

    cudaUnbindTexture(vector_x);
    cudaDeviceSynchronize();    

  } else if (iscusp) {
    dim3 block(BLOCKWIDTH_X, BLOCKWIDTH_Y);
    dim3 grid((int)ceil((float)(mat->m * mat->n * mat->p * mat->dof)/(float)BLOCKWIDTH_X / 1.0), 1);

    int shared_size = 0;
    /* Make sure y is also VECCUSP */
    ierr = PetscObjectTypeCompare((PetscObject)x,VECCUSP,&isseqgpu);CHKERRQ(ierr);
    if (!iscusp) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and y must be same type");
    }

    mat_size = mat->m * mat->n * mat->p * mat->dof;
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
      MatMultKernel<<<grid, block, shared_size, mat->stream>>>(mat->deviceData, devY, mat_size, mat->diagonals->size(), mat->deviceDiags, mat->dof);
#if _TIME
    checkCudaError(cudaStreamSynchronize(mat->stream));
    end = getclock();
    double elapsed = end - start;
    double gflops = (2.0 * mat->non_zeros / elapsed / 1e9);

    double nos = ((mat->p == 1 ? 2 : 3) * 2 + 1) * (2*mat->dof - 1);
    double nz = mat->m * mat->n * mat->p * mat->dof;
    double alt_gflops = (2.0 * nos * nz) / ((end - start)*1024*1024*1024);

#if _CSV_OUT
    fprintf(stderr, "%d,%d,%d,%d,%lf,%lf,\n", mat->m, mat->n, mat->p, mat->dof, elapsed, gflops);
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

	VecView(x,PETSC_VIEWER_STDOUT_WORLD);
	VecView(y,PETSC_VIEWER_STDOUT_WORLD);


  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_SeqSGGPU"
PetscErrorCode MatSetValuesBlocked_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatSetValuesBlocked_SeqSGGPU not implemented");

}


#undef __FUNCT__
#define __FUNCT__ "MatSetValues_SeqSGGPU"
PetscErrorCode MatSetValues_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
  int i, j;
  PetscErrorCode ierr;
  PetscBool resizegpu = PETSC_FALSE;
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;

  PetscFunctionBegin;
  SGTrace;

	//PetscPrintf(PETSC_COMM_WORLD,"MatSetValues_SeqSGGPU\n");

  // Handle each element
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      // Compute the diagonal and offset into the diagonal storage
      // for the element
      int row = irow[i];
      int col = icol[j];
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


#undef __FUNCT__  
#define __FUNCT__ "MatFDColoringView_Private"
PetscErrorCode MatFDColoringView_Private(MatFDColoring fd)
{
  PetscErrorCode ierr;
  PetscBool      flg = PETSC_FALSE;
  PetscViewer    viewer;

  PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"MatFDColoringView_Private\n");


  ierr = PetscViewerASCIIGetStdout(((PetscObject)fd)->comm,&viewer);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_fd_coloring_view",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = MatFDColoringView(fd,viewer);CHKERRQ(ierr);
  }
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_fd_coloring_view_info",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
    ierr = MatFDColoringView(fd,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  }
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_fd_coloring_view_draw",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = MatFDColoringView(fd,PETSC_VIEWER_DRAW_(((PetscObject)fd)->comm));CHKERRQ(ierr);
    ierr = PetscViewerFlush(PETSC_VIEWER_DRAW_(((PetscObject)fd)->comm));CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatSetUp_SeqSGGPU"
PetscErrorCode MatSetUp_SeqSGGPU(Mat A)
{

  PetscFunctionBegin;
  SGTrace;

	PetscPrintf(PETSC_COMM_WORLD,"MatSetUP_SeqSGGPU\n");


  //  ierr =  MatSeqSGGPUSetPreallocation(A,PETSC_DEFAULT,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatSeqSGGPUSetPreallocation"
PetscErrorCode MatSeqSGGPUSetPreallocation(Mat A,PetscInt stencil_type, PetscInt dof)
{
  PetscErrorCode ierr;
  Mat_SeqSGGPU *mat = (Mat_SeqSGGPU*)A->data;
  PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"MatSeqSGGPUSetPreallocation\n");


  mat->stencil_type = stencil_type;
  mat->dof = dof;
  if(A->preallocated)PetscFunctionReturn(0);
  PetscValidHeaderSpecific(A,MAT_CLASSID,1);
  
  ierr = PetscTryMethod(A,"MatSeqSGGPUSetPreallocation_C",(Mat,PetscInt,const PetscInt []),(A,0,0));CHKERRQ(ierr);
  A->preallocated=PETSC_TRUE;
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatSeqSGGPUSetPreallocation_SeqSGGPU"
extern PetscErrorCode MatSeqSGGPUSetPreallocation_SeqSGGPU(Mat A,PetscInt nz, const PetscInt nnz[])
{
  PetscErrorCode ierr;
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;
  PetscInt dim,diag_size,size,num_diags,i,vecsize;

	PetscPrintf(PETSC_COMM_WORLD,"MateqSGGPUSetPreallocation_SeqSGGPU\n");


  ierr = PetscLayoutSetBlockSize(A->rmap,1);CHKERRQ(ierr);
  ierr = PetscLayoutSetBlockSize(A->cmap,1);CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(A->rmap);CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(A->cmap);CHKERRQ(ierr);

  dim = A->stencil.dim;
  if (mat->dof > 1) {
    dim--;
  }
  mat->m = mat->n = mat->p = 1;
  mat->dim = dim;
  if (mat->dim > 0) mat->m = A->stencil.dims[dim-1];
  if (mat->dim > 1) mat->n = A->stencil.dims[dim-2];
  if (mat->dim > 2) mat->p = A->stencil.dims[dim-3];

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


#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_SeqSGGPU"
PetscErrorCode MatZeroEntries_SeqSGGPU(Mat A)
{
  Mat_SeqSGGPU *mat = (Mat_SeqSGGPU*)A->data;
  PetscInt size;
  PetscFunctionBegin;
  SGTrace;
  
  size = mat->diag_starts->size() * (mat->m * mat->n * mat->p * mat->dof * mat->dof);
  memset(mat->hostData, 0, size * sizeof(PetscScalar));
  
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_SeqSGGPU"
PetscErrorCode MatGetDiagonal_SeqSGGPU(Mat A, Vec v)
{
  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatGetDiagonal_SeqSGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_SeqSGGPU"
PetscErrorCode MatDiagonalScale_SeqSGGPU(Mat A, Vec ll, Vec rr)
{
  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatDiagonalScale_SeqSGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRow_SeqSGGPU"
PetscErrorCode MatGetRow_SeqSGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v)
{
  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatGetRow_SeqSGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_SeqSGGPU"
PetscErrorCode MatRestoreRow_SeqSGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v)
{
  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatRestoreRow_SeqSGGPU not implemented");
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_SeqSGGPU"
PetscErrorCode MatGetRowMaxAbs_SeqSGGPU(Mat A, Vec v, PetscInt idx[])
{
  PetscFunctionBegin;
  SGTrace;
  SETERRQ(PETSC_COMM_SELF,0,"MatGetRowMaxAbs_SeqSGGPU not implemented");
}

#undef __FUNCT__
#define __FUNCT__ "MatView_SeqSGGPU_ASCII"
PetscErrorCode MatView_SeqSGGPU_ASCII(Mat A, PetscViewer viewer)
{
  Mat_SeqSGGPU *a = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode ierr;
  PetscInt nrows,ndiag,dof,i,j,iblock,col,index,offset;
  std::map<int, int> &diag_starts = *(a->diag_starts);
  

  PetscFunctionBegin;  



  cudaDeviceSynchronize();
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscObjectPrintClassNamePrefixType((PetscObject)A,viewer,"Matrix Object");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"MatView_SeqSGGPU_ASCII still in development\n");

  nrows = a->m * a->n * a->p * a->dof;
  ndiag = a->diagonals->size();
  dof = a->dof;

  ierr = PetscViewerASCIIPrintf(viewer,"offsets: \n"); CHKERRQ(ierr);
  for (std::map<int, int>::iterator I = diag_starts.begin(),
         E = diag_starts.end(); I != E; ++I) {
    PetscViewerASCIIPrintf(viewer,"- Diag %d:%d\n", I->first, I->second);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n"); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"hostData:\n"); CHKERRQ(ierr);
  for (i=0;i<nrows;i++) {
    ierr = PetscViewerASCIIPrintf(viewer,"row %2.2D:",i); CHKERRQ(ierr);
    for (j=0;j<ndiag*dof;j++) {
      ierr = PetscViewerASCIIPrintf(viewer," %4G ",a->hostData[i+j*nrows]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n\n");



  for (iblock=0;iblock<nrows/dof;iblock++)  {
    for (i=iblock*dof;i<(iblock+1)*dof;i++) {
      ierr = PetscViewerASCIIPrintf(viewer,"row %D:",i);CHKERRQ(ierr);
      for (std::map<int, int>::iterator I = a->diag_starts->begin(),
             E = a->diag_starts->end(); I != E; ++I) {
        /* Ignore 0 padding */
        offset = I->first;
        if (offset + iblock < 0) {
          continue;
        }
        if (offset + iblock >= (nrows/dof)) {
          break;
        }
        
        for (j=0;j<dof;j++) {
          index = i + I->second + j*nrows; // column-major
          col = offset*dof+(iblock*dof) + j;
#if defined(PETSC_USE_COMPLEX)
          if (PetscImaginaryPart(a->hostData[index]) > 0.0) {
            ierr = PetscViewerASCIIPrintf(viewer," (%D, %G + %G i)",col,PetscRealPart(a->hostData[index]),PetscImaginaryPart(a->hostData[index]));CHKERRQ(ierr);
          } else if (PetscImaginaryPart(a->hostData[index]) < 0.0) {
            ierr = PetscViewerASCIIPrintf(viewer," (%D, %G - %G i)",col,PetscRealPart(a->hostData[index]),-PetscImaginaryPart(a->hostData[index]));CHKERRQ(ierr);
          } else {
            ierr = PetscViewerASCIIPrintf(viewer," (%D, %G) ",col,PetscRealPart(a->hostData[index]));CHKERRQ(ierr);
          }
#else
          ierr = PetscViewerASCIIPrintf(viewer," (%D, %G) ",col,a->hostData[index]);CHKERRQ(ierr);
#endif
        }
      }
      ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    }
  }

  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatView_SeqSGGPU"
PetscErrorCode MatView_SeqSGGPU(Mat A, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool isascii;
  PetscFunctionBegin;

  SGTrace;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = MatView_SeqSGGPU_ASCII(A,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatAssemblyBegin_SeqSGGPU"
PetscErrorCode MatAssemblyBegin_SeqSGGPU(Mat A, MatAssemblyType type)
{
  PetscFunctionBegin;
  SGTrace;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyEnd_SeqSGGPU"
PetscErrorCode MatAssemblyEnd_SeqSGGPU(Mat A, MatAssemblyType type)
{
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;
  PetscInt size;
  PetscFunctionBegin;

#if _TRACE
  printf("[SeqSGGPU] MatAssemblyEnd_SeqSGGPU\n");

  for (std::map<int, int>::iterator I = mat->diag_starts->begin(),
       E = mat->diag_starts->end(); I != E; ++I) {
    printf("- Diag %d:\n", I->first);
    for (int i = 0; i < mat->dof; ++i) {
      for (int j = 0; j < mat->dof * mat->m * mat->n * mat->p; ++j) {
        int offset = i * mat->dof * mat->m * mat->n * mat->p + j;
        printf(" %lf ", mat->hostData[offset + I->second]);
      }
      printf("\n");
    }
  }
#endif
  size = mat->diag_starts->size()*mat->m*mat->n*mat->p*mat->dof*mat->dof;

  checkCudaError(cudaMemcpyAsync(mat->deviceDiags, &(*mat->diagonals)[0], sizeof(int) * mat->diagonals->size(), cudaMemcpyHostToDevice, mat->stream));

  checkCudaError(cudaMemcpy(mat->deviceData, mat->hostData, sizeof(PetscScalar) * size, cudaMemcpyHostToDevice));

  cudaDeviceSynchronize();
  PetscFunctionReturn(0);

}



#undef __FUNCT__  
#define __FUNCT__ "MatFDColoringApply_SeqSGGPU"
PetscErrorCode  MatFDColoringApply_SeqSGGPU(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx)
{
  PetscErrorCode (*f)(void*,Vec,Vec,void*) = (PetscErrorCode (*)(void*,Vec,Vec,void *))coloring->f;
  PetscErrorCode ierr;
  PetscInt       k,start,end,l,row,col,srow,**vscaleforrow;
  PetscScalar    dx,*y,*xx,*w3_array;
  PetscScalar    *vscale_array;
  PetscReal      epsilon = coloring->error_rel,umin = coloring->umin,unorm; 
  Vec            w1,w2,w3;
  void           *fctx = coloring->fctx;
  PetscBool      flg = PETSC_FALSE;
  PetscInt       ctype=coloring->ctype,N,col_start=0,col_end=0;
  Vec            x1_tmp;

  PetscFunctionBegin;    

	PetscPrintf(PETSC_COMM_WORLD,"MatFDColoringApply_SeqSGGPU\n");


  PetscValidHeaderSpecific(J,MAT_CLASSID,1);
  PetscValidHeaderSpecific(coloring,MAT_FDCOLORING_CLASSID,2);
  PetscValidHeaderSpecific(x1,VEC_CLASSID,3);
  if (!f) SETERRQ(((PetscObject)J)->comm,PETSC_ERR_ARG_WRONGSTATE,"Must call MatFDColoringSetFunction()");

  ierr = PetscLogEventBegin(MAT_FDColoringApply,coloring,J,x1,0);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_fd_coloring_dont_rezero",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscInfo(coloring,"Not calling MatZeroEntries()\n");CHKERRQ(ierr);
  } else {
    PetscBool  assembled;
    ierr = MatAssembled(J,&assembled);CHKERRQ(ierr);
    if (assembled) {
      ierr = MatZeroEntries(J);CHKERRQ(ierr);
    }
  }

  x1_tmp = x1; 
  if (!coloring->vscale){ 
    ierr = VecDuplicate(x1_tmp,&coloring->vscale);CHKERRQ(ierr);
  }

  if (coloring->htype[0] == 'w') { /* tacky test; need to make systematic if we add other approaches to computing h*/
    ierr = VecNorm(x1_tmp,NORM_2,&unorm);CHKERRQ(ierr);
  }

  if (!coloring->w3) {
    /*
    ierr = VecDestroy(&coloring->w1); CHKERRQ(ierr);
    ierr = VecDestroy(&coloring->w2); CHKERRQ(ierr);
    ierr = VecDuplicate(x1_tmp,&coloring->w1);CHKERRQ(ierr);
     ierr = VecDuplicate(x1_tmp,&coloring->w2);CHKERRQ(ierr);
    ierr = PetscLogObjectParent(coloring,coloring->w1);CHKERRQ(ierr);
     ierr = PetscLogObjectParent(coloring,coloring->w2);CHKERRQ(ierr);*/
    ierr = VecDuplicate(x1_tmp,&coloring->w3);CHKERRQ(ierr);
    ierr = PetscLogObjectParent(coloring,coloring->w3);CHKERRQ(ierr);
  }
  w1 = coloring->w1;
  w2 = coloring->w2;
  w3 = coloring->w3;
  ierr = VecGetOwnershipRange(w1,&start,&end);CHKERRQ(ierr); /* OwnershipRange is used by ghosted x! */

  /* Set w1 = F(x1) */
  if (!coloring->fset) {
    ierr = PetscLogEventBegin(MAT_FDColoringFunction,0,0,0,0);CHKERRQ(ierr);
    ierr = (*f)(sctx,x1_tmp,w1,fctx);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(MAT_FDColoringFunction,0,0,0,0);CHKERRQ(ierr);
  } else {
    coloring->fset = PETSC_FALSE;
  }


    /* Compute all the local scale factors, including ghost points */
  ierr = VecGetLocalSize(x1_tmp,&N);CHKERRQ(ierr);
  ierr = VecGetArray(x1_tmp,&xx);CHKERRQ(ierr);
  ierr = VecGetArray(coloring->vscale,&vscale_array);CHKERRQ(ierr);
  if (ctype == IS_COLORING_GHOSTED){
    col_start = 0; col_end = N;
  } else if (ctype == IS_COLORING_GLOBAL){
    xx = xx - start;
    vscale_array = vscale_array - start;
    col_start = start; col_end = N + start;
  }
  for (col=col_start; col<col_end; col++){ 
    /* Loop over each local column, vscale[col] = 1./(epsilon*dx[col]) */      
    if (coloring->htype[0] == 'w') {
      dx = 1.0 + unorm;
    } else {
      dx  = xx[col];
    }
    if (dx == (PetscScalar)0.0) dx = 1.0;
#if !defined(PETSC_USE_COMPLEX)
    if (dx < umin && dx >= 0.0)      dx = umin;
    else if (dx < 0.0 && dx > -umin) dx = -umin;
#else
    if (PetscAbsScalar(dx) < umin && PetscRealPart(dx) >= 0.0)     dx = umin;
    else if (PetscRealPart(dx) < 0.0 && PetscAbsScalar(dx) < umin) dx = -umin;
#endif
    dx               *= epsilon;
    vscale_array[col] = (PetscScalar)1.0/dx;
  } 
  if (ctype == IS_COLORING_GLOBAL)  vscale_array = vscale_array + start;      
  ierr = VecRestoreArray(coloring->vscale,&vscale_array);CHKERRQ(ierr);
  if (ctype == IS_COLORING_GLOBAL){
    ierr = VecGhostUpdateBegin(coloring->vscale,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(coloring->vscale,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  }
    
  if (coloring->vscaleforrow) {
    vscaleforrow = coloring->vscaleforrow;
  } else SETERRQ(((PetscObject)J)->comm,PETSC_ERR_ARG_NULL,"Null Object: coloring->vscaleforrow");

  /*
    Loop over each color
  */
  ierr = VecGetArray(coloring->vscale,&vscale_array);CHKERRQ(ierr);
  for (k=0; k<coloring->ncolors; k++) { 
    coloring->currentcolor = k;
    ierr = VecCopy(x1_tmp,w3);CHKERRQ(ierr);
    ierr = VecGetArray(w3,&w3_array);CHKERRQ(ierr);
    if (ctype == IS_COLORING_GLOBAL) w3_array = w3_array - start;
    /*
      Loop over each column associated with color 
      adding the perturbation to the vector w3.
    */
    for (l=0; l<coloring->ncolumns[k]; l++) {
      col = coloring->columns[k][l];    /* local column of the matrix we are probing for */
      if (coloring->htype[0] == 'w') {
        dx = 1.0 + unorm;
      } else {
        dx  = xx[col];
      }
      if (dx == (PetscScalar)0.0) dx = 1.0;
#if !defined(PETSC_USE_COMPLEX)
      if (dx < umin && dx >= 0.0)      dx = umin;
      else if (dx < 0.0 && dx > -umin) dx = -umin;
#else
      if (PetscAbsScalar(dx) < umin && PetscRealPart(dx) >= 0.0)     dx = umin;
      else if (PetscRealPart(dx) < 0.0 && PetscAbsScalar(dx) < umin) dx = -umin;
#endif
      dx            *= epsilon;
      if (!PetscAbsScalar(dx)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Computed 0 differencing parameter");
      w3_array[col] += dx;

    } 
    if (ctype == IS_COLORING_GLOBAL) w3_array = w3_array + start;
    ierr = VecRestoreArray(w3,&w3_array);CHKERRQ(ierr);

    /*
      Evaluate function at w3 = x1 + dx (here dx is a vector of perturbations)
                           w2 = F(x1 + dx) - F(x1)
    */
    ierr = PetscLogEventBegin(MAT_FDColoringFunction,0,0,0,0);CHKERRQ(ierr);
    ierr = (*f)(sctx,w3,w2,fctx);CHKERRQ(ierr);        
    ierr = PetscLogEventEnd(MAT_FDColoringFunction,0,0,0,0);CHKERRQ(ierr);
    ierr = VecAXPY(w2,-1.0,w1);CHKERRQ(ierr); 
        
    /*
      Loop over rows of vector, putting results into Jacobian matrix
    */
    ierr = VecGetArray(w2,&y);CHKERRQ(ierr);
    for (l=0; l<coloring->nrows[k]; l++) {
      row    = coloring->rows[k][l];             /* local row index */
      col    = coloring->columnsforrow[k][l];    /* global column index */
      y[row] *= vscale_array[vscaleforrow[k][l]];
      srow   = row + start;
      ierr   = MatSetValues(J,1,&srow,1,&col,y+row,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(w2,&y);CHKERRQ(ierr);
  } /* endof for each color */
  if (ctype == IS_COLORING_GLOBAL) xx = xx + start; 
  ierr = VecRestoreArray(coloring->vscale,&vscale_array);CHKERRQ(ierr);
  ierr = VecRestoreArray(x1_tmp,&xx);CHKERRQ(ierr);
   
  coloring->currentcolor = -1;
  ierr  = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr  = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(MAT_FDColoringApply,coloring,J,x1,0);CHKERRQ(ierr);

  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_null_space_test",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = MatNullSpaceTest(J->nullsp,J,PETSC_NULL);CHKERRQ(ierr);
  }
  ierr = MatFDColoringView_Private(coloring);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetColumnIJ_SeqSGGPU"
PetscErrorCode MatGetColumnIJ_SeqSGGPU(Mat A,PetscInt oshift,PetscBool  symmetric,PetscBool  inodecompressed,PetscInt *nn, const PetscInt *ia[], const PetscInt *ja[],PetscBool  *done)
{
  Mat_SeqSGGPU     *a = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode ierr;
  PetscInt       n = A->cmap->n;
  PetscInt       ndiag = a->diagonals->size();
  PetscInt       nrows = a->m*a->n*a->p*a->dof;
  PetscInt       nz=a->dof*ndiag*nrows;
  PetscInt       iblock,i,j,col,index,colblock,offset;

  PetscFunctionBegin;  

	PetscPrintf(PETSC_COMM_WORLD,"MatGetColumnIJ_SeqSGGPU\n");


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

#undef __FUNCT__
#define __FUNCT__ "MatFDColoringCreate_SeqSGGPU"
PetscErrorCode MatFDColoringCreate_SeqSGGPU(Mat mat,ISColoring iscoloring,MatFDColoring c)
{
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

  N          = mat->cmap->N/bs;
  c->M       = mat->rmap->N/bs;  /* set total rows, columns and local rows */
  c->N       = mat->cmap->N/bs;
  c->m       = mat->rmap->N/bs;
  c->rstart  = 0;

  c->ncolors = nis;
  ierr       = PetscMalloc(nis*sizeof(PetscInt),&c->ncolumns);CHKERRQ(ierr);
  ierr       = PetscMalloc(nis*sizeof(PetscInt*),&c->columns);CHKERRQ(ierr); 
  ierr       = PetscMalloc(nis*sizeof(PetscInt),&c->nrows);CHKERRQ(ierr);
  ierr       = PetscMalloc(nis*sizeof(PetscInt*),&c->rows);CHKERRQ(ierr);
  ierr       = PetscMalloc(nis*sizeof(PetscInt*),&c->columnsforrow);CHKERRQ(ierr);

  ierr = MatGetColumnIJ(mat,0,PETSC_FALSE,PETSC_FALSE,&ncols,&ci,&cj,&done);CHKERRQ(ierr);
  if (!done) SETERRQ1(((PetscObject)mat)->comm,PETSC_ERR_SUP,"MatGetColumnIJ() not supported for matrix type %s",((PetscObject)mat)->type_name);

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
  ierr = MatRestoreColumnIJ(mat,0,PETSC_FALSE,PETSC_FALSE,&ncols,&ci,&cj,&done);CHKERRQ(ierr);

  ierr = PetscFree(rowhit);CHKERRQ(ierr);
  ierr = PetscFree(columnsforrow);CHKERRQ(ierr);

  /* Optimize by adding the vscale, and scaleforrow[][] fields */
  /*
       see the version for MPIAIJ
  */
  ierr = VecCreateGhost(((PetscObject)mat)->comm,mat->rmap->n,PETSC_DETERMINE,0,PETSC_NULL,&c->vscale);CHKERRQ(ierr);
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



