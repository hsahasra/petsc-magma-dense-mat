/// SGGPU Matrix Type

#define PETSCMAT_DLL

#include "petsc-private/matimpl.h"
#include "sggpu.h"

// Direct access to seqgpu vector type
#include "../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"

#include <stdio.h>
#include <cuda.h>
#include <sys/time.h>

// C++ library headers
#include <map>

// Debugging flags
#define _TRACE 0
#define _TIME 0
#define _CSV_OUT 0

// Hard-coded block size
#define BLOCKWIDTH_X 512
#define BLOCKWIDTH_Y 1

// Prototypes
PetscErrorCode MatDestroy_SeqSGGPU(Mat A);
PetscErrorCode MatSetGrid_SeqSGGPU(Mat B, PetscInt m, PetscInt n, PetscInt p);
PetscErrorCode MatMult_SeqSGGPU(Mat mat, Vec x, Vec y);
PetscErrorCode MatSetValuesBlocked_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetValues_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is);
PetscErrorCode MatSetStencil_SeqSGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof);
PetscErrorCode MatSetUpPreallocation_SeqSGGPU(Mat mat);
PetscErrorCode MatZeroEntries_SeqSGGPU(Mat A);
PetscErrorCode MatGetDiagonal_SeqSGGPU(Mat A, Vec v);
PetscErrorCode MatDiagonalScale_SeqSGGPU(Mat A, Vec ll, Vec rr);
PetscErrorCode MatGetRow_SeqSGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v);
PetscErrorCode MatRestoreRow_SeqSGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v);
PetscErrorCode MatGetRowMaxAbs_SeqSGGPU(Mat A, Vec v, PetscInt idx[]);
PetscErrorCode MatView_SeqSGGPU(Mat A, PetscViewer viewer);
PetscErrorCode MatAssemblyBegin_SeqSGGPU(Mat A, MatAssemblyType type);
PetscErrorCode MatAssemblyEnd_SeqSGGPU(Mat A, MatAssemblyType type);


// ----------------------------------------------------------
// helper function for error checking
// pops the CUDA error stack and exits on nonzero error code
// written by: dlowell ANL-MCS
// ----------------------------------------------------------
void checkCudaError(cudaError_t err) {
  if(cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s.\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

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

extern __shared__ int shared_diagonals[];

__global__ void MatMultKernel(PetscScalar * coeff, PetscScalar * y, PetscInt mat_size, PetscInt num_diags, int * diagonals, PetscInt dof) {
  int idx = blockDim.x * blockIdx.x * 1 + threadIdx.x * 1;

  if (idx >= mat_size)
    return;

  if (threadIdx.x < num_diags) {
    shared_diagonals[threadIdx.x] = diagonals[threadIdx.x];
  }
  __syncthreads();

#if _TRACE
  if (threadIdx.x == 0) {
    printf("Diagonals:\n");
    for (int i = 0; i < num_diags; ++i) {
      printf("- %d\n", diagonals[i]);
    }
    printf("Foo: %d\n", -42);
  }
#endif

  int diag_size = mat_size * dof;

  PetscScalar yval0 = 0.0;
  //PetscScalar yval1 = 0.0;
  //PetscScalar yval2 = 0.0;
  //PetscScalar yval3 = 0.0;

  int idx0 = idx;
  //int idx1 = idx + 1;
  //int idx2 = idx + 2;
  //int idx3 = idx + 3;

  #pragma unroll 4
  for (int i = 0; i < num_diags; ++i) {
    int d = shared_diagonals[i];
    int offset0 = diag_size * i + idx0;
    int block0 = (idx0 / dof + d) * dof;
    //int offset1 = diag_size * i + idx1;
    //int block1 = (idx1 / dof + d) * dof;
    //int offset2 = diag_size * i + idx2;
    //int block2 = (idx2 / dof + d) * dof;
    //int offset3 = diag_size * i + idx3;
    //int block3 = (idx3 / dof + d) * dof;

    #pragma unroll 6
    for (int j = 0; j < dof; ++j) {
      // Get coefficient
#if _TRACE
      if (threadIdx.x == 0) {
        printf("diag: %d  offset: %d  index: %d\n", d, offset, offset + mat_size*j);
      }
#endif

      PetscScalar aval0 = coeff[offset0 + mat_size*j];
      //PetscScalar aval1 = coeff[offset1 + mat_size*j];
      //PetscScalar aval2 = coeff[offset2 + mat_size*j];
      //PetscScalar aval3 = coeff[offset3 + mat_size*j];

#if _TRACE
      if (threadIdx.x == 0) {
        printf("aval: %lf\n", aval);
      }
#endif

      // Get x value
      int this_block0 = block0 + j;
      //int this_block1 = block1 + j;
      //int this_block2 = block2 + j;
      //int this_block3 = block3 + j;

#if _TRACE
      if (threadIdx.x == 0) {
        printf("this_block: %d\n", this_block);
      }
#endif

      /*bool in_bounds0 = this_block0 >= 0 && this_block0 < mat_size;
      bool in_bounds1 = this_block1 >= 0 && this_block1 < mat_size;
      bool in_bounds2 = this_block2 >= 0 && this_block2 < mat_size;
      bool in_bounds3 = this_block3 >= 0 && this_block3 < mat_size;

      PetscScalar xval0 = in_bounds0 ? fetch_double(vector_x, this_block0) : 0.0;
      PetscScalar xval1 = in_bounds1 ? fetch_double(vector_x, this_block1) : 0.0;
      PetscScalar xval2 = in_bounds2 ? fetch_double(vector_x, this_block2) : 0.0;
      PetscScalar xval3 = in_bounds3 ? fetch_double(vector_x, this_block3) : 0.0;*/

      PetscScalar xval0 = fetch_double(vector_x, this_block0);
      //PetscScalar xval1 = fetch_double(vector_x, this_block1);
      //PetscScalar xval2 = fetch_double(vector_x, this_block2);
      //PetscScalar xval3 = fetch_double(vector_x, this_block3);

#if _TRACE
      if (threadIdx.x == 0) {
        printf("xval: %lf\n", xval);
      }
#endif

      yval0 += aval0 * xval0;
      //yval1 += aval1 * xval1;
      //yval2 += aval2 * xval2;
      //yval3 += aval3 * xval3;
    }
  }

  y[idx0] = yval0;
  //y[idx1] = yval1;
  //y[idx2] = yval2;
  //y[idx3] = yval3;
}

//===-- Host Code --------------------------------------------------------===//


// Matrix function table
static struct _MatOps MatOps_Values = {
/*0*/ MatSetValues_SeqSGGPU,MatGetRow_SeqSGGPU,MatRestoreRow_SeqSGGPU,MatMult_SeqSGGPU,0,
/*5*/0,0,0,0,0,
/*10*/0,0,0,0,0,
/*15*/0,0,MatGetDiagonal_SeqSGGPU,MatDiagonalScale_SeqSGGPU,0,
/*20*/MatAssemblyBegin_SeqSGGPU,MatAssemblyEnd_SeqSGGPU,0,MatZeroEntries_SeqSGGPU,0,
/*25*/0,0,0,0,MatSetUpPreallocation_SeqSGGPU,
/*30*/0,0,0,0,0,
/*35*/0,0,0,0,0,
/*40*/0,0,0,0,0,
/*45*/0,0,0,0,0,
/*50*/0,0,0,0,0,
/*55*/0,0,0,MatSetValuesBlocked_SeqSGGPU,0,
/*60*/MatDestroy_SeqSGGPU,MatView_SeqSGGPU,0,0,0,
/*65*/0,0,MatSetValues_SeqSGGPU,0,MatGetRowMaxAbs_SeqSGGPU,
/*70*/0,0,0,0,0,
/*75*/0,0,0,0,0,
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
/*135*/0,0,0,0,MatSetStencil_SeqSGGPU,
/*140*/MatSetGrid_SeqSGGPU
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

#if _TRACE
  printf("[SeqSGGPU] MatCreate_SeqSGGPU\n");
#endif

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

#if _TRACE
  printf("[SeqSGGPU] MatDestroy_SeqSGGPU\n");
#endif

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
#if _TRACE
  printf("[SeqSGGPU] MatSetGrid_SeqSGGPU\n");
#endif

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
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatMult_SeqSGGPU\n");
#endif

  // Initialize y to zero
  ierr = VecSet(y, 0.0); CHKERRQ(ierr);

  // NOTE: The seqgpu vector type is not really working here...

  //const VecType vec_type;

  // Get device pointer for X
  /*ierr = PetscObjectGetType((PetscObject)x, &vec_type); CHKERRQ(ierr);
  if (!strcmp(vec_type, "seqgpu")) {
    // We have a GPU vector type, so just use the existing pointer
    deviceX = ((Vec_SeqGPU*)x->data)->devptr;
  } else {
    fprintf(stderr, "Non-GPU vector types are not implemented!");
    exit(EXIT_FAILURE);
  }

  // Get device pointer for Y
  ierr = PetscObjectGetType((PetscObject)y, &vec_type); CHKERRQ(ierr);
  if (!strcmp(vec_type, "seqgpu")) {
    // We have a GPU vector type, so just use the existing pointer
    deviceY = ((Vec_SeqGPU*)y->data)->devptr;
  } else {
    fprintf(stderr, "Non-GPU vector types are not implemented!");
    exit(EXIT_FAILURE);
  }*/

  PetscScalar * hostX;
  PetscScalar * hostY;

  int mat_size = mat->m * mat->n * mat->p * mat->dof;

  ierr = VecGetArray(x, &hostX); CHKERRQ(ierr);
  ierr = VecGetArray(y, &hostY); CHKERRQ(ierr);

  checkCudaError(cudaMemcpyAsync(mat->deviceX, hostX, mat_size * sizeof(PetscScalar), cudaMemcpyHostToDevice, mat->stream));


  // Bind X to device texture
  checkCudaError(cudaBindTexture(0, vector_x, mat->deviceX, mat_size * sizeof(PetscScalar)));

#if _TRACE
  printf("Host diagonals:\n");
  for (int i = 0; i < mat->diagonals->size(); ++i) {
    printf("- %d\n", (*mat->diagonals)[i]);
  }
#endif

  // Invoke
  dim3 block(BLOCKWIDTH_X, BLOCKWIDTH_Y);
  dim3 grid((int)ceil((float)(mat->m * mat->n * mat->p * mat->dof)/(float)BLOCKWIDTH_X / 1.0), 1);

  int shared_size = mat->diagonals->size() * sizeof(int);

  //cudaFuncSetCacheConfig(MatMultKernel, cudaFuncCachePreferL1);


#if _TIME
  double start, end;
  start = getclock();
#endif
  MatMultKernel<<<grid, block, shared_size, mat->stream>>>(mat->deviceData, mat->deviceY, mat_size, mat->diagonals->size(), mat->deviceDiags, mat->dof);
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

  checkCudaError(cudaMemcpyAsync(hostY, mat->deviceY, mat_size * sizeof(PetscScalar), cudaMemcpyDeviceToHost, mat->stream));

  // Cleanup
  cudaUnbindTexture(vector_x);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValuesBlocked_SeqSGGPU"
PetscErrorCode MatSetValuesBlocked_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetValuesBlocked_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetValues_SeqSGGPU"
PetscErrorCode MatSetValues_SeqSGGPU(Mat A, PetscInt nrow, const PetscInt irow[], PetscInt ncol, const PetscInt icol[], const PetscScalar y[], InsertMode is)
{
  int i, j;
  PetscErrorCode ierr;
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;

  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetValues_SeqSGGPU\n");
#endif

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

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetStencil_SeqSGGPU"
PetscErrorCode MatSetStencil_SeqSGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof)
{
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetStencil_SeqSGGPU  (%p)\n", A);
#endif

  if (dim < 1 || dim > 3) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Dim must be between 1 and 3.");
  }

  mat->m = dims[0];
  if (dim > 1) {
    mat->n = dims[1];
    if (dim > 2) {
      mat->p = dims[2];
    } else {
      mat->p = 1;
    }
  } else {
    mat->n = 1;
    mat->p = 1;
  }

  mat->dof = dof;

#if _TRACE
  printf("- m: %d  n: %d  p: %d  dof: %d\n", mat->m, mat->n, mat->p, mat->dof);
#endif

  // It appears that we are responsible for pre-allocating
  if (!A->preallocated) {
    ierr = MatSetUpPreallocation_SeqSGGPU(A); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSetUpPreallocation_SeqSGGPU"
PetscErrorCode MatSetUpPreallocation_SeqSGGPU(Mat A)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatSetUpPreallocation_SeqSGGPU\n");
#endif

  A->preallocated = PETSC_TRUE;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatZeroEntries_SeqSGGPU"
PetscErrorCode MatZeroEntries_SeqSGGPU(Mat A)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatZeroEntries_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_SeqSGGPU"
PetscErrorCode MatGetDiagonal_SeqSGGPU(Mat A, Vec v)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatGetDiagonal_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatDiagonalScale_SeqSGGPU"
PetscErrorCode MatDiagonalScale_SeqSGGPU(Mat A, Vec ll, Vec rr)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatDiagonalScale_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRow_SeqSGGPU"
PetscErrorCode MatGetRow_SeqSGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatGetRow_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatRestoreRow_SeqSGGPU"
PetscErrorCode MatRestoreRow_SeqSGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatRestoreRow_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetRowMaxAbs_SeqSGGPU"
PetscErrorCode MatGetRowMaxAbs_SeqSGGPU(Mat A, Vec v, PetscInt idx[])
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatGetRowMaxAbs_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatView_SeqSGGPU"
PetscErrorCode MatView_SeqSGGPU(Mat A, PetscViewer viewer)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatView_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyBegin_SeqSGGPU"
PetscErrorCode MatAssemblyBegin_SeqSGGPU(Mat A, MatAssemblyType type)
{
  PetscFunctionBegin;
#if _TRACE
  printf("[SeqSGGPU] MatAssemblyBegin_SeqSGGPU\n");
#endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatAssemblyEnd_SeqSGGPU"
PetscErrorCode MatAssemblyEnd_SeqSGGPU(Mat A, MatAssemblyType type)
{
  Mat_SeqSGGPU * mat = (Mat_SeqSGGPU*)A->data;
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

  // Create GPU buffer
  if (mat->deviceData) {
    cudaFree(mat->deviceData);
  }
  int size = mat->diag_starts->size() * mat->m * mat->n * mat->p * mat->dof * mat->dof;
  checkCudaError(cudaMalloc(&mat->deviceData, sizeof(PetscScalar) * size));

  // Copy data to device
  checkCudaError(cudaMemcpy(mat->deviceData, mat->hostData, sizeof(PetscScalar) * size, cudaMemcpyHostToDevice));


  int mat_size = mat->m * mat->n * mat->p * mat->dof;

  // We know the expected size of x, y, so go ahead and allocate them now
  checkCudaError(cudaMalloc(&mat->deviceX, mat_size * sizeof(PetscScalar)));
  checkCudaError(cudaMalloc(&mat->deviceY, mat_size * sizeof(PetscScalar)));

  // We also know how many diagonals we have, and their indices
  checkCudaError(cudaMalloc(&mat->deviceDiags, sizeof(int) * mat->diagonals->size()));
  checkCudaError(cudaMemcpyAsync(mat->deviceDiags, &(*mat->diagonals)[0], sizeof(int) * mat->diagonals->size(), cudaMemcpyHostToDevice, mat->stream));

  PetscFunctionReturn(0);
}

