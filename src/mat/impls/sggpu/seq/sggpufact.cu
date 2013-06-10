/// SGGPU Matrix Type

#define PETSCMAT_DLL

#include "petsc-private/matimpl.h"
#include "../src/mat/impls/sggpu/seq/sggpu.h"
#include "../src/mat/impls/aij/seq/aij.h"  

// Direct access to seqgpu vector type
//#include "../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"

// Interop with CUSP vector
#include "../src/vec/vec/impls/seq/seqcusp/cuspvecimpl.h"

#include <stdio.h>
#include <cuda.h>
#include <sys/time.h>

// C++ library headers
#include <map>

// Hard-coded block size
#define BLOCKWIDTH_X 128
#define BLOCKWIDTH_Y 1

// ROW BLOCKING
#define USE_ROWBLOCKS 1
#define ROWBLOCK 1

// Debugging flags
#define _TRACE 0
#define _TIME 0

// Hard-coded sizes for cuda shared memory arrays
#define MAXBLOCKSIZE 256
#define BETAMAX 512
#define MAXNUMDIAGS  48

// max dof
#define MAXDOF 32

// macro to return element e of diagonal d in the
// array hd which has nd diagonals of size ds
// assuming dof==1
#if USE_ROWBLOCKS
#define TRI_MATRIX_ENTRY1(hd,d,e,nd,ds) hd[ (e)*((nd/2)+1) + (d) ]
#else
#define TRI_MATRIX_ENTRY1(hd,d,e,nd,ds) hd[ (d)*(ds) + (e) ]
#endif

// ----------------------------------------------------------
// helper functions from sggpu.cu
// ----------------------------------------------------------
EXTERN_C_BEGIN
void checkCudaError(cudaError_t err);
EXTERN_C_END
double getclock();

// ----------------------------------------------------------
// helper functions from sggpusse.c
// ----------------------------------------------------------
extern "C"
PetscErrorCode MatSolve_SeqSGGPU_sse4(int m, int n, int p, int dof, int dim, PetscScalar *hostData,
				      int numDiags, int *diagOffsets, PetscScalar *x, const PetscScalar *b);

extern "C"
PetscErrorCode MatSolve_SeqSGGPU_sse4_nochunk(int m, int n, int p, int dof, int dim,
					      PetscScalar *hostData,
					      int numDiags, int *diagOffsets,
					      PetscScalar *x, const PetscScalar *b);

PetscBool dumpMat;
PetscBool dumpVec;
int       exNumber;

//===-- CUDA Device Code -------------------------------------------------===//

texture<int2, 1> vector_rhs;

static __inline__ __device__ double fetch_double(texture<int2, 1> tex, int i)
{
  int2 v = tex1Dfetch(tex, i);
  return __hiloint2double(v.y, v.x);
}

__global__ void MatInvertDiagBlocks_Kernel(PetscScalar * coeff, PetscInt mat_size, PetscInt centerDiagIndex, PetscInt dof) {
  // idx is the index of a dof*dof block on the center diagonal of the matrix
  // the input parameter mat_size is the number of such blocks on the diagonal
  int idx = blockDim.x * blockIdx.x * 1 + threadIdx.x * 1;
  if (idx >= mat_size)
    return;

  // idx0 is the index within the subarray of elements that make up the center
  // diagonal of the [0,0] element of the dof*dof block this thread is inverting
  int idx0 = centerDiagIndex + idx * dof;
  //can we use __shared__ memory for this?
  PetscScalar block_inv[MAXBLOCKSIZE];

  int stripeSize = mat_size * dof;
  int i,j,k;

  // invert lower triangular matrix in block
  // diagonal elements are all 1.0 but are not stored
  block_inv[0] = 1.0;
  for ( i = 1; i < dof; i++ ) {
    block_inv[dof*i+i] = 1.0;
    for ( j = 0; j < i; j++ ) {
      // dot product of i'th row of block with j'th col of inv is 0
      PetscScalar dot = -coeff[idx0 + j*stripeSize + i];
      for ( k = j+1; k < i; k++ )
 	dot -= coeff[idx0 + k*stripeSize + i] * block_inv[j*dof + k];
      block_inv[j*dof + i] = dot;
    }
  }
  
  for ( j = 0; j < dof-1; j++ )
    for ( i = j+1; i < dof; i++ )
      coeff[idx0 + j*stripeSize + i] = block_inv[j*dof + i];

  // invert upper triangular matrix in block
  block_inv[dof*dof-1] = 1.0 / coeff[idx0 + (dof-1)*stripeSize + (dof-1)];
  for ( i = dof-2; i >= 0; i-- ) {
    block_inv[dof*i+i] = 1.0 / coeff[idx0 + i*stripeSize + i];
    for ( j = i+1; j < dof; j++ ) {
      // dot product of i'th row of block with j'th col of inv is 0
      PetscScalar dot = 0;
      for ( k = i+1; k <= j; k++ )
 	dot -= coeff[idx0 + k*stripeSize + i] * block_inv[j*dof + k];
      block_inv[j*dof + i] = dot / coeff[idx0 + i*stripeSize + i];
    }
  }
  
  for ( j = 0; j < dof; j++ )
    for ( i = 0; i <=j; i++ )
      coeff[idx0 + j*stripeSize + i] = block_inv[j*dof + i];
}

//---
//
// This kernel is invoked with input parameters
//     chunk     -- indicates the current portion of the vector of unknowns being solved for
//     chunkSize -- indicates how many unknowns are solved by one invocation of this kernel
//     Note that the unit of a chunk is a block of size dof.
//
// The unknowns already solved for in previous invocations have indices in the range [ 0 , ... , chunk*chunkSize*dof )
// Unknowns being solved for in this invocation have indices in the range [ chunk*chunkSize*dof , ... , (chunk+1)*chunkSize*dof )
// A given thread is solving for the unknown with index idx = chunk*chunkSize*dof + blockDim.x * blockIdx.x + threadIdx.x
//
// The diagonals that come into play for this thread are determined as follows:
// First note that the vector of unknowns is divided into blocks of size dof.
// Let blockIdx = idx / dof;  -- this is the index within the vector of unknowns of
//                               the dof-size block the working thread's unknown belongs to.
// Suppose a diagonal has offset -d  (this is a lower triangular solve so all diagonals have negative offset)
//
// Then the diagonal contributes to the equation involving the working thread's unknown in two possible ways:
// 1. It contains a block with block indices [blockIdx,col] where 0 <= col < chunk*chunkSize. This means it
//    contributes to the current equation with factors involving unknowns already solved for.
//    Since col = blockIdx - d, to determine if a diagonal falls into this category we need to determine if
//                                  d <= blockIdx < d + chunk*chunkSize
//                                                OR
//                                  blockIdx - chunk*chunkSize < d <= blockIdx
//
// 2. It contains a block with block indices [blockIdx,col] where  chunk*chunkSize <= col < (chunk+1)*chunkSize.
//    This means it contributes to the current equation with factors involving unknows that are being solved
//    for concurrently. To fall into this category we must have
//                                  chunk*chunkSize + d <= blockIdx < (chunk+1)*chunkSize + d
//                                                OR
//                                  d <= blockIdx - chunk*chunkSize
//---
__global__ void MatSolveLowerKernel(PetscScalar * coeff, PetscScalar * y, PetscInt mat_size,
				    int * diagonals, PetscInt center_diag,
				    PetscInt chunk, PetscInt chunkSize, PetscInt dof) {
  int idx = blockDim.x * blockIdx.x * 1 + threadIdx.x * 1;

  if (idx >= mat_size)
    return;

  int locidx0 = idx;
  int idx0 = locidx0 + chunk*chunkSize*dof;
  if (idx0 >= mat_size)
    return;

  int blockIdx = idx0 / dof;
  int offset0;
  int block0; 
  __shared__ double beta[BETAMAX];
  __shared__ double temp[BETAMAX];

  // Initialize the solution to corresponding RHS value
  PetscScalar rhs_val0 = fetch_double(vector_rhs, idx0);
  beta[locidx0] = rhs_val0;
  // do we need to sync threads here?

  int diag_size = mat_size * dof;

  // determine the diagonals that contribute already solved values to the RHS for this thread
  PetscInt min_solved_diag = 0, max_solved_diag = center_diag;
  while ( (min_solved_diag < max_solved_diag) && (diagonals[min_solved_diag] < -blockIdx) )
    ++min_solved_diag;
  while ( (0 <= max_solved_diag) && (diagonals[max_solved_diag] >= chunk*chunkSize - blockIdx) )
    --max_solved_diag;

  // multiply matrix*solved and subtract from RHS
  for (int i = min_solved_diag; i <= max_solved_diag; ++i) {
    int d = diagonals[i];

    offset0 = diag_size * i + idx0;
    block0 = (idx0 / dof + d) * dof;

    for (int j = 0; j < dof; ++j) {
      PetscScalar aval0 = coeff[offset0 + mat_size*j];
      PetscScalar yval0 = y[block0 + j];
      beta[locidx0] -= aval0 * yval0;
    }
  }

  // solve the system in the (chunkSize*dof)X(chunkSize*dof) block on the diagonal
  // there is not much parallelism in this loop
  // only certain threads have work to do
  offset0 = diag_size * center_diag + idx0;
  block0 = locidx0 / dof;
  int min_unsolved_diag = 0;
  while ( (min_unsolved_diag < center_diag) && (diagonals[min_unsolved_diag] < chunk*chunkSize - blockIdx) )
    ++min_unsolved_diag;

  for ( int i = 0; i < chunkSize; i++ ) {

    // solve for the next dof unknowns -- dof*dof block on diagonal is stored inverted
    // threads involved: i*dof,...,(i+1)*dof-1
    // the lower triangular elements on the center diagonal have been inverted
    // and the 1's on the diagonal are not stored.
    if ( block0 == i ) {
      temp[locidx0] = beta[locidx0];
      for ( int j = 0; j < locidx0-(i*dof); j++ ) {
 	PetscScalar aval0 = coeff[offset0 + mat_size*j];
 	PetscScalar bval0 = beta[i*dof + j];
 	temp[locidx0] += aval0 * bval0;
      }
      beta[locidx0] = temp[locidx0];
    }
 
    // update rhs for each diagonal
    // threads involved: (i-d)*dof,...,(i-d+1)*dof-1
    // for every non-zero diagonal d in the (chunkSize*dof)X(chunkSize*dof) block on the diagonal
    for ( int di = min_unsolved_diag; di < center_diag; di++ ) {
      int d = diagonals[di];
      // is this thread in a diagonal block?
      if ( block0 == (i-d) ) {
   	PetscInt offset1 = diag_size * di + idx0;
   	for ( int j = 0; j < dof; j++ ) {
   	  PetscScalar aval1 = coeff[offset1 + mat_size*j];
   	  PetscScalar bval1 = beta[i*dof + j];
   	  beta[locidx0] -= aval1 * bval1;
   	}	
      }
    }
  }
  
  y[idx0] = beta[locidx0];
  //y[idx0] = max_solved_diag;
}

//---
//
// This kernel is invoked with input parameters
//     chunk     -- indicates the current portion of the vector of unknowns being solved for
//     chunkSize -- indicates how many unknowns are solved by one invocation of this kernel
//     Note that the unit of a chunk is a block of size dof.
//
// The unknowns already solved for in previous invocations have indices in the range [ 0 , ... , chunk*chunkSize*dof )
// Unknowns being solved for in this invocation have indices in the range [ chunk*chunkSize*dof , ... , (chunk+1)*chunkSize*dof )
// A given thread is solving for the unknown with index idx = chunk*chunkSize*dof + blockDim.x * blockIdx.x + threadIdx.x
//
// The diagonals that come into play for this thread are determined as follows:
// First note that the vector of unknowns is divided into blocks of size dof.
// Let blockIdx = idx / dof;  -- this is the index within the vector of unknowns of
//                               the dof-size block the working thread's unknown belongs to.
// Suppose a diagonal has offset d  (this is an upper triangular solve so all diagonals have positive offset)
//
// Then the diagonal contributes to the equation involving the working thread's unknown in two possible ways:
// 1. It contains a block with block indices [blockIdx,col] where (chunk+1)*chunkSize <= col < mat_size/dof. This means it
//    contributes to the current equation with factors involving unknowns already solved for.
//    Since col = blockIdx + d, to determine if a diagonal falls into this category we need to determine if
//
//                            (chunk+1)*chunkSize - blockIdx <= d < mat_size/dof - blockIdx
//
// 2. It contains a block with block indices [blockIdx,col] where  chunk*chunkSize <= col < (chunk+1)*chunkSize.
//    This means it contributes to the current equation with factors involving unknows that are being solved
//    for concurrently. To fall into this category we must have
//                                  chunk*chunkSize - d <= blockIdx < (chunk+1)*chunkSize - d
//                                                OR
//                                  d < (chunk+1)*chunkSize - blockIdx
//---
__global__ void MatSolveUpperKernel(PetscScalar * coeff, PetscScalar * y, // PetscScalar * rhs,
				    PetscInt grid_size,
				    int * diagonals, PetscInt center_diag,
				    PetscInt chunk, PetscInt chunkSize, PetscInt dof) {
  int idx = blockDim.x * blockIdx.x * 1 + threadIdx.x * 1;

  PetscInt mat_size = grid_size * dof;

  if (idx >= mat_size)
    return;

  int locidx0 = idx;
  int idx0 = locidx0 + chunk*chunkSize*dof;
  if (idx0 >= mat_size)
    return;

  int blockIdx = idx0 / dof;
  int offset0;
  int block0; 
  __shared__ double beta[BETAMAX];
  __shared__ double temp[BETAMAX];

  // Initialize the solution to corresponding RHS value
  PetscScalar rhs_val0 = y[idx0]; //fetch_double(vector_rhs, idx0);
  beta[locidx0] = rhs_val0;
  // do we need to sync threads here?

  int diag_size = mat_size * dof;

  // determine the diagonals that contribute already solved values to the RHS for this thread
  PetscInt min_solved_diag = center_diag, max_solved_diag = 2*center_diag;
  while ( (min_solved_diag <= max_solved_diag) && (diagonals[min_solved_diag] < (chunk+1)*chunkSize - blockIdx) )
    ++min_solved_diag;
  while ( (center_diag <= max_solved_diag) && (diagonals[max_solved_diag] >= grid_size - blockIdx) )
    --max_solved_diag;

  // multiply matrix*solved and subtract from RHS
  for (int i = min_solved_diag; i <= max_solved_diag; ++i) {
    int d = diagonals[i];

    offset0 = diag_size * i + idx0;
    block0 = (idx0 / dof + d) * dof;

    for (int j = 0; j < dof; ++j) {
      PetscScalar aval0 = coeff[offset0 + mat_size*j];
      PetscScalar yval0 = y[block0 + j];
      beta[locidx0] -= aval0 * yval0;
    }
  }

  // solve the system in the (chunkSize*dof)X(chunkSize*dof) block on the diagonal
  // there is not much parallelism in this loop
  // only certain threads have work to do
  offset0 = diag_size * center_diag + idx0;
  block0 = locidx0 / dof;
  int max_unsolved_diag = 2*center_diag;
  while ( (max_unsolved_diag > center_diag) && (diagonals[max_unsolved_diag] >= (chunk+1)*chunkSize - blockIdx) )
    --max_unsolved_diag;

  for ( int i = chunkSize-1; i >= 0; i-- ) {

    // solve for the next dof unknowns -- dof*dof block on diagonal is stored inverted
    // threads involved: i*dof,...,(i+1)*dof-1
    if ( block0 == i ) {
      temp[locidx0] = 0;
      for ( int j = locidx0-(i*dof); j < dof; j++ ) {
 	PetscScalar aval0 = coeff[offset0 + mat_size*j];
 	PetscScalar bval0 = beta[i*dof + j];
 	temp[locidx0] += aval0 * bval0;
      }
      beta[locidx0] = temp[locidx0];
    }
 
    // update rhs for each diagonal
    // threads involved: (i-d)*dof,...,(i-d+1)*dof-1
    // for every non-zero diagonal d in the (chunkSize*dof)X(chunkSize*dof) block on the diagonal
    for ( int di = max_unsolved_diag; di > center_diag; di-- ) {
      int d = diagonals[di];
      // is this thread in a diagonal block?
      if ( block0 == (i-d) ) {
   	PetscInt offset1 = diag_size * di + idx0;
   	for ( int j = 0; j < dof; j++ ) {
   	  PetscScalar aval1 = coeff[offset1 + mat_size*j];
   	  PetscScalar bval1 = beta[i*dof + j];
   	  beta[locidx0] -= aval1 * bval1;
   	}	
      }
    }
  }
  
  // y[idx0] = rhs_val0;
  y[idx0] = beta[locidx0];
  // y[idx0] = max_unsolved_diag;
}


//===-- Host Code --------------------------------------------------------===//

static PetscErrorCode InvertFactoredDiagBlocks( Mat_SeqSGGPU *a, PetscScalar *data )
{
  PetscErrorCode ierr = 0;
  PetscInt dof = a->dof, numBlocks = a->m * a->n * a->p;
  PetscInt numElements = numBlocks * dof;
  PetscInt i,j,k, offset;
  PetscInt          num_diags = a->diagonals->size();

  PetscScalar block_inv[MAXBLOCKSIZE];

  for ( int row = 0; row < numBlocks; row++ ) {

    offset = (num_diags/2) * numElements * dof + row*dof;

    // invert lower triangular matrix in block
    // diagonal elements are all 1.0 but are not stored
    block_inv[0] = 1.0;
    for ( i = 1; i < dof; i++ ) {
      block_inv[dof*i+i] = 1.0;
      for ( j = 0; j < i; j++ ) {
	// dot product of i'th row of block with j'th col of inv is 0
	PetscScalar dot = -data[offset + j*numElements + i];
	for ( k = j+1; k < i; k++ )
	  dot -= data[offset + k*numElements + i] * block_inv[j*dof + k];
	block_inv[j*dof + i] = dot;
      }
    }
  
    for ( j = 0; j < dof-1; j++ )
      for ( i = j+1; i < dof; i++ )
	data[offset + j*numElements + i] = block_inv[j*dof + i];

    // invert upper triangular matrix in block
    block_inv[dof*dof-1] = 1.0 / data[offset + (dof-1)*numElements + (dof-1)];
    for ( i = dof-2; i >= 0; i-- ) {
      block_inv[dof*i+i] = 1.0 / data[offset + i*numElements + i];
      for ( j = i+1; j < dof; j++ ) {
	// dot product of i'th row of block with j'th col of inv is 0
	PetscScalar dot = 0;
	for ( k = i+1; k <= j; k++ )
	  dot -= data[offset + k*numElements + i] * block_inv[j*dof + k];
	block_inv[j*dof + i] = dot * block_inv[dof*i+i];
      }
    }
  
    for ( j = 0; j < dof; j++ )
      for ( i = 0; i <=j; i++ )
	data[offset + j*numElements + i] = block_inv[j*dof + i];
  }

  PetscFunctionReturn(ierr);
}

static PetscErrorCode MatGetBlockIJ_SeqSGGPU( Mat_SeqSGGPU * mat, PetscInt i, PetscInt j,
					      PetscScalar** blockPtr )
{
  PetscErrorCode ierr = 0;

  *blockPtr = PETSC_NULL;

  PetscInt numBlockRows = mat->m * mat->n * mat->p;
  if ( (i >= 0) && (i < numBlockRows) && (j >= 0) && (j < numBlockRows) ) {
    PetscInt blockIndex = -1;
    std::map<int, int> &diag_starts = *(mat->diag_starts);
    std::map<int, int>::iterator I = diag_starts.find(j-i);
    if (I != diag_starts.end()) {
      blockIndex = I->second + i*mat->dof;
      *blockPtr = &(mat->hostData[blockIndex]);
    }
  }

  PetscFunctionReturn(ierr);
}


//------------------------------------
// updateMainDiagonalBlock
// does LU decomp on a dofxdof block
// on the main diagonal
//
// Example for dof=3, topRow=0:
//
// INPUT:
// a  b  c   quotients = (*, d/a, g/a)
// d  e  f             = (*, q, r)
// g  h  i
//
// OUTPUT:
// a  b     c
// q  e-qb  f-qc
// r  h-rb  i-rc
//-------------------------------------
static PetscErrorCode MatUpdateMainDiagonalBlock_SeqSGGPU( PetscInt numElems, PetscInt dof,
							   PetscScalar *quotients, PetscScalar *block, int topRow )
{
  PetscErrorCode ierr = 0;
  PetscScalar quot;

  for ( int i = topRow+1; i < dof; i++ ) {
    quot = quotients[i];
    block[topRow*numElems + i] = quot;
    for ( int j = topRow+1; j < dof; j++ )
      block[j*numElems + i] -= quot*block[j*numElems+topRow];
  }

  PetscFunctionReturn(ierr);
}

//------------------------------------------
// MatUpdateTopRowBlock
// extends the LU decomp from a block on the
// main diagonal to a block in the same row
//
// Example for dof=3 and topRow=0:
//
// INPUT:
// a  b  c   quotients = (*, q, r)
// d  e  f
// g  h  i
//
// OUTPUT:
// a     b     c
// d-qa  e-qb  f-qc
// g-ra  h-rb  i-rc
//-------------------------------------
static PetscErrorCode MatUpdateTopRowBlock_SeqSGGPU( PetscInt numElems, PetscInt dof,
						     PetscScalar *quotients, PetscScalar *block, int topRow )
{
  PetscErrorCode ierr = 0;
  PetscScalar quot;

  for ( int i = topRow+1; i < dof; i++ ) {
    quot = quotients[i];
    for ( int j = 0; j < dof; j++ )
      block[j*numElems + i] -= quot*block[j*numElems + topRow];
  }

  return ierr;
}

//--------------------------------------------
// MatUpdateLeftColBlock
// extends the LU decomp from a block on the
// main diagonal to a block in the same column
//
// Example for dof=3 and startCol=0:
//
// INPUT:
// topRowVals = (x, y, z)
//
// a  b  c   quotients = (a/x, d/x, g/x)
// d  e  f             = (p, q, r)
// g  h  i
//
// OUTPUT:
// p  b-py  c-pz
// q  e-qy  f-qz
// r  h-ry  i-rz
//-------------------------------------
static PetscErrorCode MatUpdateLeftColBlock_SeqSGGPU( PetscInt numElems, PetscInt dof,
						      PetscScalar *quotients, PetscScalar *topRowVals,
						      PetscScalar *block, int startCol )
{
  PetscErrorCode ierr = 0;

  PetscScalar quot;

  for ( int i = 0; i < dof; i++ ) {
    quot = quotients[i];
    block[startCol*numElems + i] = quot;
    for ( int j = startCol+1; j < dof; j++ )
      block[j*numElems + i] -= quot*topRowVals[j];
  }

  PetscFunctionReturn(ierr);
}

//--------------------------------------------
// MatUpdateGeneralBlock
// extends the LU decomp from a block on the
// main diagonal to a non-zero block below and
// to the right
//
// Example:
//
// INPUT:
// topRowVals    = (x, y, z)
//
// a  b  c   quotients = (p, q, r)
// d  e  f             
// g  h  i
//
// OUTPUT:
// a-px  b-py  c-pz
// b-qx  e-qy  f-qz
// c-rx  h-ry  i-rz
//-------------------------------------
static PetscErrorCode MatUpdateGeneralBlock_SeqSGGPU( PetscInt numElems, PetscInt dof,
						      PetscScalar *quotients, PetscScalar *topRowVals,
						      PetscScalar *block )
{
  PetscErrorCode ierr = 0;

  PetscScalar quot;

  for ( int i = 0; i < dof; i++ ) {
    quot = quotients[i];
    for ( int j = 0; j < dof; j++ )
      block[j*numElems + i] -= quot*topRowVals[j];
  }

  PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "MatLUFactorNumeric_SeqSGGPU"
PetscErrorCode MatLUFactorNumeric_SeqSGGPU(Mat B,Mat A,const MatFactorInfo *info)
{
  Mat              C=B;
  Mat_SeqSGGPU       *a=(Mat_SeqSGGPU*)A->data,*b=(Mat_SeqSGGPU *)C->data;
  PetscInt dof = a->dof, numBlocks = a->m * a->n * a->p;
  PetscInt numElements = numBlocks * dof;
  PetscInt diagSize = numElements * dof;
  PetscScalar *topRowValsDiag, *topRowValsOffDiag, *quotients;
  PetscErrorCode ierr;
  PetscScalar *mainDiagonalBlock, *topRowBlock, *leftColBlock, *generalBlock;
#if USE_ROWBLOCKS
  PetscScalar *blockedHostData, *pBHD;
  PetscInt coord;
#endif
  PetscInt i,j,k, offset;
  PetscInt          num_diags = b->diagonals->size();
  PetscInt          size = diagSize * (num_diags+1);
  PetscInt mainDiagStartIndex = numElements * dof * (num_diags/2);
  std::vector<int>::iterator J, I;

  PetscFunctionBegin;

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif

  // Copy the diagonals from A to B
  int diag_offset = 0;
  for ( J = b->diagonals->begin();
	J != b->diagonals->end(); J++ ) {
    // first store the start index of this diagonal in b->diag_starts
    int d = *J;

    std::map<int, int>::iterator I = a->diag_starts->find(d);
    // if this diag is in A, copy the data into B
    if (I != a->diag_starts->end()) {
      ierr = PetscMemcpy( &(b->hostData[diag_offset]), &(a->hostData[I->second]),
			  diagSize*sizeof(PetscScalar)); CHKERRQ(ierr);
    }
    // otherwise, zero it out
    else {
      ierr = PetscMemzero(&(b->hostData[diag_offset]),
			  diagSize*sizeof(PetscScalar));CHKERRQ(ierr);
    }
    
    diag_offset += diagSize;
    if ( d == 0 )
      diag_offset += diagSize;
  }

  if (b->deviceData) {
    cudaFree(b->deviceData);
  }
  checkCudaError(cudaMalloc(&b->deviceData, sizeof(PetscScalar) * size));
  //checkCudaError(cudaMemset(b->deviceData,0.0,sizeof(PetscScalar)*size));

  if (b->deviceDiags) {
    cudaFree(b->deviceDiags);
  }
  checkCudaError(cudaMalloc(&b->deviceDiags, sizeof(int) * b->diagonals->size()));

  checkCudaError(cudaMemcpyAsync(b->deviceDiags, &(*b->diagonals)[0], sizeof(int) * b->diagonals->size(), cudaMemcpyHostToDevice, b->stream));
  checkCudaError(cudaMemcpy(b->deviceData, b->hostData, sizeof(PetscScalar) * size, cudaMemcpyHostToDevice));
  cudaDeviceSynchronize();


  ierr = PetscMalloc(dof * sizeof(PetscScalar), &topRowValsDiag); CHKERRQ(ierr);
  ierr = PetscMalloc(dof * sizeof(PetscScalar), &topRowValsOffDiag); CHKERRQ(ierr);
  ierr = PetscMalloc(dof * sizeof(PetscScalar), &quotients); CHKERRQ(ierr);

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("factor numeric preamble time %lf\n",elapsed);
  t_start = getclock();
#endif

  // iterate over blocks in main diagonal
  for ( int row = 0; row < numBlocks; row++ ) {
    
    //ierr = MatGetBlockIJ_SeqSGGPU( b, row, row, &mainDiagonalBlock );
    mainDiagonalBlock = &(b->hostData[mainDiagStartIndex + row*dof]);

    // we need to factor block(row,row) and extend this to
    //    -- any nonzero block(row,row+d) on a super-diagonal
    //    -- any nonzero block(row+d,row) on a subdiagonal
    //    -- any nonzero block(row+d1,row+d2)
    //
    // the outer loop here is over the rows of all blocks
    // that need to be updated
    for ( int k = 0; k < dof; k++ ) {

      // copy the values in the current row of the main diagonal block(row,row)
      // these will be used to update any nonzero block(row+d,row)
      for ( int j = k; j < dof; j++ ) {
	topRowValsDiag[j] = mainDiagonalBlock[j*numElements + k];
      }

      quotients[k] = 1.0;
      for ( int j = k+1; j < dof; j++ )
	quotients[j] = mainDiagonalBlock[k*numElements +j] / mainDiagonalBlock[k*numElements + k];

      // update block(row,row)[k+1..dof-1][k+1..dof-1] using the quotients from column k
      ierr = MatUpdateMainDiagonalBlock_SeqSGGPU( numElements, dof, quotients, mainDiagonalBlock, k );

      // update the blocks to the right using the same quotients from the k'th column of block(row,row)
      for (int i = 0; i < num_diags; ++i) {
	int d = (*b->diagonals)[i];
	if ( (d > 0) && (row+d < numBlocks) ) {
	  //ierr = MatGetBlockIJ_SeqSGGPU( b, row, row+d, &topRowBlock );
	  topRowBlock = &(b->hostData[(i+1)*dof*numElements + row*dof]);

	  if ( topRowBlock )
	    ierr = MatUpdateTopRowBlock_SeqSGGPU( numElements, dof, quotients, topRowBlock, k );
	}
      }

      for (int i = 0; i < num_diags; ++i) {
	int dSub = (*b->diagonals)[i];
	if ( (dSub < 0) && (row-dSub < numBlocks) ) {
	  //ierr = MatGetBlockIJ_SeqSGGPU( b, row-dSub, row, &leftColBlock );
	  leftColBlock = &(b->hostData[i*dof*numElements + (row-dSub)*dof]);
	  if ( leftColBlock != 0 ) {
	    for ( int j = 0; j < dof; j++ )
	      quotients[j] = leftColBlock[k*numElements +j] / mainDiagonalBlock[k*numElements + k];

	    ierr = MatUpdateLeftColBlock_SeqSGGPU( numElements, dof, quotients, topRowValsDiag, leftColBlock, k );

	    for (int ii = 0; ii < num_diags; ++ii) {
	      int dSuper = (*b->diagonals)[ii];
	      if ( (dSuper > 0) && (row+dSuper < numBlocks) ) {
		//ierr = MatGetBlockIJ_SeqSGGPU( b, row, row+dSuper, &topRowBlock );
		topRowBlock = &(b->hostData[(ii+1)*dof*numElements + row*dof]);

		if ( topRowBlock != 0 ) {
		  ierr = MatGetBlockIJ_SeqSGGPU( b, row-dSub, row+dSuper, &generalBlock );
		  if ( generalBlock != 0 ) {
		    for ( int j = 0; j < dof; j++ )
		      topRowValsOffDiag[j] = topRowBlock[j*numElements + k];
		    
		    ierr = MatUpdateGeneralBlock_SeqSGGPU( numElements, dof, quotients, topRowValsOffDiag, generalBlock);
		  }
		}
	      }
	    }
	  }
	}
      }
    }      
  }

  // Now split the center diagonal
  PetscScalar *lowerMainDiag, *upperMainDiag;
  //ierr = MatGetBlockIJ_SeqSGGPU( b, 0, 0, &lowerMainDiag );
  lowerMainDiag = &(b->hostData[mainDiagStartIndex]);
  upperMainDiag = lowerMainDiag + diagSize;
  ierr = PetscMemzero(upperMainDiag,diagSize*sizeof(PetscScalar));CHKERRQ(ierr);

  for ( int row = 0; row < numBlocks; row++ ) {
    for ( int j = 0; j < dof; j++ )
      for ( int k = 0; k <= j; k++ ) {
	upperMainDiag[j*numElements+k] = lowerMainDiag[j*numElements+k];
	if ( j == k )
	  lowerMainDiag[j*numElements+k] = 1.0;
	else
	  lowerMainDiag[j*numElements+k] = 0.0;
      }
    lowerMainDiag += dof;
    upperMainDiag += dof;
  }

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("factor numeric main loop time %lf\n",elapsed);
  t_start = getclock();
#endif

  C->ops->solve = MatSolve_SeqSGGPU;

  //=================================
  //invert the blocks on the diagonal
  //=================================
  PetscScalar block_inv[MAXBLOCKSIZE];

  for ( int row = 0; row < numBlocks; row++ ) {

    offset = (b->diagonals->size()/2) * numElements * dof + row*dof;

    // invert lower triangular matrix in block
    // diagonal elements are all 1.0 but are not stored
    block_inv[0] = 1.0;
    for ( i = 1; i < dof; i++ ) {
      block_inv[dof*i+i] = 1.0;
      for ( j = 0; j < i; j++ ) {
	// dot product of i'th row of block with j'th col of inv is 0
	PetscScalar dot = -b->hostData[offset + j*numElements + i];
	for ( k = j+1; k < i; k++ )
	  dot -= b->hostData[offset + k*numElements + i] * block_inv[j*dof + k];
	block_inv[j*dof + i] = dot;
      }
    }
  
    for ( j = 0; j < dof-1; j++ )
      for ( i = j+1; i < dof; i++ )
	b->hostData[offset + j*numElements + i] = block_inv[j*dof + i];

    // invert upper triangular matrix in block
    offset += numElements*dof;
    block_inv[dof*dof-1] = 1.0 / b->hostData[offset + (dof-1)*numElements + (dof-1)];
    for ( i = dof-2; i >= 0; i-- ) {
      block_inv[dof*i+i] = 1.0 / b->hostData[offset + i*numElements + i];
      for ( j = i+1; j < dof; j++ ) {
	// dot product of i'th row of block with j'th col of inv is 0
	PetscScalar dot = 0;
	for ( k = i+1; k <= j; k++ )
	  dot -= b->hostData[offset + k*numElements + i] * block_inv[j*dof + k];
	block_inv[j*dof + i] = dot * block_inv[dof*i+i];
      }
    }
  
    for ( j = 0; j < dof; j++ )
      for ( i = 0; i <=j; i++ )
	b->hostData[offset + j*numElements + i] = block_inv[j*dof + i];
  }


#if USE_ROWBLOCKS
  // reformat the matrix so that a block of ROWBLOCK rows is stored consecutively
  ierr = PetscMalloc( size * sizeof(PetscScalar), &blockedHostData); CHKERRQ(ierr);
  //  ierr = PetscMemcpy( blockedHostData, b->hostData,
  //	      size*sizeof(PetscScalar)); CHKERRQ(ierr);

  pBHD = blockedHostData;
  // FIRST STORE THE LOWER TRIANGULAR MATRIX
  for ( i = 0; i < numBlocks; i++ ) {
    for ( j = 0; j < (num_diags/2)+1; j++ ) {
      for ( k = 0; k < dof; k++ ) {
	offset = i*dof + j*diagSize + k*numElements;
	for ( coord = 0; coord < dof; coord++ )
	  *(pBHD++) = b->hostData[offset+coord];
      }
    }
  }
  // NOW STORE THE UPPER TRIANGULAR MATRIX
  for ( i = 0; i < numBlocks; i++ ) {
    for ( j = (num_diags/2)+1; j < num_diags+1; j++ ) {
      for ( k = 0; k < dof; k++ ) {
	offset = i*dof + j*diagSize + k*numElements;
	for ( coord = 0; coord < dof; coord++ )
	  *(pBHD++) = b->hostData[offset+coord];
      }
    }
  }


  PetscFree(b->hostData);
  b->hostData = blockedHostData;
#endif

  // this will copy hostData to deviceData
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  int mat_size = b->m * b->n * b->p;
  PetscInt centerDiag = ( b->diagonals->size() / 2 ) * mat_size * b->dof * b->dof;


#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("factor numeric finalize time %lf\n",elapsed);
#endif

  PetscFunctionReturn(0); 
}

PetscErrorCode MatILUFactor_SeqSGGPU(Mat inA,IS row,IS col,const MatFactorInfo *info)
{
  PetscErrorCode ierr;
  Mat            outA = inA;

  ierr = MatLUFactorNumeric_SeqSGGPU(outA,inA,info);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

extern PetscErrorCode MatDuplicateNoCreate_SeqSGGPU(Mat,Mat,MatDuplicateOption,PetscBool );
extern PetscErrorCode MatGetRow_SeqSGGPU(Mat A, PetscInt row, PetscInt * nz, PetscInt **idx , PetscScalar ** v);
extern PetscErrorCode MatRestoreRow_SeqSGGPU(Mat A, PetscInt row, PetscInt *nz, PetscInt **idx, PetscScalar **v);
extern PetscErrorCode MatSetStencil_SeqSGGPU(Mat A, PetscInt dim, const PetscInt dims[], const PetscInt starts[], PetscInt dof);

#undef __FUNCT__  
#define __FUNCT__ "MatILUFactorSymbolic_SeqSGGPU"
PetscErrorCode MatILUFactorSymbolic_SeqSGGPU(Mat fact,Mat A,IS isrow,IS iscol,const MatFactorInfo *info)
{
  //Mat_SeqSGGPU         *a = (Mat_SeqSGGPU*)A->data; //,*b;
  IS                    isicol;
  PetscErrorCode        ierr;
  PetscInt              iluLevel = (PetscInt)info->levels;
  Mat_SeqSGGPU         *a = (Mat_SeqSGGPU*)A->data, *b;
  std::vector<int>     *iluDiagonals, newDiagonals;
  int                   j, nz, numDiag;
  std::vector<int>::iterator J, I, II;

  PetscFunctionBegin;
  if (A->rmap->n != A->cmap->n) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Must be square matrix, rows %D columns %D",A->rmap->n,A->cmap->n);

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif

  ierr = ISInvertPermutation(iscol,PETSC_DECIDE,&isicol);CHKERRQ(ierr);

  // duplicate the diagonals of the input matrix,
  // then add new diagonals for each level of ilu(k)
  // a new diagonal at offset (d1+d2) is added for
  // every d1>0 and d2<0 in the current set of diagonals.
  // maybe keep some extra information around about the
  // elements of the new diagonals which are known to be 0?
  iluDiagonals = new std::vector<int>( *(a->diagonals) );
  numDiag = a->diagonals->size();

  for ( j = 0; j < iluLevel; j++ ) {
    
    for ( I = iluDiagonals->begin();
	  I != iluDiagonals->end(); I++ ) {
      int d1 = *I;
      if ( d1 < 0 ) {
	for ( J = iluDiagonals->begin();
	      J != iluDiagonals->end(); J++ ) {
	  int d2 = *J;
	  if ( d2 > 0 ) {
	    II = find( iluDiagonals->begin(), iluDiagonals->end(), d1+d2 );
	    if ( II == iluDiagonals->end() ) {
	      II = find( newDiagonals.begin(), newDiagonals.end(), d1+d2 );
	      if ( II == newDiagonals.end() ) {
		newDiagonals.push_back(d1+d2);
	      }
	    }
	  }
	}
      }
    }

    numDiag += newDiagonals.size();
    iluDiagonals->insert( iluDiagonals->end(), newDiagonals.begin(), newDiagonals.end() );
    newDiagonals.clear();

#if _TRACE
    printf("Level %d diagonals:  ", j+1);
    for ( J = iluDiagonals->begin();
	  J != iluDiagonals->end(); J++ )
      printf(" %d ", *J);
    printf("\n");
#endif

  }

  sort( iluDiagonals->begin(), iluDiagonals->end() );
#if _TRACE
  printf("Sorted level %d diagonals:  ", iluLevel);
  for ( J = iluDiagonals->begin();
	J != iluDiagonals->end(); J++ )
    printf(" %d ", *J);
  printf("\n");
#endif

  // Duplicate the input matrix A so that the hostData of fact will be the same
  nz = a->m * a->n * a->p * a->dof;
  ierr = MatSetSizes(fact,nz,nz,nz,nz);CHKERRQ(ierr);
  MatSetType(fact,MATSEQSGGPU);
  b = (Mat_SeqSGGPU*)fact->data;
  int diagSize = nz * a->dof;
  int num_diags = iluDiagonals->size();
  b->stpoints = num_diags;
  b->dim = a->dim;
  b->dof = a->dof;
  b->m = a->m;
  b->n = a->n;
  b->p = a->p;
  
  // The center diagonal (offset=0) is twice as long as the others because
  // it stores the diagonal blocks for both L and U.
  // We want to store these separately because cuda kernels need to have
  // as little logic in them as possible.
  int size = diagSize * (num_diags+1);
  int diag_offset = 0;
  ierr = PetscMalloc( size * sizeof(PetscScalar), &b->hostData ); CHKERRQ(ierr);

  b->diagonals = iluDiagonals;

  fact->factortype             = MAT_FACTOR_ILU;
  fact->info.factor_mallocs    = 0;
  fact->info.fill_ratio_given  = info->fill;
  fact->info.fill_ratio_needed = 1.0;
  fact->ops->lufactornumeric   = MatLUFactorNumeric_SeqSGGPU;

  for ( J = b->diagonals->begin();
	J != b->diagonals->end(); J++ ) {
    // store the start index of this diagonal in b->diag_starts
    int d = *J;
    (*b->diag_starts)[d] = diag_offset;
    diag_offset += diagSize;
    if ( d == 0 )
      diag_offset += diagSize;
  }

  // ierr    = PetscMalloc((fact->rmap->n+1)*sizeof(PetscScalar),&b->solve_work);CHKERRQ(ierr);
  ierr    = PetscObjectReference((PetscObject)isrow);CHKERRQ(ierr);
  ierr    = PetscObjectReference((PetscObject)iscol);CHKERRQ(ierr);

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("factor symbolic time %lf\n",elapsed);
#endif

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatSolve_SeqSGGPU"
//---
//
// Choice Version
//
//---
PetscErrorCode MatSolve_SeqSGGPU(Mat A,Vec bb,Vec xx)
{
  Mat_SeqSGGPU        *a = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode ierr;

  if ( a->dof == 1 ) {
#if USE_ROWBLOCKS
    ierr = MatSolve_SeqSGGPU_cpu_1_nochunk(A,bb,xx);
#else
    ierr = MatSolve_SeqSGGPU_cpu_1(A,bb,xx);
#endif
    PetscFunctionReturn(0);
  }
  else if (a->dof != 4 ) {
    ierr = MatSolve_SeqSGGPU_cpu_simple(A,bb,xx);
    PetscFunctionReturn(0);
  }
  else {
    PetscBool          isseqgpu;
    PetscScalar       *x;
    const PetscScalar *b;
    int diagOffsets[MAXNUMDIAGS];
    

    if ( dumpVec ) {
      char  sgvecFname[64];
      dumpVec = PETSC_FALSE;
      PetscViewer rhsViewer;
      sprintf( sgvecFname, "sgvec-ex%d-%d-%d-%d-%d.bin", exNumber, a->m, a->n, a->p, a->dof );
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sgvecFname,FILE_MODE_WRITE, &rhsViewer);CHKERRQ(ierr);
      ierr = VecView(bb,rhsViewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&rhsViewer);CHKERRQ(ierr);
    }

    ierr = PetscObjectTypeCompare((PetscObject)xx,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
    if (isseqgpu) {
      Vec_SeqGPU *xData = (Vec_SeqGPU*) xx->data;
      if (xData->syncState==VEC_GPU) {
	ierr = VecCopyOverD2H(xx,xData->cpuptr);CHKERRQ(ierr);
	xData->syncState=VEC_SYNCHED;
      }
    }
    ierr = PetscObjectTypeCompare((PetscObject)bb,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
    if (isseqgpu) {
      Vec_SeqGPU *bData = (Vec_SeqGPU*) bb->data;
      if (bData->syncState==VEC_GPU) {
	ierr = VecCopyOverD2H(bb,bData->cpuptr);CHKERRQ(ierr);
	bData->syncState=VEC_SYNCHED;
      }
    }

    ierr = VecGetArray(xx,&x);CHKERRQ(ierr);
    ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);

    int di = 0;
    for (std::vector<int>::iterator it = a->diagonals->begin() ; it != a->diagonals->end(); ++it)
      diagOffsets[di++] = *it;
    //printf("calling sse4 code\n");
#if USE_ROWBLOCKS
    ierr = MatSolve_SeqSGGPU_sse4_nochunk(a->m,a->n,a->p,a->dof,a->dim,a->hostData,
					  a->diagonals->size(),
					  diagOffsets, x, b );CHKERRQ(ierr);
#else
    ierr = MatSolve_SeqSGGPU_sse4(a->m,a->n,a->p,a->dof,a->dim,a->hostData,
				  a->diagonals->size(),
				  diagOffsets, x, b );CHKERRQ(ierr);
#endif

    ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
}

#undef __FUNCT__  
#define __FUNCT__ "MatSolve_SeqSGGPU_gpu"
//---
//
// GPU Version
//
//---
PetscErrorCode MatSolve_SeqSGGPU_gpu(Mat A,Vec bb,Vec xx)
{
  Mat_SeqSGGPU      *a = (Mat_SeqSGGPU*)A->data;
  PetscBool         isseqcusp,isseqgpu,ismpicusp,iscusp;
  PetscErrorCode    ierr;
  PetscInt          n=A->rmap->n, dof = a->dof;
  PetscInt          numDiags = a->diagonals->size();
  PetscInt          matSize, chunkSize, chunk, numChunks;
  PetscScalar       *devX, *devB;

  PetscFunctionBegin;
  if (!n) PetscFunctionReturn(0);

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif

  if ( dumpVec ) {
    dumpVec = PETSC_FALSE;
    PetscViewer rhsViewer;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"sgvec.full.bin",FILE_MODE_WRITE, &rhsViewer);CHKERRQ(ierr);
    ierr = VecView(bb,rhsViewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&rhsViewer);CHKERRQ(ierr);
  }

  // set chunkSize -- this determines how many rows of the solve to do as a unit
  chunkSize = a->m; //dim==3 ? a->p : dim==2 ? a->n : 1;

  CUSPARRAY * rhsgpu;
  CUSPARRAY * xgpu;

  ierr = PetscObjectTypeCompare((PetscObject)xx,VECSEQCUSP,&isseqcusp);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)xx,VECMPICUSP,&ismpicusp);CHKERRQ(ierr);
  iscusp = (isseqcusp || ismpicusp) ? PETSC_TRUE : PETSC_FALSE;
  ierr = PetscObjectTypeCompare((PetscObject)xx,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
  if (isseqgpu) {
    dim3 block(BLOCKWIDTH_X, BLOCKWIDTH_Y);
    dim3 grid((int)ceil((float)(a->m * a->n * a->p * a->dof)/(float)BLOCKWIDTH_X / 1.0), 1);

    Vec_SeqGPU *xd = (Vec_SeqGPU*) xx->data;
    /* Make sure bb is also VECSEQGPU */
    ierr = PetscObjectTypeCompare((PetscObject)bb,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
    if (!isseqgpu) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and b must be same type");
    }
    Vec_SeqGPU *bd = (Vec_SeqGPU*) bb->data;
    /* synch up x */
    if (xd->syncState==VEC_CPU) {
      ierr = VecCopyOverH2D(xx,xd->cpuptr);CHKERRQ(ierr);
      xd->syncState=VEC_SYNCHED;
    }
    /* Get device pointer for X */
    devX = xd->devptr;
    devB = bd->devptr;
    /* Bind X to device texture */
    matSize = a->m * a->n * a->p * a->dof;
    
    checkCudaError(cudaBindTexture(0, vector_rhs, devB, matSize * sizeof(PetscScalar)));    
  }


  else if (iscusp) {
    dim3 block(BLOCKWIDTH_X, BLOCKWIDTH_Y);
    dim3 grid((int)ceil((float)(a->m * a->n * a->p * a->dof)/(float)BLOCKWIDTH_X / 1.0), 1);

    /* Make sure y is also VECCUSP */
    ierr = PetscObjectTypeCompare((PetscObject)bb,VECSEQCUSP,&isseqcusp);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)bb,VECMPICUSP,&ismpicusp);CHKERRQ(ierr);
    iscusp = (isseqcusp || ismpicusp) ? PETSC_TRUE : PETSC_FALSE;
    if (!iscusp) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and b must be same type");
    }
    matSize = a->m * a->n * a->p * a->dof;
    ierr = VecCUSPGetArrayWrite(xx, &xgpu); CHKERRQ(ierr);
    ierr = VecCUSPGetArrayRead(bb, &rhsgpu); CHKERRQ(ierr);
    devX = thrust::raw_pointer_cast(&(*xgpu)[0]);
    devB = thrust::raw_pointer_cast(&(*rhsgpu)[0]);

    /* Bind X to device texture */
    checkCudaError(cudaBindTexture(0, vector_rhs, devB, matSize * sizeof(PetscScalar)));
  }

  else {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Vec must be seqgpu or cusp type");
  }

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("gpu solve preamble time %lf\n",elapsed);
  t_start = getclock();
#endif

  dim3 block(chunkSize*dof, BLOCKWIDTH_Y);
  dim3 grid(1, 1);
  // dim3 grid((int)ceil((float)(matSize)/(float)BLOCKWIDTH_X / 1.0), 1);
  int shared_size = 2*chunkSize*dof*sizeof(PetscScalar);
  numChunks = ((matSize/dof) + (chunkSize-1)) / chunkSize;
#if _TRACE
  printf("matSize = %d, chunkSize = %d, numChunks = %d, dof = %d\n", matSize, chunkSize, numChunks, dof);
#endif

  for ( chunk = 0; chunk < numChunks; chunk++ ) {
#if _TRACE
    //    printf("do cuda lower triangular solve for chunk %d\n", chunk);
#endif
    MatSolveLowerKernel<<<grid, block, shared_size, a->stream>>>(a->deviceData, devX, matSize, a->deviceDiags,
								 a->diagonals->size() / 2, chunk, chunkSize, dof);
  }

#if _TIME
  checkCudaError(cudaStreamSynchronize(a->stream));
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("gpu solve lower time %lf\n",elapsed);
  t_start = getclock();
#endif

  ierr = WaitForGPU() ; CHKERRCUSP(ierr);
  cudaUnbindTexture(vector_rhs);
  ierr = VecCUSPRestoreArrayRead(bb, &rhsgpu); CHKERRQ(ierr);

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("gpu solve unbind time %lf\n",elapsed);
  t_start = getclock();
#endif

#if _TRACE
  printf("done cuda lower triangular solve\n");
  PetscScalar *gpuCheckY;
  ierr = PetscMalloc(matSize * sizeof(PetscScalar), &gpuCheckY); CHKERRQ(ierr);
  checkCudaError(cudaMemcpy(gpuCheckY, a->deviceY, matSize * sizeof(PetscScalar), cudaMemcpyDeviceToHost));
#endif

  // checkCudaError(cudaBindTexture(0, vector_rhs, a->deviceY, matSize * sizeof(PetscScalar)));

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("gpu solve bind 2 time %lf\n",elapsed);
  t_start = getclock();
#endif

  for ( chunk = numChunks-1; chunk >= 0; chunk-- ) {
#if _TRACE
    //    printf("do cuda upper triangular solve for chunk %d\n", chunk);
#endif
    MatSolveUpperKernel<<<grid, block, shared_size, a->stream>>>(a->deviceData, devX, matSize, a->deviceDiags,
								 a->diagonals->size() / 2, chunk, chunkSize, dof);
  }


#if _TIME
  checkCudaError(cudaStreamSynchronize(a->stream));
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("gpu solve upper time  %lf\n", elapsed);
  t_start = getclock();
#endif

#if _TRACE
  printf("done cuda upper triangular solve\n");
#endif
  cudaUnbindTexture(vector_rhs);
  ierr = VecCUSPRestoreArrayWrite(xx, &xgpu); CHKERRQ(ierr);

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("gpu solve vec restore time %lf\n",elapsed);
  t_start = getclock();
#endif

  //ierr = PetscLogFlops(2*a->nz - A->cmap->n);CHKERRQ(ierr);

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("gpu cleanup time %lf\n",elapsed);
#endif

  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "MatSolve_SeqSGGPU_cpu_1_nochunk"
//---
//
// CPU Version for dof==1
//
//---
PetscErrorCode MatSolve_SeqSGGPU_cpu_1_nochunk(Mat A,Vec bb,Vec xx)
{
  Mat_SeqSGGPU        *a = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode    ierr;
  PetscInt          n=A->rmap->n;
  PetscInt          numDiags = a->diagonals->size();
  PetscScalar       *x;
  const PetscScalar *b, *uTri;
  PetscInt          offset1;
  PetscScalar       matElt,vecElt,sum;
  PetscInt          gridSize = a->m * a->n * a->p;
  PetscInt          numElements = gridSize;
  PetscInt          startDiag,endDiag;
  PetscBool         isseqgpu;
  PetscInt          tempDiags[MAXNUMDIAGS];

  PetscFunctionBegin;
  if (!n) PetscFunctionReturn(0);

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif

  if ( dumpVec ) {
    dumpVec = PETSC_FALSE;
    char  sgvecFname[64];
    PetscViewer rhsViewer;
    sprintf( sgvecFname, "sgvec-ex%d-%d-%d-%d-1.bin",  exNumber,a->m, a->n, a->p );
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sgvecFname,FILE_MODE_WRITE, &rhsViewer);CHKERRQ(ierr);
    ierr = VecView(bb,rhsViewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&rhsViewer);CHKERRQ(ierr);
  }  

  for ( int di = 0; di < a->diagonals->size(); di++ )
    tempDiags[di] = (*a->diagonals)[di];

  ierr = PetscObjectTypeCompare((PetscObject)xx,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
  if (isseqgpu) {
    ierr = PetscObjectTypeCompare((PetscObject)bb,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
    if (!isseqgpu) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and b must be same type");
    }
    Vec_SeqGPU *bd = (Vec_SeqGPU*) bb->data;
    /* synch up b */
    if (bd->syncState==VEC_GPU) {
      ierr = VecCopyOverD2H(bb,bd->cpuptr);CHKERRQ(ierr);
      bd->syncState=VEC_SYNCHED;
    }
  }
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);

  // LOWER TRIANGULAR SOLVE
  startDiag = 0;
  endDiag = numDiags / 2;
  x[0] = b[0];
  for ( int block = 1; block < gridSize; block++ ) {
    sum = b[block];
    for ( int di = startDiag; di < endDiag; di++ ) {
      offset1 = block + tempDiags[di];
      if ( offset1 > -1 ) {
	matElt = TRI_MATRIX_ENTRY1(a->hostData,di,block,numDiags,numElements);
	vecElt = x[offset1];
	sum -= vecElt*matElt;
      }
    }
    x[block] = sum;
  }
 
  // UPPER TRIANGULAR SOLVE
  startDiag = endDiag + 1;
  endDiag = numDiags;
  uTri = a->hostData + startDiag * numElements;
  for ( int block = gridSize-1; block >= 0; block-- ) {
    sum = x[block];
    for ( int di = startDiag; di < endDiag; di++ ) {
      offset1 = block + tempDiags[di];
      if ( offset1 < gridSize ) {
	matElt = TRI_MATRIX_ENTRY1(uTri, di-startDiag+1, block, numDiags, numElements);
	vecElt = x[offset1];
	sum -= vecElt*matElt;
      }
    }
    x[block] = TRI_MATRIX_ENTRY1(uTri, 0, block, numDiags, numElements) * sum;
  }

  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("cpu solve time %lf\n",elapsed);
#endif
  //ierr = PetscLogFlops(2*a->nz - A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatSolve_SeqSGGPU_cpu_1"
//---
//
// CPU Version for dof==1
//
//---
PetscErrorCode MatSolve_SeqSGGPU_cpu_1(Mat A,Vec bb,Vec xx)
{
  Mat_SeqSGGPU        *a = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode    ierr;
  PetscInt          n=A->rmap->n;
  PetscInt          numDiags = a->diagonals->size();
  PetscScalar       *x;
  const PetscScalar *b;
  PetscScalar       sca1, sca2;
  PetscInt          gridSize = a->m * a->n * a->p;
  PetscInt          numElements = gridSize;
  PetscInt          chunkSize, chunk, numChunks;
  PetscBool         isseqgpu;
  PetscInt          tempDiags[MAXNUMDIAGS];

  PetscFunctionBegin;
  if (!n) PetscFunctionReturn(0);

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif


  if ( dumpVec ) {
    dumpVec = PETSC_FALSE;
    char  sgvecFname[64];
    PetscViewer rhsViewer;
    sprintf( sgvecFname, "sgvec-ex%d-%d-%d-%d-1.bin",  exNumber,a->m, a->n, a->p );
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sgvecFname,FILE_MODE_WRITE, &rhsViewer);CHKERRQ(ierr);
    ierr = VecView(bb,rhsViewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&rhsViewer);CHKERRQ(ierr);
  }  

  // set chunkSize -- this determines how many rows of the solve to do as a unit
  chunkSize = a->m; //dim==3 ? a->p : dim==2 ? a->n : 1;
  for ( int di = 0; di < a->diagonals->size(); di++ )
    tempDiags[di] = (*a->diagonals)[di];

  // x is the input RHS, temp is the result of lower triangular solve
  ierr = PetscObjectTypeCompare((PetscObject)xx,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
  if (isseqgpu) {
    ierr = PetscObjectTypeCompare((PetscObject)bb,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
    if (!isseqgpu) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and b must be same type");
    }
    Vec_SeqGPU *bd = (Vec_SeqGPU*) bb->data;
    /* synch up b */
    if (bd->syncState==VEC_GPU) {
      ierr = VecCopyOverD2H(bb,bd->cpuptr);CHKERRQ(ierr);
      bd->syncState=VEC_SYNCHED;
    }
  }
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);


  // LOWER TRIANGULAR SOLVE
  // determine the diagonals that contribute already solved values to the RHS for this thread
  PetscInt center_diag = numDiags / 2;
  PetscInt max_solved_diag = center_diag;

  // initialize the RHS with the components of b
  // PetscMemcpy( a->iluSolveVec, b, numElements*sizeof(PetscScalar) );    
  // PetscMemcpy( x, b, numElements*sizeof(PetscScalar) );

  // PetscScalar sum[MAXBLOCKSIZE];

  numChunks = (numElements + (chunkSize-1)) / chunkSize;
#if _TRACE
  printf("matSize = %d, chunkSize = %d, numChunks = %d, dof = 1\n", numElements, chunkSize, numChunks);
#endif

  // proceed a chunk of rows at a time
  // result of LT solve is stored in a->iluSolveVec
  while ( (0 <= max_solved_diag) && (tempDiags[max_solved_diag] > -chunkSize ) )
    --max_solved_diag;

  for ( chunk = 0; chunk < numChunks; chunk++ ) {
    // offset into portion of iluSolveVec being solved for this chunk
    int offset2 = chunk*chunkSize;

    PetscMemcpy( x + offset2, b + offset2, chunkSize*sizeof(PetscScalar) );

    //---------------------------------------------
    // Update the RHS with elements of the solution
    // vector solved in previous chunks
    //---------------------------------------------
    for ( int di = 0; di < center_diag; di++ ) {  //for ( int di = 0; di <= max_solved_diag; di++ ) {
      int d = tempDiags[di];
      int startBlock = chunk*chunkSize + d >= 0 ? 0 : -(chunk*chunkSize + d);
      int endBlock = chunkSize > -d ? -d : chunkSize;
      int offset1 = chunk*chunkSize + d;

      // offset into matrix coefficients to top of stripe in chunk
      int offset0 = di*numElements + chunk*chunkSize;
      // offset into already solved portion of iluSolveVec needed for this chunk/diagonal combination 
	
      for ( int block = startBlock; block < endBlock; block++ ) {
	sca1 = a->hostData[offset0 + block];
	sca2 = x[offset1 + block];
	x[offset2 + block] -= sca1 * sca2;
      }
    }

    for ( int block = 0; block < chunkSize; block++ ) {
      int offset2 = chunk*chunkSize + block;

      // update blocks below this diag block
      for ( int di = max_solved_diag+1; di < center_diag; di++ ) {
	int d = tempDiags[di];
	// is this thread in a diagonal block?
	if ( block - d < chunkSize ) {
	  int offset0 = di*numElements + chunk*chunkSize + block - d;
	  int offset1 = chunk*chunkSize + block - d;
	  x[offset1] -= a->hostData[offset0] * x[offset2];
	}
      }

    }

  }

  //-----------------------
  // UPPER TRIANGULAR SOLVE
  //-----------------------
  // PetscMemcpy( x, iluSolveVec, numElements*sizeof(PetscScalar) );    
  
  // proceed a chunk of rows at a time
  PetscInt min_solved_diag = center_diag+1;
  while ( (min_solved_diag < numDiags) && (tempDiags[min_solved_diag] < chunkSize ) )
    ++min_solved_diag;

  for ( chunk = numChunks-1; chunk >= 0; chunk-- ) {
    // offset into portion of x being solved for this chunk
    int offset2 = chunk*chunkSize;

    /* for ( int di = min_solved_diag; di < numDiags; di++ ) { */
    for ( int di = center_diag+1; di < numDiags; di++ ) {
      int d = tempDiags[di];
      int startBlock = d < chunkSize ? chunkSize - d : 0;
      int endBlock = (chunk+1)*chunkSize -1 + d < gridSize ? chunkSize - 1 : gridSize - chunk*chunkSize - 1 - d;
      // offset into already solved portion of iluSolveVec needed for this chunk/diagonal combination 
      int offset1 = chunk*chunkSize + d;
      int offset0 = (di+1)*numElements + chunk*chunkSize;
	
      for ( int block = endBlock; block >= startBlock; block-- )
	x[offset2 + block] -= a->hostData[offset0 + block] * x[offset1 + block];
    }

    for ( int block = chunkSize-1; block >= 0; block-- ) {
      int offset0 = (center_diag+1)*numElements + chunk*chunkSize + block;
      int offset2 = chunk*chunkSize + block;
      // solve diagBlock * x = beta
      //PetscMemcpy( beta, &(x[offset2]), dof*sizeof(PetscScalar) );
      x[offset2] *= a->hostData[offset0];

      // update blocks above this diag block
      for ( int di = min_solved_diag-1; di > center_diag; di-- ) {
	int d = tempDiags[di];
	// is this thread in a diagonal block?
	if ( block - d >= 0 ) {
	  int offset0 = (di+1)*numElements + chunk*chunkSize + block - d;
	  int offset1 = chunk*chunkSize + block;
	  int offset2 = chunk*chunkSize + block - d;
	  x[offset2] -= a->hostData[offset0] * x[offset1];
	}
      }

    }
  }

  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);


#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("cpu solve time %lf\n",elapsed);
#endif
  //ierr = PetscLogFlops(2*a->nz - A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatSolve_SeqSGGPU_cpu_simple"
//---
//
// Simple CPU Version
//
//---
PetscErrorCode MatSolve_SeqSGGPU_cpu_simple(Mat A,Vec bb,Vec xx)
{
  Mat_SeqSGGPU        *a = (Mat_SeqSGGPU*)A->data;
  PetscErrorCode    ierr;
  PetscInt          j, k, n=A->rmap->n, dof = a->dof;
  PetscInt          numDiags = a->diagonals->size();
  PetscScalar       *x;
  const PetscScalar *b;
  PetscScalar       sca1, sca2;
  PetscScalar       beta_scal, beta[MAXDOF];
  PetscInt          gridSize = a->m * a->n * a->p;
  PetscInt          numElements = gridSize * dof;
  PetscInt          chunkSize, chunk, numChunks;
  PetscBool         isseqgpu;
  PetscInt          tempDiags[MAXNUMDIAGS];
  int               offset0, offset1, offset2;

  PetscFunctionBegin;
  if (!n) PetscFunctionReturn(0);

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif


  if ( dumpVec ) {
    dumpVec = PETSC_FALSE;
    char  sgvecFname[64];
    PetscViewer rhsViewer;
    sprintf( sgvecFname, "sgvec-ex%d-%d-%d-%d-%d.bin", exNumber, a->m, a->n, a->p, a->dof );
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sgvecFname,FILE_MODE_WRITE, &rhsViewer);CHKERRQ(ierr);
    ierr = VecView(bb,rhsViewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&rhsViewer);CHKERRQ(ierr);
  }  

  // set chunkSize -- this determines how many rows of the solve to do as a unit
  chunkSize = a->m; //dim==3 ? a->p : dim==2 ? a->n : 1;
  for ( int di = 0; di < a->diagonals->size(); di++ )
    tempDiags[di] = (*a->diagonals)[di];

  // x is the input RHS, temp is the result of lower triangular solve
  ierr = PetscObjectTypeCompare((PetscObject)xx,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
  if (isseqgpu) {
    ierr = PetscObjectTypeCompare((PetscObject)bb,VECSEQGPU,&isseqgpu);CHKERRQ(ierr);
    if (!isseqgpu) {
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Both x and b must be same type");
    }
    Vec_SeqGPU *bd = (Vec_SeqGPU*) bb->data;
    /* synch up b */
    if (bd->syncState==VEC_GPU) {
      ierr = VecCopyOverD2H(bb,bd->cpuptr);CHKERRQ(ierr);
      bd->syncState=VEC_SYNCHED;
    }
  }
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);

  // LOWER TRIANGULAR SOLVE
  // determine the diagonals that contribute already solved values to the RHS for this thread
  PetscInt center_diag = numDiags / 2;

  numChunks = ((numElements/dof) + (chunkSize-1)) / chunkSize;
#if _TRACE
  printf("matSize = %d, chunkSize = %d, numChunks = %d, dof = %d\n", numElements, chunkSize, numChunks, dof);
#endif

  for ( chunk = 0; chunk < numChunks; chunk++ ) {

    //PetscMemcpy( x + chunk*chunkSize*dof, b + chunk*chunkSize*dof, chunkSize*dof*sizeof(PetscScalar) );

    if ( dof == 1 ) {

      for ( int block = 0; block < chunkSize; block++ ) {

	// offset2 is offset into portion of x being solved for this block
	offset2 = chunk*chunkSize + block;
	beta_scal = b[offset2];

	for ( int di = 0; di < center_diag; di++ ) {  //for ( int di = 0; di <= max_solved_diag; di++ ) {
	  int d = tempDiags[di];
	  // offset into already solved portion of iluSolveVec needed for this chunk/diagonal combination 
	  int offset1 = chunk*chunkSize + block + d;

	  if ( offset1 >= 0 ) {
	    sca2 = x[offset1];
	    // offset into matrix coefficients to top of stripe in chunk
	    int offset0 = di*numElements + chunk*chunkSize + block;
	    
	    sca1 = a->hostData[offset0];
	    beta_scal -= sca1 * sca2;
	  }
	}

	x[offset2] = beta_scal;
      }

    }

    else {
      for ( int block = 0; block < chunkSize; block++ ) {

	// offset2 is offset into portion of x being solved for this block
	offset2 = (chunk*chunkSize + block)*dof;
	PetscMemcpy( beta, b + offset2, dof*sizeof(PetscScalar) );

	for ( int di = 0; di < center_diag; di++ ) {  //for ( int di = 0; di <= max_solved_diag; di++ ) {
	  int d = tempDiags[di];
	  // offset into already solved portion of iluSolveVec needed for this chunk/diagonal combination 
	  int offset1 = (chunk*chunkSize + block + d)*dof;

	  if ( offset1 >= 0 )
	    for ( int rowCoord = 0; rowCoord < dof; rowCoord++ ) {
	      sca2 = x[offset1 + rowCoord];

	      // offset into matrix coefficients to top of stripe in chunk
	      int offset0 = (di*dof + rowCoord)*numElements + (chunk*chunkSize + block)*dof;
	
	      for ( int colCoord = 0; colCoord < dof; colCoord++ ) {
		// x[offset2 + block*dof + colCoord] -= a->hostData[ offset0 + block*dof + colCoord ] * x[offset1 + block*dof + rowCoord];
		sca1 = a->hostData[ offset0 + colCoord ];
		beta[colCoord] -= sca1 * sca2;
	      }
	    }
	}

	offset0 = center_diag*dof*numElements + (chunk*chunkSize + block)*dof;

	// solve diagBlock * x = beta
	for ( j = 0; j < dof; j++ )
	  x[offset2+j] = beta[j];
	for ( j = 0; j < dof-1; j++ ) {
	  for ( k = j+1; k < dof; k++ )
	    x[offset2 + k] += ( a->hostData[offset0 + j*numElements + k] * beta[j] );
	}

      }

    }

  }

  //-----------------------
  // UPPER TRIANGULAR SOLVE
  //-----------------------  

  // proceed a chunk of rows at a time
  for ( chunk = numChunks-1; chunk >= 0; chunk-- ) {

    if ( dof == 1 ) {

      for ( int block = chunkSize-1; block >= 0; block-- ) {

	// offset2 is offset into portion of x being solved for this block
	offset2 = chunk*chunkSize + block;
	beta_scal = x[offset2];

	for ( int di = center_diag+1; di < numDiags; di++ ) {
	  int d = tempDiags[di];
	  // offset into already solved portion of iluSolveVec needed for this chunk/diagonal combination 
	  offset1 = chunk*chunkSize + block + d;

	  if ( offset1 < numElements ) {
	    sca2 = x[offset1];
	    // offset into matrix coefficients to top of stripe in chunk
	    offset0 = (di+1)*numElements + chunk*chunkSize + block;
	    
	    sca1 = a->hostData[offset0];
	    beta_scal -= sca1 * sca2;
	  }
	}

	offset0 = (center_diag+1)*numElements + chunk*chunkSize + block;
	x[offset2] = a->hostData[offset0] * beta_scal;
      }
    }

    else {
      for ( int block = chunkSize-1; block >= 0; block-- ) {
	// offset into portion of x being solved for this chunk
	offset2 = (chunk*chunkSize + block)*dof;
	PetscMemcpy( beta, x + offset2, dof*sizeof(PetscScalar) );

	for ( int di = center_diag+1; di < numDiags; di++ ) {
	  int d = tempDiags[di];
	  offset1 = (chunk*chunkSize + block + d)*dof;

	  if ( offset1 < numElements )
	    for ( int rowCoord = 0; rowCoord < dof; rowCoord++ ) {

	      // offset into matrix coefficients to top of stripe in chunk
	      offset0 = ((di+1)*dof + rowCoord)*numElements + (chunk*chunkSize + block)*dof;
	      
	      for ( int colCoord = 0; colCoord < dof; colCoord++ ) {
		beta[colCoord] -= a->hostData[ offset0 + colCoord ] * x[offset1 + rowCoord];
	      }
	    }
	}

	offset0 = (center_diag+1)*dof*numElements + (chunk*chunkSize + block)*dof;

	// solve diagBlock * x = beta
	for ( j = 0; j < dof; j++ ) {
	  x[offset2 + j] = a->hostData[offset0 + j*numElements + j] * beta[j];
	  for ( k = 0; k < j; k++ )
	    x[offset2 + k] += ( a->hostData[offset0 + j*numElements + k] * beta[j] );
	}
      }

    }
  }

  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);


#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("cpu solve time %lf\n",elapsed);
#endif
  //ierr = PetscLogFlops(2*a->nz - A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "MatGetFactor_seqsggpu_petsc"
PetscErrorCode MatGetFactor_seqsggpu_petsc(Mat A,MatFactorType ftype,Mat *B)
{
  PetscInt           n = A->rmap->n;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  dumpMat = PETSC_FALSE;
  dumpVec = PETSC_FALSE;

  if ( dumpMat ) {
    dumpMat = PETSC_FALSE;
    PetscViewer sgmatViewer;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"sgmat.full.bin",FILE_MODE_WRITE, &sgmatViewer);CHKERRQ(ierr);
    ierr = MatView(A,sgmatViewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&sgmatViewer);CHKERRQ(ierr);
  }

  ierr = MatCreate(((PetscObject)A)->comm,B);CHKERRQ(ierr);
  ierr = MatSetSizes(*B,n,n,n,n);CHKERRQ(ierr);
  ierr = MatSetType(*B,MATSEQSGGPU);CHKERRQ(ierr);
  if (!(*B)->preallocated) {
    Mat_SeqSGGPU* a = (Mat_SeqSGGPU*)A->data;
    PetscInt       dims[3], *starts;
    dims[0] = a->m;
    dims[1] = a->n;
    dims[2] = a->p;
    starts = (PetscInt*)malloc(sizeof(PetscInt)*a->dim);
    ierr = MatSetStencil(*B,a->dim,dims,starts,a->dof); CHKERRQ(ierr);
    ierr = MatSeqSGGPUSetPreallocation(*B,0,a->dof);CHKERRQ(ierr);
  }

  if (ftype == MAT_FACTOR_LU || ftype == MAT_FACTOR_ILU || ftype == MAT_FACTOR_ILUDT){
    (*B)->ops->ilufactorsymbolic = MatILUFactorSymbolic_SeqSGGPU;
    // (*B)->ops->lufactorsymbolic  = MatLUFactorSymbolic_SeqSGGPU;
  } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Factor type not supported");
  
  (*B)->factortype = ftype;
  
  PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatConvertLU_SeqSGGPU_SeqAIJ"
/*
  MatConvertLU_SeqSGGPU_SeqAIJ - Converts from an LU-factored sggpu format to two seqaij.
  John Eisenlohr*/
PetscErrorCode MatConvertLU_SeqSGGPU_SeqAIJ(Mat A, Mat *AIJ_L, Mat *AIJ_U ){
  //printf(".................MatConvertLU_SeqSGGPU_SeqAIJ() called\n");
  PetscFunctionBegin;
  Mat_SeqSGGPU *a = (Mat_SeqSGGPU *) A->data;
  PetscScalar *hostDataCopy, *hostDataOrig;
  Mat B,C;
  PetscInt i,j, kl, ku;
  PetscInt m = A->rmap->n,n = A->cmap->n;
  PetscScalar *vals, *vals_l, *vals_u;
  PetscInt *cols, *cols_l, *cols_u;
  PetscErrorCode ierr;
  PetscInt *nnz_l, *nnz_u, nnz;
  PetscInt blockRow, rowInBlock;
  PetscInt dof = a->dof, numBlocks = a->m * a->n * a->p;
  PetscInt numDiags = a->diagonals->size();

  // For efficiency, the diagonal blocks of both the
  // lower and upper factors of A have been inverted
  // but we want to return the factors with univerted diagonal blocks.
  // So we copy the host data from A, invert the diagonals and
  // temporarily replace A's host data with the modified data
  ierr = PetscMalloc(numBlocks*dof*dof*numDiags*sizeof(PetscScalar),&hostDataCopy);CHKERRQ(ierr);
  ierr = PetscMemcpy( hostDataCopy, a->hostData, numBlocks*dof*dof*numDiags*sizeof(PetscScalar) );CHKERRQ(ierr);
  ierr = InvertFactoredDiagBlocks(a,hostDataCopy);CHKERRQ(ierr);
  hostDataOrig = a->hostData;
  a->hostData = hostDataCopy;

  // arrays to hold number of zeros in lower and upper portions of each row
  ierr = PetscMalloc(m*sizeof(PetscInt),&nnz_l); CHKERRQ(ierr);
  ierr = PetscMalloc(m*sizeof(PetscInt),&nnz_u); CHKERRQ(ierr);

  // workspace for holding the upper and lower portions of each row
  // these are reused for each row
  ierr = PetscMalloc(n*sizeof(PetscInt),&cols_l); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(PetscInt),&cols_u); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(PetscScalar),&vals_l); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(PetscScalar),&vals_u); CHKERRQ(ierr);

  // count non-zeros in upper and lower triangles
  // there is surely a closed-form way to do this
  for(i=0;i<m;i++){
    blockRow = i / a->dof;
    nnz_l[i]=0;
    nnz_u[i]=0;
    for (j = 0; j < a->diagonals->size(); ++j) {
      int d = (*a->diagonals)[j];
      if ( ((d + blockRow) >= 0) && ((d + blockRow) < numBlocks) ) {
	if ( d < 0 ) // block is below diagonal
	  nnz_l[i] += a->dof;
	else if ( d > 0 ) // block is above diagonal
	  nnz_u[i] += a->dof;
	else { // block is on diagonal
	  rowInBlock = i % a->dof;
	  nnz_l[i] += rowInBlock + 1;
	  nnz_u[i] += (a->dof - rowInBlock);
	}
      }
    }
  }

  printf("m: %d, n: n: %d\n",m,n);
  ierr = MatCreateSeqAIJ(((PetscObject)A)->comm,m,n,PETSC_NULL,nnz_l,&B);CHKERRQ(ierr);
  ierr = MatCreateSeqAIJ(((PetscObject)A)->comm,m,n,PETSC_NULL,nnz_u,&C);CHKERRQ(ierr);

  // get each row from sggpu, divide into L and R pieces,
  // set values in L and R matrices
  for(i=0;i<m;i++){
      nnz = 0;
      ierr = MatGetRow_SeqSGGPU(A,i, &nnz, &cols,&vals); CHKERRQ(ierr);
      kl = 0;
      ku = 0;
      for ( j = 0; j < nnz; j++ ) {
	if ( cols[j] < i ) {
	  vals_l[kl] = vals[j];
	  cols_l[kl++] = cols[j];
	}
	else {
	  vals_u[ku] = vals[j];
	  cols_u[ku++] = cols[j];
	}
      }
      // add the 1 on the diagonal of L
      vals_l[kl] = 1.0;
      cols_l[kl] = i;

      ierr = MatRestoreRow_SeqSGGPU(A,i,&nnz,&cols,&vals);CHKERRQ(ierr);

      ierr = MatSetValues(B,1,&i,nnz_l[i],cols_l,vals_l,INSERT_VALUES);
      ierr = MatSetValues(C,1,&i,nnz_u[i],cols_u,vals_u,INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr=PetscFree(nnz_l);CHKERRQ(ierr);
  ierr=PetscFree(nnz_u);CHKERRQ(ierr);
  ierr=PetscFree(cols_l);CHKERRQ(ierr);
  ierr=PetscFree(cols_u);CHKERRQ(ierr);
  ierr=PetscFree(vals_l);CHKERRQ(ierr);
  ierr=PetscFree(vals_u);CHKERRQ(ierr);
  a->hostData = hostDataOrig;
  ierr=PetscFree(hostDataCopy);CHKERRQ(ierr);

  *AIJ_L = B;
  *AIJ_U = C;

  PetscFunctionReturn(0);
}
EXTERN_C_END
