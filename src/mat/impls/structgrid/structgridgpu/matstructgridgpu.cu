/*  -------------------------------------------------------------------- 

     This file extends structgrid data type to make use of GPUS. The new data type
     is structgridgpu. The implementation of the new datatype emulates the seqaijcusp
     implementation which is an extension to aij matrix format. 
     Author: Chekuri S. Choudary, RNET
 */


#define PETSCMAT_DLL
#include "../src/mat/impls/structgrid/structgridgpu/matstructgridgpu.h"

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>
#include "../src/mat/impls/structgrid/matstructgrid.h"

#include "private/matimpl.h"
#include "matstructgridgpu.h"

#define _DBGFLAG 0
//#define PRINT
//block size is 1x256. 
#define BLOCKWIDTH_X 256		
#define BLOCKWIDTH_Y 1   


// ----------------------------------------------------------
// hardcodiing the shared memory size this should be set
// to give maximum performance, however should be
// replaced soon with a more flexable dynamically allocated
// shared memory scheme
// written by: dlowell ANL-MCS
// ----------------------------------------------------------
#define SHDSIZE 4


#define size 64
// Structure for Constant Device memory
// storing constants and indices and index limits
// stencile size is hard coded
// written by: dlowell ANL-MCS
// -----------------------------------------------
#define STLSIZE 64
struct Stencilparams{
  int m;
  int n;
  int p;
  int vecsize_x;
  int vecsize_y;
  int matsize;
  int nos;
  int dof;
  int lda1;
  int lda2;
  int lda3;
  int idx[STLSIZE];
  int idy[STLSIZE];
  int idz[STLSIZE];
  int tile_x;
  int tile_y;
  int tile_z;
  int tsizex;
  int tsizey;
  int tsizez;
};//836 bytes

__constant__ Stencilparams devparams;//device memory


static double* devA;
static PetscScalar* d_coeff;
static double* devX;
static double* devY;




// ----------------------------------------------------------
// helper function for error checking
// pops the CUDA error stack and exits on nonzero error code
// written by: dlowell ANL-MCS
// ----------------------------------------------------------
void checkCUDAError(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) ); 
    exit(EXIT_FAILURE); 
  }
} 

//------------------------------------------------------
// general timer function using unix system call
// dlowell ANL-MCS
//------------------------------------------------------
double getclock(){
  struct timezone tzp;
  struct timeval tp;
  gettimeofday (&tp, &tzp);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}


/*  --------------------------------------------------------------------
     This function destroys the matrix of type structgridgpu. It first 
     deallocates the memory on GPU and then calls the MatDestroy_SeqSG 
     function.
     Author: Chekuri S. Choudary, RNET
 */

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatDestroy_SeqSGGPU"
PetscErrorCode  MatDestroy_SeqSGGPU(Mat B)
{
  //printf("Call to MatDestroy_SeqSGGPU(Mat B)\n");
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (B->valid_GPU_matrix != PETSC_CUSP_UNALLOCATED) 
  {
    if (devA) cudaFree(devA);
    if (d_coeff) cudaFree(d_coeff);
    if(devY) cudaFree(devY);
    if(devX) cudaFree(devX);
  }

  B->valid_GPU_matrix = PETSC_CUSP_UNALLOCATED;

  ierr             = MatDestroy_SeqSG(B);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
EXTERN_C_END



/*  --------------------------------------------------------------------
     This function creates a datatype of structgridgpu. It first creates a
     structgrid datatype and overrides the matrix multiplication method.
     Author: Chekuri S. Choudary, RNET
 */

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatCreate_SeqSGGPU"
PetscErrorCode  MatCreate_SeqSGGPU(Mat B)
{

  printf("Call to MatCreate_SeqSGGPU(Mat B)\n");

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr             = MatCreate_SeqSG(B);CHKERRQ(ierr);
  B->ops->mult     = MatMult_SeqSGGPU;
  B->ops->destroy  = MatDestroy_SeqSGGPU;

  ierr = PetscObjectChangeTypeName((PetscObject)B,MATSTRUCTGRIDGPU);CHKERRQ(ierr);
  B->valid_GPU_matrix = PETSC_CUSP_UNALLOCATED;
  PetscFunctionReturn(0);
}
EXTERN_C_END


//---------------------------------------------------------------------
//     This function implements matrix vector multiplication for the
//     structgridgpu datatype. It calls a CUDA kernel to do matrix
//     multiplication on the GPU.
//     Author: Daniel Lowell, ANL-MCS, Chekuri S. Choudary, RNET
//---------------------------------------------------------------------
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqSGGPU"
PetscErrorCode MatMult_SeqSGGPU(Mat mat, Vec x, Vec y)
{
  int i;
  PetscErrorCode ierr;
  Mat_SeqSG * a = (Mat_SeqSG *) mat->data;
  PetscScalar * v = a->a, *xx,*yy;

  PetscFunctionBegin;
  ierr = VecSet(y,0.0); CHKERRQ(ierr);
  ierr = VecGetArray(x, &xx); CHKERRQ(ierr);
  ierr = VecGetArray(y, &yy); CHKERRQ(ierr);

  /* Call to Jeswin's version */
  ierr = SGCUDA_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,
      a->stpoints,&(mat->valid_GPU_matrix),a->dof);CHKERRQ(ierr);

  ierr = VecRestoreArray(x,&xx); CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&yy); CHKERRQ(ierr);
  ierr = PetscLogFlops(2*a->nz*a->stpoints); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
EXTERN_C_END


/*  -------------------------------------------------------------------- 
     The following is a CUDA kernel for matrix vector multiplication on 
     the GPU. The matrix is in a custom layout that facilitates better 
     memory accesses and vectorization. 
     Author: Chekuri S. Choudary, RNET
 */
/* Version with Shared memory for X only supports rectangular tiles. */
/* __global__ void MatMult_Kernel(PetscScalar * ptr_coeff, PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt *idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos)
{

int tx= blockDim.x * blockIdx.x + threadIdx.x;
int ty= blockDim.y * blockIdx.y + threadIdx.y;

int l,i;
int xdisp,ydisp,zdisp,offset;
int lda1=m*n*p,lda2=m*n,lda3=m;

__shared__ PetscScalar y_sm[256];

// initializing to the zero
y_sm[threadIdx.y*BLOCKWIDTH_X + threadIdx.x]=0;
for (l=0;l<nos;l++)
	{
	xdisp = idx[l]; ydisp = idy[l]; zdisp = idz[l]; offset = l*lda1;
	if (tx > n-1)
	{
	break; //use Break and test performance later(divergence)
	}
	if (ty > m-1)
	{
	break; //use Break and test performance later(divergence)
	}
	if (l==1 && tx==n-1 && ty==m-1)
	{
	continue;
	}
	if (l==2 && tx==0 && ty==0)
	{
	continue;
	}
	if (l==3 && ty==m-1)
	{
	continue;
	}
	if (l==4 && ty==0)
	{
	continue;
	}
	for(i=0;i<p;i++)
	y_sm[threadIdx.y*BLOCKWIDTH_X + threadIdx.x]+= (ptr_coeff[offset + i*lda2 + ty*lda3 +tx] * ptr_x[(i+zdisp)*lda2 + (ty+ydisp)*lda3 + (tx+xdisp)]);
	}

	ptr_y[ty*lda3 + tx]= y_sm[threadIdx.y*BLOCKWIDTH_X + threadIdx.x];
}
 */
/*  #define BLOCKWIDTH 8
#define BLOCKWIDTH_X 8
#define BLOCKWIDTH_Y 8
#define BLOCKWIDTH_Z 8 

  __global__ void MatMult_Kernel(PetscScalar * ptr_coeff, PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt *idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos)
{

int tx= blockDim.x * blockIdx.x + threadIdx.x;
int ty= blockDim.y * blockIdx.y + threadIdx.y;
int tz= blockDim.z * blockIdx.z + threadIdx.z;
int l,i;
int xdisp,ydisp,zdisp,offset;
int lda1=m*n*p,lda2=m*n,lda3=m;

__shared__ PetscScalar y_sm[512];

// initializing to the zero
y_sm[threadIdx.z*BLOCKWIDTH_X*BLOCKWIDTH_Y + threadIdx.y*BLOCKWIDTH_X + threadIdx.x]=0;
for (l=0;l<nos;l++)
	{
	xdisp = idx[l]; ydisp = idy[l]; zdisp = idz[l]; offset = l*lda1;
	if (tx > n-1)
	{
	break; //use Break and test performance later(divergence)
	}
	if (ty > m-1)
	{
	break; //use Break and test performance later(divergence)
	}
	if (tz > p-1)
	{
	break;
	}
	if (l==1 && tx==n-1 && ty==m-1 && tz==p-1)
	{
	continue;
	}
	if (l==2 && tx==0 && ty==0 && tz==0)
	{
	continue;
	}
	if (l==3 && ty==m-1)
	{
	continue;
	}
	if (l==4 && ty==0)
	{
	continue;
	}
	if (l==5 && tz==p-1)
	{
	continue;
	}
	if (l==6 && tz==0)
	{
	continue;
	}
	//for(i=0;i<p;i++)
	y_sm[threadIdx.z*BLOCKWIDTH_X*BLOCKWIDTH_Y + threadIdx.y*BLOCKWIDTH_X + threadIdx.x]+= (ptr_coeff[offset + tz*lda2 + ty*lda3 +tx] * ptr_x[(tz+zdisp)*lda2 + (ty+ydisp)*lda3 + (tx+xdisp)]);
	}

	ptr_y[tz*lda2+ ty*lda3 + tx]= y_sm[threadIdx.z*BLOCKWIDTH_X*BLOCKWIDTH_Y +threadIdx.y*BLOCKWIDTH_X + threadIdx.x];
}
 */   


//------------------------------------------------------------------------------------
//   These functions are used to bind and unbind the Vector x to the texture Memory.	   
//------------------------------------------------------------------------------------ 
texture<int2, 1> tex_x_double;

void unbind_x( double * x)
{
  cudaUnbindTexture(tex_x_double);
}


static __inline__ __device__ double fetch_double(texture<int2, 1> tex_x_double, int i)
{
  int2 v = tex1Dfetch(tex_x_double,i);
  return __hiloint2double(v.y, v.x);
}

//------------------------------------------------------------------------------------
// Dynamically allocating Shared Memory size	   
//------------------------------------------------------------------------------------ 

extern __shared__ PetscInt idx_sm[];

//------------------------------------------------------------------------------------
//  Below functions are SPMV kernel functions where x through the Texture Memory, offsets are accesed
//	through the Shared Memory, Y is accessed per thread from registers. Coeff accesses are from the global Memory 
//	but they are coalesced.  	   
//------------------------------------------------------------------------------------ 


__global__ void MatMul_Kernel_tex_1_DOF(PetscScalar * ptr_coeff,
    PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt *idx, PetscInt m,
    PetscInt n, PetscInt p, PetscInt nos, PetscInt DOF) {
  int tx = blockDim.x * blockIdx.x + threadIdx.x;
  int l, offset;
  int lda1 = m * n * p * DOF;
  int lda2 = m * p * DOF; //lda3=m*DOF
  PetscInt Index;
  PetscScalar y_reg = 0;

  if (threadIdx.x < nos) {
    idx_sm[threadIdx.x] = idx[threadIdx.x];
  }

  int reg2 = blockIdx.y * lda2 + tx;

  //Iterating through the Diagonals
  for (l = 0; l < nos; l++) {
    Index = reg2 + idx_sm[l];

    if (Index >= 0 && Index < lda1) {
      offset = l * lda1;
      y_reg += ptr_coeff[offset + reg2] * fetch_double(tex_x_double, Index);
    }
  }

  ptr_y[reg2] = y_reg;
}


// Jeswin's uncommitted kernel
__global__ void MatMult_Kernel(PetscScalar * ptr_coeff, PetscScalar* ptr_x,
    PetscScalar* ptr_y, PetscInt *idx, PetscInt m, PetscInt n, PetscInt p,
    PetscInt nos, PetscInt stpoints, PetscInt DOF) {
  int tx = blockDim.x * blockIdx.x + threadIdx.x;
  int l, i;
  int offset;
  int lda1 = m * n * p * DOF, lda2 = m * p * DOF, lda3 = m * DOF;
  int Index;
  double y_reg = 0;

  if (tx >= lda1)
    return;

  if (threadIdx.x < stpoints) {
    idx_sm[threadIdx.x] = idx[threadIdx.x];
  }
  __syncthreads();

  int reg2 = blockIdx.y * lda2 + tx;

  //Divergence in Computation
  for (l = 0; l < stpoints; l++) {
    Index = reg2 + idx_sm[l];
    if ((Index >= 0) && (Index < lda1)) {
      for (i = 0; i < DOF; i++) {
        offset = (l * DOF + i) * lda1;
        y_reg += ptr_coeff[offset + reg2]
            * fetch_double(tex_x_double, Index - (threadIdx.x % DOF) + i);
      }
    }
  }
  ptr_y[reg2] = y_reg;
}



// Jeswin's kernel as of HG-tip
#if 0
__global__ void MatMul_Kernel_tex(double * ptr_coeff, double* ptr_x,
    double* ptr_y, PetscInt *idx, PetscInt m, PetscInt n, PetscInt p,
    PetscInt nos, PetscInt DOF) {

  int tx = blockDim.x * blockIdx.x + threadIdx.x;
  int l, i, offset;
  int lda1 = m * n * p * DOF, lda2 = m * p * DOF; //lda3=m*DOF
  PetscInt X_Index, Index;
  double y_reg = 0;
  int BAND_SIZE = (DOF - 1) * 2 + 1;

  if (threadIdx.x < nos) {
    idx_sm[threadIdx.x] = idx[threadIdx.x];
  }

  int reg2 = blockIdx.y * lda2 + tx;

  //Iterating through the Diagonals
  for (l = 0; l < nos; l++) {
    X_Index = reg2 + idx_sm[l];

    if (X_Index >= 0 && X_Index < lda1) {
      for (i = 0; i < BAND_SIZE; i++) {
        offset = (l * BAND_SIZE + i) * lda1;

        if (i > DOF - 1) {
          Index = X_Index - (i - (DOF - 1));
          if (Index < 0) {
            continue;
          } else {
            y_reg += ptr_coeff[offset + reg2]
                * fetch_double(tex_x_double, Index);
          }
        } else {
          Index = X_Index + i;

          y_reg += ptr_coeff[offset + reg2] * fetch_double(tex_x_double, Index);
        }
      }
    }
  }
  ptr_y[reg2] = y_reg;
}
#endif

//------------------------------------------------------------------------------------
//   The function is a wrapper function which sets up the device memory, transfers
//   data to and from the device, and calls the MatMult kernel. 
//------------------------------------------------------------------------------------ 

int SGCUDA_MatMult(PetscScalar* coeff, PetscScalar* x, PetscScalar* y,
    PetscInt *idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n,
    PetscInt p, PetscInt stpoints, PetscCUSPFlag* fp, PetscInt DOF) {
  double tbegin1, tbegin2, tend1, tend2;
  static PetscInt size_coeff; 
  double tsetup,tkernel;
  static unsigned int kcalls=0;
  PetscInt size_xy, size_id; 
  static double temp=0;
  PetscScalar* d_x;
  PetscScalar* d_y;
  PetscInt *d_linear_idx;
  int nos;

  int BLOCK_SIZE;
  int cons=m*DOF;
  int cons1=m*n*DOF;

  nos = DOF*stpoints;

#if _DBGFLAG
  printf("m: %d  n: %d  p: %d  nos: %d  stpoints: %d  DOF: %d\n", m, n, p, nos, stpoints, DOF);

  for (unsigned i = 0; i < stpoints; ++i) {
    printf("Stencil Point:  %d, %d, %d\n", idx[i], idy[i], idz[i]);
  }
#endif

  PetscInt *linear_idx;
  PetscMalloc(stpoints * sizeof(PetscInt), &linear_idx);

  for (unsigned i = 0; i < stpoints; ++i) {
    linear_idx[i] = idx[i]*DOF + idy[i]*cons + idz[i]*cons1;
  }
  //----Single offset instead of using three offsets int the x,y and z direction.  
  //if (stpoints==5)
  //{
  //  idx[0]=0;
  //  idx[1]=DOF;
  //  idx[2]=cons;
  //  idx[3]=-DOF;
  //  idx[4]=-cons;
  //}
  //else if (stpoints==7)
  //{
  //  idx[0]=0;
  //  idx[1]=DOF;
  //  idx[2]=cons;
  //  idx[3]=cons1;
  //  idx[4]=-DOF;
  //  idx[5]=-cons;
  //  idx[6]=-cons1;
  //} else {
  //  printf("Bad value for stpoints\n");
  //  exit(EXIT_FAILURE);
  //}


  //------------Printing Matrices for Debugging---------------------------------
#ifdef PRINT
  printf("offset vector\n");
  for (int i=0;i<nos;i++)
  {
    printf("%d  ", idx[i]);
  }
  printf("\n");

  printf("Matrix X\n");
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<m*DOF;j++)
    {
      printf("%0.2f   ",x[i*cons+j]);
    }
    //printf("\n");
  }
  printf("\n");
  printf("Matrix A\n");
  for(int s=0;s<DOF*nos;s++)
  {
    printf("\n");
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<m*DOF;j++)
      {
        printf("%0.2f   ",coeff[s*cons1+i*cons+j]);
      }
      printf("\n");
    }

  }
#endif
  //------------------------------------------------------------------------------------------


#if(_DBGFLAG) 
  tbegin1=getclock();
#endif

  if ((*fp == PETSC_CUSP_UNALLOCATED) ||
      (*fp == PETSC_CUSP_CPU) )
  {
    if (*fp == PETSC_CUSP_UNALLOCATED)
    {
      size_coeff=nos*m*n*p*DOF*sizeof(PetscScalar);
      cudaMalloc((void**)&d_coeff,size_coeff);
      checkCUDAError("cudaMalloc (d_coeff)");

      //cudastatus0=cudaMalloc((void**)&devA,matsize);
      //if(cudastatus0!=cudaSuccess)
      //	{
      //  printf("Error in devA memory allocation:\nstatus0: %s\n",
      //	cudaGetErrorString(cudastatus0));
      //  PetscFunctionReturn(PETSC_ERR_MEM);
      //	}
    }

    cudaMemcpy(d_coeff, coeff, size_coeff, cudaMemcpyHostToDevice);
    checkCUDAError("cudaMemcpy (d_coeff)");

    //cudastatus1=cudaMemcpy(devA,A,matsize,cudaMemcpyHostToDevice);
    //if(cudastatus1!=cudaSuccess)
    //{
    // if(devA) cudaFree(devA);
    //  printf("Error in devA memory copying:\nstatus1: %s\n",
    //	cudaGetErrorString(cudastatus1));
    //  PetscFunctionReturn(PETSC_ERR_MEM);
    //}

    *fp = PETSC_CUSP_BOTH;
  }


  //size_coeff=nos*m*n*p*sizeof(PetscScalar);
  //cudaMalloc((void**)&d_coeff,size_coeff);
  //cudaMemcpy(d_coeff, coeff, size_coeff, cudaMemcpyHostToDevice);


  size_xy = m*n*p*DOF*sizeof(PetscScalar);
  cudaMalloc((void**)&d_x,size_xy);
  checkCUDAError("cudaMalloc (d_x)");
  cudaMemcpy(d_x, x, size_xy, cudaMemcpyHostToDevice);
  checkCUDAError("cudaMemcpy (d_x)");

  cudaMalloc((void**)&d_y,size_xy);
  checkCUDAError("cudaMalloc (d_y)");
  cudaMemcpy(d_y, y, size_xy, cudaMemcpyHostToDevice);
  checkCUDAError("cudaMemcpy (d_y)");

  size_id = stpoints*sizeof(PetscInt);
  cudaMalloc((void**)&d_linear_idx,size_id);
  checkCUDAError("cudaMalloc (d_idx)");
  cudaMemcpy(d_linear_idx, linear_idx, size_id, cudaMemcpyHostToDevice);
  checkCUDAError("cudaMemcpy (d_idx)");


  //Binding X to the texture Memory
  cudaBindTexture(0, tex_x_double, d_x, size_xy);
  checkCUDAError("cudaBindTexture");

#if(_DBGFLAG)
  //cudaPrintfInit();
  tend1=getclock();
  tsetup=tend1-tbegin1;
  tbegin2=getclock();
#endif


  // Kernel Setup and Configurations

  if ((DOF%2)!=0 || (DOF==6))
  {
    BLOCK_SIZE=BLOCKWIDTH_X-BLOCKWIDTH_X%DOF;
  }
  else{
    BLOCK_SIZE=BLOCKWIDTH_X;
  }

  //dim3 dimBlock(BLOCK_SIZE,BLOCKWIDTH_Y);
  //dim3 dimGrid((int)ceil((float)(m*p*DOF)/(float)BLOCK_SIZE),((int)ceil((float)(n)/(float)BLOCKWIDTH_Y)));

  dim3 dimBlock(BLOCK_SIZE,BLOCKWIDTH_Y);
  dim3 dimGrid((int)ceil((float)(m*n*p*DOF)/(float)BLOCK_SIZE), 1);

  PetscInt shared_size = stpoints * sizeof(PetscInt);

#if _DBGFLAG
  printf("Launch Bounds:\n");
  printf("dimBlock: %d, %d\n", dimBlock.x, dimBlock.y);
  printf("dimGrid:  %d, %d\n", dimGrid.x, dimGrid.y);

  for (unsigned i = 0; i < stpoints; ++i) {
    printf("linear_idx[%d] = %d\n", i, linear_idx[i]);
  }
#endif

  //if (DOF==1)
 // {
  //  MatMul_Kernel_tex_1_DOF<<<dimGrid,dimBlock,shared_size>>>(d_coeff, d_x, d_y, d_idx, m, n, p, nos,DOF);
  //}
  //else{
    MatMult_Kernel<<<dimGrid,dimBlock,shared_size>>>(d_coeff, d_x, d_y, d_linear_idx, m, n, p, nos, stpoints, DOF);
  //}

  // check if kernel execution generated and error
  //cutilCheckMsg("Kernel execution failed");

  /*
	   if (m > BLOCKWIDTH){
	   // dim3 dimBlock(BLOCKWIDTH,BLOCKWIDTH);
	   // dim3 dimGrid((int)ceil((float)m/(float)BLOCKWIDTH),((int)ceil((float)n/(float)BLOCKWIDTH)));

	   dim3 dimBlock(BLOCKWIDTH,BLOCKWIDTH,BLOCKWIDTH);
	   dim3 dimGrid((int)ceil((float)m/(float)BLOCKWIDTH),((int)ceil((float)n/(float)BLOCKWIDTH)),p/BLOCKWIDTH);

	   // cutilCheckError(cutCreateTimer(&timer));
	   // cutilCheckError(cutStartTimer(timer));

	   MatMult_Kernel<<<dimGrid,dimBlock>>>(d_coeff, d_x, d_y, d_idx, d_idy, d_idz, m, n, p, nos);

	   }
	   else
	   {
	   //dim3 dimBlock(m,n);
	   dim3 dimBlock(m,n,p);
	   dim3 dimGrid(1,1,1);

	   // cutilCheckError(cutCreateTimer(&timer));
	   // cutilCheckError(cutStartTimer(timer));

	   MatMult_Kernel<<<dimGrid,dimBlock>>>(d_coeff, d_x, d_y, d_idx, d_idy, d_idz, m, n, p, nos);


	   }

   */


#if(_DBGFLAG) 
  cudaDeviceSynchronize();
  checkCUDAError("Launch/Sync");
  tend2=getclock();
  tkernel=tend2-tbegin2;
  //cudaPrintfDisplay(stdout, true);
  //cudaPrintfEnd();
#endif

  //Read y from the Device Memory

  cudaMemcpy(y, d_y, size_xy, cudaMemcpyDeviceToHost);
  checkCUDAError("cudaMemcpy (d_y) OUT");
  /*
	  int i;
	  char *fn = "/homes/dlowell/cudaexprs/dcheck/outfile_SG1.txt";
	  FILE *fptr;
	  fptr=fopen(fn,"a");
	  if(!fptr){
	  printf("file pointer error.\n");
	  PetscFunctionReturn(PETSC_ERR_LIB);
	  }else{
	  //printf("yy->map->n: %d\n",yy->map->n);
	  for(i=0;i<m*n*p*DOF;i++){
	  //printf("printed to file: %d\n",i);
	  if(y[i]!=0.)fprintf(fptr,"%f ",y[i]);
	  }
	  fclose(fptr);
	  }
   */

#if(_DBGFLAG)
  temp+=tkernel;
  if (kcalls==0)
  {
    printf("\n Structured Grid MatrixMul Kernel Permormance for m *%d* and n size *%d* \n",m,n);
  }
  //if (kcalls==1000)
  {
    printf("\ncopy time (sec) : %f\n",tsetup);
    printf("Kernel time (sec): %f\n",tkernel);
    printf("Performance in Megaflops with for %dth Kernel call\n",kcalls);
    printf("Performance in Megaflops with copy time = %f\n",(2*nos*n*m*p*1.0e-6)/(tsetup+tkernel));
    printf("Performance in Megaflops without copy time = %f\n",(2*nos*n*m*p*1.0e-6)/tkernel);
    printf("Culmative Performance in Megaflops for *%d* calls without copy time = %f\n",kcalls,(2*nos*n*m*p*1.0e-6)/(temp/(kcalls+1)));
  }
#endif	
  kcalls++;


#ifdef PRINT
  for(int i=0;i<m*n;i++)
    printf("Y[%d]: %lf\n",i,y[i]);
#endif

  PetscFree(linear_idx);

  //Free Device Memory
  //cudaFree(d_coeff);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_linear_idx);

  return 0;
}




