
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
#include "../src/mat/impls/structgrid/matstructgrid.h"

#include "private/matimpl.h"
#include "matstructgridgpu.h"


#define size 64


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
  PetscErrorCode ierr;
  PetscFunctionBegin;
  
  ierr             = MatCreate_SeqSG(B);CHKERRQ(ierr);
  B->ops->mult     = MatMult_SeqSGGPU;
  
  ierr = PetscObjectChangeTypeName((PetscObject)B,MATSTRUCTGRIDGPU);CHKERRQ(ierr);
  B->valid_GPU_matrix = PETSC_CUSP_UNALLOCATED;
  PetscFunctionReturn(0);
}
EXTERN_C_END




/*  -------------------------------------------------------------------- 
     This function implements matrix vector multiplication for the 
     structgridgpu datatype. It calls a CUDA kernel to do matrix 
     multiplication on the GPU.  
     Author: Chekuri S. Choudary, RNET
*/
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqSGGPU"
PetscErrorCode MatMult_SeqSGGPU(Mat mat, Vec x, Vec y)
{
	PetscErrorCode ierr;
	Mat_SeqSG * a = (Mat_SeqSG *) mat->data;
	PetscScalar * v = a->a, *xx,*yy;
	
	PetscFunctionBegin;
	ierr = VecSet(y,0.0); CHKERRQ(ierr);
	ierr = VecGetArray(x, &xx); CHKERRQ(ierr);
	ierr = VecGetArray(y, &yy); CHKERRQ(ierr);

ierr = SGCUDA_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->stpoints); 
CHKERRQ(ierr);

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
__global__ void MatMult_Kernel(PetscScalar * ptr_coeff, PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt *idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos)
{
int tx=  blockDim.x * blockIdx.x + threadIdx.x;
int ty=  blockDim.y * blockIdx.y + threadIdx.y;
int l,i;
int xdisp,ydisp,zdisp,offset;
int lda1=m*n*p,lda2=m*n,lda3=m;

for (l=0;l<nos;l++)
        {
        xdisp = idx[l]; ydisp = idy[l]; zdisp = idz[l]; offset = l*lda1;
        if (l==1 && tx==size-1 && ty==size-1)
        {
        	continue;
        }
        if (l==2 && tx==0 && ty==0)
        {
        	continue;
        }
        if (l==3 && ty==size-1)
        {
        	continue;
        }
        if (l==4 && ty==0)
        {
        	continue;
        }
        for(i=0;i<p;i++)
        	ptr_y[ i*lda2 + ty*lda3 + tx]+= (ptr_coeff[offset + i*lda2 + ty*lda3 +tx] * ptr_x[(i+zdisp)*lda2 + (ty+ydisp)*lda3 + (tx+xdisp)]);
        }
}


int SGCUDA_MatMult(PetscScalar* coeff, PetscScalar* x, PetscScalar* y, PetscInt *idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos)
{

//double tbegin3, tbegin4, tend3, tend4;
PetscInt size_coeff, size_xy, size_id; 
PetscScalar* d_coeff;
PetscScalar* d_x;
PetscScalar* d_y;
PetscInt *d_idx, *d_idy, *d_idz;
PetscInt i,j;

fprintf(stdout,"%d\t%d\t%d\t%d\n",m,n,p,nos);

//loading the coeff, x, y, idx, idy, idz to device memory

  unsigned int timer1 = 0;
  //cutilCheckError(cutCreateTimer(&timer1));
  //cutilCheckError(cutStartTimer(timer1));

  fprintf(stdout,"In SGCUDA_MatMult\n");
	
//tbegin3 = rtclock();
size_coeff=nos*m*n*p*sizeof(PetscScalar);
cudaMalloc((void**)&d_coeff,size_coeff);
cudaMemcpy(d_coeff, coeff, size_coeff, cudaMemcpyHostToDevice);

size_xy = m*n*p*sizeof(PetscScalar);
cudaMalloc((void**)&d_x,size_xy); 
cudaMemcpy(d_x, x, size_xy, cudaMemcpyHostToDevice);

cudaMalloc((void**)&d_y,size_xy); 
cudaMemcpy(d_y, y, size_xy, cudaMemcpyHostToDevice);

size_id = nos*sizeof(PetscInt);
cudaMalloc((void**)&d_idx,size_id); 
cudaMemcpy(d_idx, idx, size_id, cudaMemcpyHostToDevice);

cudaMalloc((void**)&d_idy,size_id); 
cudaMemcpy(d_idy, idy, size_id, cudaMemcpyHostToDevice);

cudaMalloc((void**)&d_idz,size_id); 
cudaMemcpy(d_idz, idz, size_id, cudaMemcpyHostToDevice);

//cutilCheckError(cutStopTimer(timer1));
// kernel Configuration

dim3 dimBlock(16,16);
dim3 dimGrid((size/16),(size/16));

//Cuda Printf
//cudaPrintfInit();

//tbegin4 = rtclock();
// create and start timer
    unsigned int timer = 0;
    //cutilCheckError(cutCreateTimer(&timer));
    //cutilCheckError(cutStartTimer(timer));

	MatMult_Kernel<<<dimGrid,dimBlock>>>(d_coeff, d_x, d_y, d_idx, d_idy, d_idz, m, n, p, nos);

   // check if kernel execution generated and error
    	//cutilCheckMsg("Kernel execution failed");

   // stop and destroy timer
    	//cutilCheckError(cutStopTimer(timer));
		
//tend4 = rtclock();
//Read y from the Device Memory

cudaMemcpy(y, d_y, size_xy, cudaMemcpyDeviceToHost); 
 
// double time_sec=cutGetTimerValue(timer)/1000;
// double time_sec1=cutGetTimerValue(timer1)/1000;
   
// printf("MFLOPS: GPU Structured Grid Matrix Mult kernel : %f; time(sec): %f\n",(2*stpoints*csr_size*csr_size*1.0e-6/time_sec),time_sec);
// printf("MFLOPS: GPU Structured Grid Matrix Mult kernel setup time(sec) : %f\n",time_sec1);
    
// cutilCheckError(cutDeleteTimer(timer));
// cutilCheckError(cutDeleteTimer(timer1));
// tend3 = rtclock();
// printf("MFLOPS: GPU Structured Grid Matrix Mult kernel with copy time : %f; time: %f\n",2*stpoints*csr_size*csr_size*1.0e-6/(tend3-tbegin3),tend3-tbegin3);
// printf("MFLOPS: GPU Structured Grid Matrix Mult kernel : %f; time: %f\n",2*stpoints*csr_size*csr_size*1.0e-6/(tend4-tbegin4),tend4-tbegin4);
  
// printf("\n");
// printf("Matrix cuda y\n");

  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++) 
    {
      printf("%.2f\n",y[i*size+j]);
    }
   printf("\n");
  }
 

//Free Device Memory
cudaFree(d_coeff);
cudaFree(d_x);
cudaFree(d_y);
cudaFree(d_idx);
cudaFree(d_idy);
cudaFree(d_idz);

return 0;
}



