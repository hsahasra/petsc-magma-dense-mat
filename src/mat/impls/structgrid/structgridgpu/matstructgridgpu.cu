
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
#include "cuPrintf.cu"

// ----------------------------------------------------------
// hardcodiing the shared memory size this should be set
// to give maximum performance, however should be
// replaced soon with a more flexable dynamically allocated
// shared memory scheme
// written by: dlowell ANL-MCS
// ----------------------------------------------------------
#define SHDSIZE 32


// -----------------------------------------------
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
       int numtiles;
       int tilesize;
};//820 bytes

__constant__ Stencilparams devparams;//device memory




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
     This function creates a datatype of structgridgpu. It first creates a
     structgrid datatype and overrides the matrix multiplication method.
     Author: Chekuri S. Choudary, RNET
*/

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatCreate_SeqSGGPU"
PetscErrorCode  MatCreate_SeqSGGPU(Mat B)
{

  //printf("Call to MatCreate_SeqSGGPU(Mat B)\n");
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr             = MatCreate_SeqSG(B);CHKERRQ(ierr);
  B->ops->mult     = MatMult_SeqSGGPU;

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


        //set up parameters for constant memory
        struct Stencilparams sparams;
        for(i=0;i<a->stpoints;i++){
            sparams.idx[i]=a->idx[i];
            sparams.idy[i]=a->idy[i];
            sparams.idz[i]=a->idz[i];
        }
        sparams.m=a->m;
        sparams.n=a->n;
        sparams.p=a->p;
        VecGetLocalSize(x,&sparams.vecsize_x);
        VecGetLocalSize(y,&sparams.vecsize_y);
        sparams.nos = a->stpoints;
        sparams.dof = a->dof;
        sparams.lda1=a->m*a->n*a->p;
        sparams.lda2=a->m*a->n;
        sparams.lda3=a->m;
        sparams.matsize=a->m*a->n*a->p*a->stpoints;

        /// Debugging block .....................................................
            //int xsize,ysize;
            //printf("Matrix A ::: m: %d, n: %d, p: %d, nos: %d dof: %d nz: %d\n",
            //    a->m,a->n,a->p,a->stpoints,a->dof,a->nz);
            //VecGetLocalSize(x,&xsize);
            //VecGetLocalSize(y,&ysize);
            //printf("Amat size: %d, Xvec size: %d, Yvec size: %d\n",sparams.matsize,xsize,ysize);
            //static PetscInt count = 1;// running count of function calls
            //printf("MatMult_SeqSGGPU(Mat mat, Vec x, Vec y): %d\n",count++);
        ///....................................................................

// Call to dlowell's version
        ierr = SGCUDA_MatMult_v2(v,xx,yy,sparams); CHKERRQ(ierr);

// Call to Jeswin's version
//        ierr = SGCUDA_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->stpoints);
//        CHKERRQ(ierr);

       	ierr = VecRestoreArray(x,&xx); CHKERRQ(ierr);
	ierr = VecRestoreArray(y,&yy); CHKERRQ(ierr);
	ierr = PetscLogFlops(2*a->nz*a->stpoints); CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
EXTERN_C_END







//-------------------------------------------------------------------------------
//   This function is the matrix vector multiplication kernel
//   structgridgpu datatype. This version uses shared memory for the write
//   back vector Y. Constant memory for reused constants and indices.
//   More offloading to registers might be possible as well.
//   written by: dlowell, ANL-MCS
//-------------------------------------------------------------------------------
__global__ void MatMul_Kernel_v2(double* A, double* X, double* Y){

   // indices for global accesses
   int btid = blockDim.x*blockIdx.x+threadIdx.x;

   //indices for local accesses
   int tid = threadIdx.x;

   int i,j;

   int nos  = devparams.nos; // set to register
   int lda1 = devparams.lda1;//  "   "
   int lda2 = devparams.lda2;//  "   "
   int lda3 = devparams.lda3;//  "   "
   int numtiles = devparams.numtiles;

   int Aindex;
   int Xindex;
   int index;

   __shared__ double Ys[SHDSIZE];

   for(j=0;j<nos;j++){
       for(i=0;i<numtiles;i++){
           index = (devparams.tilesize*i)+btid;
           Ys[tid]=0.;
           Aindex = j*lda1 + index;
           Xindex = (devparams.idz[j]*lda2 + devparams.idy[j]*lda3 + devparams.idx[j]) + index;
           __syncthreads();

           if (!((j==1 && ((index%devparams.m)+(index%devparams.n))==0)||
                (j==2 && index==0)||
                (j==3 && index%(devparams.n)==0)||
                (j==4 && index<devparams.m))){
                 if(Aindex<devparams.matsize && Xindex<devparams.vecsize_x) Ys[tid]=A[Aindex]*X[Xindex];
           }//end if

           if(index<devparams.vecsize_y) Y[index]+=Ys[tid];//global write back
       }//end i-for
       __syncthreads();
   }//end j-for

}




//------------------------------------------------------------------------------------
//   This function is the wrapper function which sets up the device memory, transfers
//   data to and from the device, and calls the kernel. Error checking is done at
//   each step. Timing stats are recorded using static vars.
//   written by: Daniel Lowell, ANL-MCS
//------------------------------------------------------------------------------------
PetscErrorCode SGCUDA_MatMult_v2(PetscScalar* A, PetscScalar* X, PetscScalar* Y, struct Stencilparams P){

        // vars for testing
        int i;
        static double cumtime=0.;//cummalitive call time
        static double cumkern=0.;//cummalitive kernel time
        static unsigned int kcalls=0;//number of kernel calls

        // using CUDA device timer
	float elapsedtime;
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	double cs,ce,temp;

	cudaError_t cudastatus0,cudastatus1,
	            cudastatus2,cudastatus3,
	            cudastatus4,cudastatus5,
           	    cudastatus6,cudastatus7;


        int matsize =P.matsize*sizeof(double);
        int vecsize_x = P.vecsize_x*sizeof(double);
        int vecsize_y = P.vecsize_y*sizeof(double);

// 	allocate GPU device memory
	cs=getclock();

        // Allocate and Memcpy Structured Matrix A
	double* devA;
	cudastatus0=cudaMalloc((void**)&devA,matsize);
	cudastatus1=cudaMemcpy(devA,A,matsize,cudaMemcpyHostToDevice);
	if(cudastatus0!=cudaSuccess|cudastatus1!=cudaSuccess){
	  printf("Error in devA memory allocation:\nstatus0: %s, status1: %s\n",
  			cudaGetErrorString(cudastatus0),
			cudaGetErrorString(cudastatus1));
	  if(devA) cudaFree(devA);
          PetscFunctionReturn(PETSC_ERR_MEM);
	}

        // Allocate and Memcpy Vector X
	double* devX;
	cudastatus2=cudaMalloc((void**)&devX,vecsize_x);
	cudastatus3=cudaMemcpy(devX,X,vecsize_x,cudaMemcpyHostToDevice);
	if(cudastatus2!=cudaSuccess|cudastatus3!=cudaSuccess){
	  printf("Error in devX memory allocation:\nstatus2: %s, status3: %s.\n",
  			cudaGetErrorString(cudastatus2),
			cudaGetErrorString(cudastatus3));
	  if(devA) cudaFree(devA);
	  if(devX) cudaFree(devX);
          PetscFunctionReturn(PETSC_ERR_MEM);
	}

        // Allocate and Memset(0.) Vector Y
	double* devY;
	cudastatus4=cudaMalloc((void**)&devY,vecsize_y);
	cudastatus5=cudaMemset(devY,0.0,vecsize_y);
	if(cudastatus4!=cudaSuccess|cudastatus5!=cudaSuccess){
	  printf("Error in devY memory allocation:\nstatus4: %s, status5: %s\n",
	  		cudaGetErrorString(cudastatus4),
			cudaGetErrorString(cudastatus5));
	  if(devA) cudaFree(devA);
	  if(devY) cudaFree(devY);
          if(devX) cudaFree(devX);
          PetscFunctionReturn(PETSC_ERR_MEM);
	}

        //Set up blocks and thread numbers
        P.tilesize = floor(49152.0/(float)sizeof(double));///max shared elements
        if(P.vecsize_y < P.tilesize) P.tilesize = P.vecsize_y;
        P.numtiles = ceil((float)P.vecsize_y/P.tilesize);

        int blocks = ceil((float)P.tilesize/(float)SHDSIZE);//number of blocks
        int threads = SHDSIZE;//number of threads per block

        dim3 dimGrid(blocks,1);
	dim3 dimBlock(threads,1);

        static unsigned char dbgflag = 1;
        if(dbgflag){
           printf("blocks: %d, threads: %d, Tnum: %d, Tsize: %d,  P.vecsize_y: %d\n",
                        blocks,threads,P.numtiles,P.tilesize,P.vecsize_y);
                dbgflag=0;
        }

        // update constant memory with structured grid parameters
	cudastatus6=cudaMemcpyToSymbol("devparams",&P,sizeof(Stencilparams));
	if(cudastatus6!=cudaSuccess){
	  printf("Error in symbol copy: status6: %s.\n",
	  		cudaGetErrorString(cudastatus6));
	  if(devA) cudaFree(devA);
	  if(devY) cudaFree(devY);
	  if(devX) cudaFree(devX);
          PetscFunctionReturn(PETSC_ERR_MEM);
	}


	ce=getclock();//end setup timer
	temp=ce-cs;

        cudaPrintfInit();//start cuda printf environ.
	cudaEventRecord(start,0);//begin recording kernel
	MatMul_Kernel_v2<<<dimGrid,dimBlock>>>(devA,devX,devY);
        checkCUDAError("CUDA Kernel launch...");//check for failure
        cudaPrintfDisplay(stdout, true);//choose output
        cudaPrintfEnd();//kill cuda printf environ
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop); // event barrier
	cudaEventElapsedTime(&elapsedtime,start,stop);
        cudaEventDestroy(start);
	cudaEventDestroy(stop);


        // Copy back Vector Y from Kernel
	cs=getclock();
	cudastatus7=cudaMemcpy(Y,devY,vecsize_y,cudaMemcpyDeviceToHost);
	if(cudastatus7!=cudaSuccess){
          printf("Error on copy back Y, kernel status: %s\nExiting...\n\n",cudaGetErrorString(cudastatus7));
	  if(devA) cudaFree(devA);
	  if(devY) cudaFree(devY);
	  if(devX) cudaFree(devX);
          PetscFunctionReturn(PETSC_ERR_MEM);
        }
        //for(i=0;i<P.lda1;i++)printf("Y[%d]: %lf\n",i,Y[i]);//for verification

        //Free device memory
	if(devA) cudaFree(devA);
	if(devY) cudaFree(devY);
	if(devX) cudaFree(devX);

	ce=getclock();
	temp+=1000.0*(ce-cs);
        cumkern+=(elapsedtime);
        cumtime+=(elapsedtime)+temp;
        kcalls++;
        printf("Kernel call #: %d\n",kcalls);
        printf("Cumilative function call time: %lf msec.\n", cumtime);
        printf("Cumilative kernel time: %lf msec.\n", cumkern);
        printf("Current call: setup+copyback: %lf msec., kernel time: %f msec.\n",temp,elapsedtime);
        printf("-------------------------------------------------------------------------\n\n");

        PetscFunctionReturn(0);
}









































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

//fprintf(stdout,"%d\t%d\t%d\t%d\n",m,n,p,nos);

//loading the coeff, x, y, idx, idy, idz to device memory

  unsigned int timer1 = 0;
  //cutilCheckError(cutCreateTimer(&timer1));
  //cutilCheckError(cutStartTimer(timer1));

  //  fprintf(stdout,"In SGCUDA_MatMult\n");
	
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
if (m > 16){
dim3 dimBlock(16,16);
dim3 dimGrid((m/16),(n/16));

    // cutilCheckError(cutCreateTimer(&timer));
    // cutilCheckError(cutStartTimer(timer));

MatMult_Kernel<<<dimGrid,dimBlock>>>(d_coeff, d_x, d_y, d_idx, d_idy, d_idz, m, n, p, nos);

}
else
{
dim3 dimBlock(m,n);
dim3 dimGrid(1,1);
   
    // cutilCheckError(cutCreateTimer(&timer));
    // cutilCheckError(cutStartTimer(timer));

MatMult_Kernel<<<dimGrid,dimBlock>>>(d_coeff, d_x, d_y, d_idx, d_idy, d_idz, m, n, p, nos);


}

//Cuda Printf
//cudaPrintfInit();

//tbegin4 = rtclock();
// create and start timer
    //unsigned int timer = 0;
    //cutilCheckError(cutCreateTimer(&timer));
    //cutilCheckError(cutStartTimer(timer));

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

//  for(i=0;i<m;i++)
//  {
//    for(j=0;j<n;j++) 
//    {
//      printf("%.2f\n",y[i*n+j]);
//    }
//   printf("\n");
//  }


//Free Device Memory
cudaFree(d_coeff);
cudaFree(d_x);
cudaFree(d_y);
cudaFree(d_idx);
cudaFree(d_idy);
cudaFree(d_idz);

return 0;
}




