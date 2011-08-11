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

#define _DBGFLAG 1

// ----------------------------------------------------------
// hardcodiing the shared memory size this should be set
// to give maximum performance, however should be
// replaced soon with a more flexable dynamically allocated
// shared memory scheme
// written by: dlowell ANL-MCS
// ----------------------------------------------------------
#define SHDSIZE 4


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
  printf("Call to MatDestroy_SeqSGGPU(Mat B)\n");
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
            /*int xsize,ysize;
            //printf("Matrix A ::: m: %d, n: %d, p: %d, nos: %d dof: %d nz: %d\n",
            //    a->m,a->n,a->p,a->stpoints,a->dof,a->nz);
            //VecGetLocalSize(x,&xsize);
            //VecGetLocalSize(y,&ysize);
            //printf("Amat size: %d, Xvec size: %d, Yvec size: %d\n",sparams.matsize,xsize,ysize);
            */
            //static PetscInt count = 1;// running count of function calls
            //printf("MatMult_SeqSGGPU(Mat mat, Vec x, Vec y): %d\n",count++);
        ///....................................................................


// Call to dlowell's version
      ierr = SGCUDA_MatMult_v2(v,xx,yy,sparams,&(mat->valid_GPU_matrix));
	// CHKERRQ(ierr);

// Call to Jeswin's version
     //   ierr = SGCUDA_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->stpoints,&(mat->valid_GPU_matrix));CHKERRQ(ierr);

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

   //indices for local accesses
   int tbtx, tbty, tbtz;
   int ix,iy,iz,j;

   int nos  = devparams.nos; // set to register
   int lda1 = devparams.lda1;//  "   "
   int lda2 = devparams.lda2;//  "   "
   int lda3 = devparams.lda3;//  "   "
   int tilex = devparams.tile_x*devparams.tsizex;
   int tiley = devparams.tile_y*devparams.tsizey;
   int tilez = devparams.tile_z*devparams.tsizez;

   int Aindex;
   int Xindex;
   int index;


   //Min shared mem byte count: 2*(gridDim.x*gridDim.y*gridDim.z)*SHDSIZE^3*8 byte double
   __shared__ double Ys[SHDSIZE][SHDSIZE][SHDSIZE];
   __shared__ double As[SHDSIZE][SHDSIZE][SHDSIZE];

//------------------------------------------------------------------------

   for(iz=0;iz<tilez;iz+=devparams.tsizez){// tiles Z loop
        tbtz = blockDim.z*blockIdx.z+threadIdx.z + iz;

   for(iy=0;iy<tiley;iy+=devparams.tsizey){// tiles Y loop
        tbty = blockDim.y*blockIdx.y+threadIdx.y + iy;

   for(ix=0;ix<tilex;ix+=devparams.tsizex){// tiles X loop
        tbtx = blockDim.x*blockIdx.x+threadIdx.x + ix;
        __syncthreads();

        //cuPrintf("tid xyz: %d, %d, %d\n",tbtx,tbty,tbtz);
        //initialize current return-tile
        Ys[threadIdx.z][threadIdx.y][threadIdx.x]=0.;

        //adjusted index for global access
        index = tbtz*lda2 + tbty*lda3 + tbtx;

        __syncthreads();

//......STENCIL...........................................
        for(j=0;j<nos;j++){//loop over stencil pattern

           Aindex=j*lda1+index;//set up Aindex and read from global A into As tile
           if(Aindex<devparams.matsize) As[threadIdx.z][threadIdx.y][threadIdx.x]=A[Aindex];//needs to be coalesced
           else As[threadIdx.z][threadIdx.y][threadIdx.x]=0.;

           __syncthreads();

           if(!((j==1 && tbtx==devparams.m-1 && tbty==devparams.n-1)||
                (j==2 && tbtx==0 && tbty==0)|| (j==3 && tbty==devparams.n-1)||
                (j==4 && tbty==0))){

              //set up Xindex for element-wise operation using stencil pattern
              Xindex=(devparams.idz[j]*lda2 + devparams.idy[j]*lda3 + devparams.idx[j]) + index;
              if(Xindex<devparams.vecsize_x)
                 Ys[threadIdx.z][threadIdx.y][threadIdx.x]+=As[threadIdx.z][threadIdx.y][threadIdx.x]*X[Xindex];

           }//end if
           __syncthreads();

        }//end j-for

        __syncthreads();
        if(index<devparams.vecsize_y) Y[index]=Ys[threadIdx.z][threadIdx.y][threadIdx.x];//global write back

   }//end ix-for
   }//end iy-for
   }//end iz-for
}//end kernel_v2






//   int Xoffset;
              //Xoffset=(devparams.idz[j]*lda2+devparams.idy[j]*lda3+devparams.idx[j]);
              //  cuPrintf("Xindex: %d, Xoffset: %d, Xsize: %d idx: %d, idy: %d, idz: %d, j: %d nos: %d\n",Xindex,Xoffset,devparams.vecsize_x,devparams.idx[j],devparams.idy[j],devparams.idz[j],j,nos);


//------------------------------------------------------------------------------------
//   This function is the wrapper function which sets up the device memory, transfers
//   data to and from the device, and calls the kernel. Error checking is done at 
//   each step. Timing stats are recorded using static vars.
//   written by: Daniel Lowell, ANL-MCS
//------------------------------------------------------------------------------------
PetscErrorCode SGCUDA_MatMult_v2(PetscScalar* A, PetscScalar* X, 
PetscScalar* Y, struct Stencilparams P, PetscCUSPFlag* fp){

        // vars for testing
        int i;
        static double cumktime=0.;//cummalitive kernel time
        static double cumtime=0.;//cummalitive call time
        static unsigned int kcalls=0;//number of kernel calls
	double cs,ce,temp;
	float elapsedtime;        // using CUDA device timer
	cudaEvent_t start,stop;

        static unsigned char allocflag = 1;
        static double maxshared;
        static int bx,by,bz;//number of blocks in 3-D
        static int tx,ty,tz;//number of threads ber block in 3-D
        static int maxblocks_xy;
        static int maxblocks_z;
        static dim3 dimGrid;
	static dim3 dimBlock;
        static int xytile;

	cudaError_t cudastatus0,cudastatus1,
	            cudastatus2,cudastatus3,
	            cudastatus4,cudastatus5,
           	    cudastatus6,cudastatus7;


        //size in bytes to be allocated onto device
        int matsize =P.matsize*sizeof(double);
        int vecsize_x = P.vecsize_x*sizeof(double);
        int vecsize_y = P.vecsize_y*sizeof(double);


        if(_DBGFLAG){//create CUDA events for timer
           cudaEventCreate(&start);
           cudaEventCreate(&stop);
        }


	if(_DBGFLAG) cs=getclock();

        //Allocate and Memcpy Structured Matrix A
	//The matrix remains the same throughout one iteration
        //of the linear solver. The following uses a flag
        //defined in the base class to check the status of the
        //matrix. The matrix is copied to the GPU only if
        //it has been changed on the CPU side
        //This feature added by Chekuri S. Choudary

        if ((*fp == PETSC_CUSP_UNALLOCATED) || (*fp == PETSC_CUSP_CPU)){
		if (*fp == PETSC_CUSP_UNALLOCATED){
	   	   cudastatus0=cudaMalloc((void**)&devA,matsize);
	   	   if(cudastatus0!=cudaSuccess){
                        printf("Error in devA memory allocation:\nstatus0: %s\n",
  			cudaGetErrorString(cudastatus0));
          	        PetscFunctionReturn(PETSC_ERR_MEM);
		   }
		}

           	cudastatus1=cudaMemcpy(devA,A,matsize,cudaMemcpyHostToDevice);
	   	if(cudastatus1!=cudaSuccess){
		  if(devA) cudaFree(devA);
		  printf("Error in devA memory copying:\nstatus1: %s\n",
  			cudaGetErrorString(cudastatus1));
          	  PetscFunctionReturn(PETSC_ERR_MEM);
		}

	       *fp = PETSC_CUSP_BOTH;
	}


        //Allocate device memory for X and Y, and shape grid and blocks
        if(allocflag){
                cudastatus2=cudaMalloc((void**)&devX,vecsize_x);//allocate X on device
	        if(cudastatus2!=cudaSuccess){
                        printf("Error in devX memory allocation: %s\n",cudaGetErrorString(cudastatus2));
	                if(devA) cudaFree(devA);
                        PetscFunctionReturn(PETSC_ERR_MEM);
                }

                cudastatus3=cudaMalloc((void**)&devY,vecsize_y);//allocate Y on device
                if(cudastatus3!=cudaSuccess){
                        printf("Error in devY memory allocation: %s\n",cudaGetErrorString(cudastatus3));
	                if(devA) cudaFree(devA);
	                if(devX) cudaFree(devX);
                        PetscFunctionReturn(PETSC_ERR_MEM);
	        }



                //Set up blocks and thread numbers
                maxshared = 49152.0/(double)(2.0*sizeof(double));
                if(P.p==1){
                    xytile = pow(maxshared,0.5);//square blocks
                    maxblocks_z = 1;
                }else{
                    temp=maxshared/P.p;//lop off z
                    xytile = pow(temp,0.5);//xyblocks
                    maxblocks_z=ceil((float)SHDSIZE/(float)P.p);
                }
                maxblocks_xy = xytile/SHDSIZE;


                //Set up blocks and thread numbers for columns
                if(P.m <= SHDSIZE){
                       tx = P.m;
                       bx = 1;
                       P.tile_x = 1;
                       P.tsizex=1;
                }else{
                       tx = SHDSIZE;
                       bx = ceil((float)P.m/(float)SHDSIZE);//create enough blocks
                       if(bx>maxblocks_xy){                    //too many blocks created
                          bx = maxblocks_xy;                   //set to max number of blocks allowed
                          P.tile_x=ceil((float)P.m/(float)(bx*SHDSIZE));//number of tiles
                          P.tsizex=bx*SHDSIZE;              //tilesize is block-thread coverage
                       }else{
                          P.tile_x=1;
                          P.tsizex=1;
                       }
                }

                //Set up blocks and thread numbers for rows
                if(P.n <= SHDSIZE){
                       ty = P.n;
                       by = 1;
                       P.tile_y = 1;
                       P.tsizey=1;
                }else{
                       ty = SHDSIZE;
                       by = ceil((float)P.n/(float)SHDSIZE);
                       if(by > maxblocks_xy){
                          by = maxblocks_xy;
                          P.tile_y=ceil((float)P.n/(float)(by*SHDSIZE));
                          P.tsizey=by*SHDSIZE;
                       }else{
                          P.tile_y=1;
                          P.tsizey=1;
                       }
                }

                //Set up blocks and thread numbers for z
                if(P.p <= SHDSIZE){
                       tz = P.p;
                       bz = 1;
                       P.tile_z = 1;
                       P.tsizez=1;
                }else{
                       tz = SHDSIZE;
                       bz = ceil((float)P.p/(float)SHDSIZE);
                       if(bz > maxblocks_z){
                          bz = maxblocks_z;
                          P.tile_z=ceil((float)P.p/(float)(bz*SHDSIZE));
                          P.tsizez=bz*SHDSIZE;
                       }else{
                          P.tile_z=1;
                          P.tsizez=1;
                       }
                }

                //set grid shape
                dimGrid.x = bx;
                dimGrid.y = by;
                dimGrid.z = bz;

                //set block shape
                dimBlock.x = tx;
                dimBlock.y = ty;
                dimBlock.z = tz;

                // update constant memory with structured grid parameters
	        cudastatus6=cudaMemcpyToSymbol("devparams",&P,sizeof(Stencilparams));
	        if(cudastatus6!=cudaSuccess){
                        printf("Error in symbol copy to device: %s.\n",cudaGetErrorString(cudastatus6));
                        if(devA) cudaFree(devA);
	                if(devY) cudaFree(devY);
	                if(devX) cudaFree(devX);
                        PetscFunctionReturn(PETSC_ERR_MEM);
	        }

                //toggle off allocation flag
                allocflag = 0;

        }//end allocflag-if

        //grid and block shape & device config. debugging.....................................
        unsigned int sharebytes = 2*tx*ty*tz*bx*by*bz*sizeof(double);
        static unsigned char dbgflag = 1;
        if(dbgflag){
           printf("(m, n, p, nos): (%d, %d, %d, %d)\n",P.m,P.n,P.p,P.nos);
           printf("MAXBLOCKS: %d maxshared: %lf\n",maxblocks_xy+maxblocks_z,maxshared);
           printf("blocks: (%d, %d, %d), threads per block: %d\n", bx,by,bz,tx*ty*tz );
           printf("Shared elements occupied: %0.3f SharedOccupied in Bytes: %d\n",sharebytes/49152.0,sharebytes);
           printf("Blocks x,y,z: (%d, %d, %d)\n",bx,by,bz);
           printf("Blocks*ThreadsPer size x,y,z: (%d, %d, %d)\n",bx*tx,by*ty,bz*tz);
           printf("Tiles (x,y,z): (%d, %d, %d)\n",P.tile_x,P.tile_y,P.tile_z);
           printf("Tile Size (x,y,z): (%d, %d, %d)\n",P.tsizex,P.tsizey,P.tsizez);
           dbgflag=0;
        }
        //...End config. debug section.................................................................




       //copy over values of X to device memory
	cudastatus4=cudaMemcpy(devX,X,vecsize_x,cudaMemcpyHostToDevice);
	if(cudastatus4!=cudaSuccess){
                printf("Error in devX memory copy to device: status: %s\n",cudaGetErrorString(cudastatus4));
	        if(devA) cudaFree(devA);
	        if(devX) cudaFree(devX);
	        if(devY) cudaFree(devY);
                PetscFunctionReturn(PETSC_ERR_MEM);
	}

/*//probably an unnecessary step.
        // memset to 0. Vector Y on device
	cudastatus5=cudaMemset(devY,0.0,vecsize_y);
	if(cudastatus5!=cudaSuccess){
                printf("Error in devY memset to device: %s\n",cudaGetErrorString(cudastatus5));
	        if(devA) cudaFree(devA);
                if(devY) cudaFree(devY);
                if(devX) cudaFree(devX);
                PetscFunctionReturn(PETSC_ERR_MEM);
	}
*/

        //toggle timer and debug settings
        if(_DBGFLAG){
                ce=getclock();//end setup timer
	        temp=ce-cs;
                cudaPrintfInit();//start cuda printf environ.
	        cudaEventRecord(start,0);//begin recording kernel
        }

        //Launch the kernel..........................................
	MatMul_Kernel_v2<<<dimGrid,dimBlock>>>(devA,devX,devY);
        checkCUDAError("CUDA Kernel launch...");//check for failure
        //...........................................................

        //toggle timer and debug settings
        if(_DBGFLAG){
                cudaPrintfDisplay(stdout, true);//choose output
                cudaPrintfEnd();//kill cuda printf environ
	        cudaEventRecord(stop,0);
	        cudaEventSynchronize(stop); // event barrier
	        cudaEventElapsedTime(&elapsedtime,start,stop);
                cudaEventDestroy(start);
	        cudaEventDestroy(stop);
        }

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



        if(_DBGFLAG){
          //for(i=0;i<P.lda1;i++)printf("Y[%d]: %lf\n",i,Y[i]);//for verification
	  ce=getclock();
	  temp+=ce-cs;
          cumktime+=(elapsedtime/1000);
          cumtime+=(elapsedtime/1000)+temp;
          kcalls++;
          printf("Kernel call #: %d\n",kcalls);
          printf("setup+copyback: %f sec.\nelapsed time: %f sec.\ntotal call time: %f sec.\n",
                  temp,elapsedtime/1000,(elapsedtime/1000)+temp);
          printf("Cum. kernel time: %lf sec.\n", cumktime);
          printf("Cum. call time (with setup): %lf sec.\n", cumtime);
          printf(".........................................\n\n");
        }//end _DBGFLAG-if

        PetscFunctionReturn(0);
}





/*  -------------------------------------------------------------------- 
     The following is a CUDA kernel for matrix vector multiplication on 
     the GPU. The matrix is in a custom layout that facilitates better 
     memory accesses and vectorization. 
     Author: Chekuri S. Choudary, RNET
*/
/*  __global__ void MatMult_Kernel(PetscScalar * ptr_coeff, PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt *idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos)
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
 */ 
#define BLOCKWIDTH 16
 
 __global__ void MatMult_Kernel(PetscScalar * ptr_coeff, PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt* idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos)
{

int tx= blockDim.x * blockIdx.x + threadIdx.x;
int ty= blockDim.y * blockIdx.y + threadIdx.y;
int l,i;
int xdisp,ydisp,zdisp,offset;
int lda1=m*n*p,lda2=m*n,lda3=m;
__shared__ double y_sm[256];
__shared__ double x_sm[324];

// copying a Tile from Y into the shared Memory
y_sm[threadIdx.y*BLOCKWIDTH + threadIdx.x]=0;

// Copying a tile x into the shared Memory with 2 steps.

// Copying without the Ghost Cells  
x_sm[(threadIdx.y+1)*(BLOCKWIDTH+2) + (threadIdx.x+1)]=ptr_x[ty*lda3 + tx];


// Copying the Ghost Cells
if (tx!=0 || ty!=0 || tx != n-1 || ty != m-1)
{
	if (threadIdx.x==0)
	x_sm[(threadIdx.y+1)*(BLOCKWIDTH+2) + threadIdx.x]=ptr_x[ty*lda3 + tx-1];

	if (threadIdx.y==0)
	x_sm[(threadIdx.y)* (BLOCKWIDTH+2) + threadIdx.x+1]=ptr_x[(ty-1)*lda3 + tx];

	if (threadIdx.x==BLOCKWIDTH-1)
	x_sm[(threadIdx.y+1)*(BLOCKWIDTH+2) + threadIdx.x + 2]=ptr_x[ty*lda3 + tx + 1];

	if (threadIdx.y==BLOCKWIDTH-1)
	x_sm[(threadIdx.y+2)*(BLOCKWIDTH+2) + threadIdx.x +1]=ptr_x[(ty+1)*lda3 + tx];
}
__syncthreads();

// if (tx==2 && ty==2)
// {
// cuPrintf("\nPrinting the X from Shared Memory \n ");

// for (int j=0;j<324;j++)
// {
// if(j % 16 ==0)
// {
// cuPrintf("\n");
// }
// cuPrintf("%f  ",  x_sm[j]);
// }
// }
// MATMUL
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
	y_sm[threadIdx.y*BLOCKWIDTH + threadIdx.x]+= (ptr_coeff[offset + i*lda2 + ty*lda3 +tx] * x_sm[(i+zdisp)*lda2 + (threadIdx.y+ydisp +1)*(BLOCKWIDTH+2) + (threadIdx.x+xdisp+1)]); //forgetting Z currently.. I have to Fix it.
	}
	// removing i tempararily
	ptr_y[ty*lda3 + tx]= y_sm[threadIdx.y*BLOCKWIDTH + threadIdx.x];
}
  
int SGCUDA_MatMult(PetscScalar* coeff, PetscScalar* x, PetscScalar* y, PetscInt *idx, PetscInt* idy, 
PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos, PetscCUSPFlag* fp)
{

double tbegin1, tbegin2, tend1, tend2;
static PetscInt size_coeff; 
double tsetup,tkernel;
static unsigned int kcalls=0;
PetscInt size_xy, size_id; 
static double temp=0;
PetscScalar* d_x;
PetscScalar* d_y;
PetscInt *d_idx, *d_idy, *d_idz;

  //unsigned int timer1 = 0;
  //cutilCheckError(cutCreateTimer(&timer1));
  //cutilCheckError(cutStartTimer(timer1));

  //  fprintf(stdout,"In SGCUDA_MatMult\n");
	
      if(_DBGFLAG) tbegin1=getclock();
	  
	  if ((*fp == PETSC_CUSP_UNALLOCATED) ||
	  (*fp == PETSC_CUSP_CPU) )
	{
		if (*fp == PETSC_CUSP_UNALLOCATED)
		{
		size_coeff=nos*m*n*p*sizeof(PetscScalar);	
		cudaMalloc((void**)&d_coeff,size_coeff);
	
	   	//cudastatus0=cudaMalloc((void**)&devA,matsize);
	   	//if(cudastatus0!=cudaSuccess)
		//	{
		//  printf("Error in devA memory allocation:\nstatus0: %s\n",
  		//	cudaGetErrorString(cudastatus0));
          	//  PetscFunctionReturn(PETSC_ERR_MEM);
		//	}
		}
	
		cudaMemcpy(d_coeff, coeff, size_coeff, cudaMemcpyHostToDevice);
	
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

if(_DBGFLAG) 
	{
		tend1=getclock();
		tsetup=tend1-tbegin1;
		tbegin2=getclock();
	}

//cutilCheckError(cutStopTimer(timer1));
// kernel Configuration
if (m > BLOCKWIDTH){
dim3 dimBlock(BLOCKWIDTH,BLOCKWIDTH);
dim3 dimGrid((int)ceil((float)m/(float)BLOCKWIDTH),((int)ceil((float)n/(float)BLOCKWIDTH)));

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

	if(_DBGFLAG) 
	{
		cudaThreadSynchronize();
		tend2=getclock();
		tkernel=tend2-tbegin2;
	}

	
	
   // check if kernel execution generated and error
    	//cutilCheckMsg("Kernel execution failed");

   // stop and destroy timer
    	//cutilCheckError(cutStopTimer(timer));
		
//Read y from the Device Memory

cudaMemcpy(y, d_y, size_xy, cudaMemcpyDeviceToHost); 
 
// double time_sec=cutGetTimerValue(timer)/1000;
// double time_sec1=cutGetTimerValue(timer1)/1000;
   
// printf("MFLOPS: GPU Structured Grid Matrix Mult kernel : %f; time(sec): %f\n",(2*stpoints*csr_size*csr_size*1.0e-6/time_sec),time_sec);
// printf("MFLOPS: GPU Structured Grid Matrix Mult kernel setup time(sec) : %f\n",time_sec1);
    
// cutilCheckError(cutDeleteTimer(timer));
// cutilCheckError(cutDeleteTimer(timer1));
if(_DBGFLAG)
{
temp+=tkernel;
	if (kcalls==0)
			{
			printf("\n Structured Grid MatrixMul Kernel Permormance for m *%d* and n size *%d* \n",m,n);
			}
	if (kcalls==1000)
		{
		printf("\ncopy time (sec) : %f\n",tsetup);
		printf("Kernel time (sec): %f\n",tkernel);
		printf("Performance in Megaflops with for %dth Kernel call\n",kcalls);
		printf("Performance in Megaflops with copy time = %f\n",(2*nos*n*m*p*1.0e-6)/(tsetup+tkernel));
		printf("Performance in Megaflops without copy time = %f\n",(2*nos*n*m*p*1.0e-6)/tkernel);
		printf("Culmative Performance in Megaflops for *%d* calls without copy time = %f\n",kcalls,(2*nos*n*m*p*1.0e-6)/(temp/(kcalls+1)));
		}
}
kcalls++;
//Free Device Memory
//cudaFree(d_coeff);
cudaFree(d_x);
cudaFree(d_y);
cudaFree(d_idx);
cudaFree(d_idy);
cudaFree(d_idz);

return 0;
}




