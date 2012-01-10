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

//block size is 1x256. 
#define BLOCKWIDTH_X 4		
#define BLOCKWIDTH_Y 1   

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
	sparams.lda3=a->m*a->dof;
	sparams.lda2=sparams.lda3*a->n;
        sparams.lda1=sparams.lda2*a->p;
        sparams.matsize=sparams.lda1*a->stpoints;

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
    //ierr = SGCUDA_MatMult_v2(v,xx,yy,sparams,&(mat->valid_GPU_matrix));
	// CHKERRQ(ierr);

// Call to Jeswin's version
    ierr = SGCUDA_MatMult(v,xx,yy,a->idx,a->idy,a->idz,a->m,a->n,a->p,a->stpoints,&(mat->valid_GPU_matrix),a->dof);CHKERRQ(ierr);

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
   //int offset = 0;

   int Aindex;
   int Xindex;
   int index;
   int offset;

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

        //cuPrintf("tid xyz: %d, %d, %d\n",tbtx,tbty,tbtz);
        //initialize current return-tile
        Ys[threadIdx.z][threadIdx.y][threadIdx.x]=0.;

        //adjusted index for global access
        index = tbtz*lda2 + tbty*lda3 + tbtx;

//......STENCIL...........................................
        for(j=0;j<nos;j++){//loop over stencil pattern
	   offset= j*lda1;
           Aindex=offset+index;//set up Aindex and read from global A into As tile
           if(Aindex<devparams.matsize) As[threadIdx.z][threadIdx.y][threadIdx.x]=A[Aindex];//needs to be coalesced
           else As[threadIdx.z][threadIdx.y][threadIdx.x]=0.;

           __syncthreads();

           //set up Xindex for element-wise operation using stencil pattern
           Xindex=(devparams.idz[j]*lda2 + devparams.idy[j]*lda3 + devparams.idx[j]) + index;
           if(Xindex<devparams.vecsize_x && Xindex>=0){
              Ys[threadIdx.z][threadIdx.y][threadIdx.x]+=As[threadIdx.z][threadIdx.y][threadIdx.x]*X[Xindex];
           }

        }//end j-for

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

//probably an unnecessary step.
        // memset to 0. Vector Y on device
	cudastatus5=cudaMemset(devY,0.0,vecsize_y);
	if(cudastatus5!=cudaSuccess){
                printf("Error in devY memset to device: %s\n",cudaGetErrorString(cudastatus5));
	        if(devA) cudaFree(devA);
                if(devY) cudaFree(devY);
                if(devX) cudaFree(devX);
                PetscFunctionReturn(PETSC_ERR_MEM);
	}


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

	//for(i=0;i<P.lda1;i++)printf("Y[%d]: %lf\n",i,Y[i]);//for verification


        if(_DBGFLAG){
          for(i=0;i<P.lda1;i++)printf("Y[%d]: %lf\n",i,Y[i]);//for verification
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


 //Version with Shared memory for X only supports rectangular tiles.
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
//  Below functions are SPMV kernel functions where x through the Texture Memory, offsets are accesed
//	through the Shared Memory, Y is accessed per thread from registers. Coeff accesses are from the global Memory 
//	but they are coalesced.  	   
//------------------------------------------------------------------------------------ 
#define stpoints 5 // I have to fix this 

__global__ void MatMul_Kernel_tex_1_DOF(PetscScalar * ptr_coeff, PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt *idx, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos,PetscInt DOF)
	{
		
		__shared__ float idx_sm[stpoints];
		
		int tx= blockDim.x * blockIdx.x + threadIdx.x;
		int l,offset;
		int lda1=m*n*p*DOF,lda2=m*p*DOF;  //lda3=m*DOF
		PetscInt Index;
		PetscScalar y_reg=0;
		
		if (threadIdx.x < stpoints)
			{
			idx_sm[threadIdx.x]=idx[threadIdx.x];
			}
		int reg2=blockIdx.y*lda2+tx;
		
		//Iterating through the Diagonals
		for (l=0;l<stpoints;l++)
			{
			Index =reg2 + idx_sm[l];
				
			if (Index >= 0 && Index <lda1)
				{
				offset = l*lda1;
				y_reg+= ptr_coeff[offset + reg2] * fetch_double(tex_x_double,Index);
				
			
				  if (threadIdx.y==0){
							cuPrintf("l= %d ptr_coeff= %f X= %f Index =%d y_sm=%f \n",l,ptr_coeff[offset + reg2],tex1Dfetch(tex_x_double,Index),Index, y_reg);
						}  
					
				}
			}
							
				
				ptr_y[reg2]= y_reg;
	}

	
__global__ void MatMul_Kernel_tex(double * ptr_coeff, double* ptr_x, double* ptr_y, PetscInt *idx, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos,PetscInt DOF)
	{
		
		__shared__ float idx_sm[stpoints];
		
		int tx= blockDim.x * blockIdx.x + threadIdx.x;
		int l,i,offset;
		int lda1=m*n*p*DOF,lda2=m*p*DOF; //lda3=m*DOF
		PetscInt X_Index,Index;
		double y_reg=0;
		int BAND_SIZE=(DOF-1)*2+1;
		
		if (threadIdx.x < stpoints)
			{
			idx_sm[threadIdx.x]=idx[threadIdx.x];
			}
		
		
		int reg2=blockIdx.y*lda2+tx;
		
		//Iterating through the Diagonals
		for (l=0;l<stpoints;l++)
			{
			X_Index =reg2 + idx_sm[l];
				
			if (X_Index >= 0 && X_Index <lda1)
				{
				for (i=0;i<BAND_SIZE;i++)
					{
					offset = (l*BAND_SIZE+i)*lda1;
									
						if (i > DOF-1)
							{
							Index =X_Index-(i-(DOF-1));
							if (Index < 0)
								{
								continue;
								}
							else{
								y_reg+= ptr_coeff[offset + reg2] * fetch_double(tex_x_double,Index) ;
								}
							}	
						else {
							Index=X_Index+i;
							
							y_reg+= ptr_coeff[offset + reg2] * fetch_double(tex_x_double,Index);
							}
			
				  /* if (threadIdx.y==0){
							cuPrintf("l= %d ptr_coeff= %f X= %f X_Index =%d Index =%d y_sm=%f \n",l,ptr_coeff[offset + reg2],tex1Dfetch(tex_x_double,Index),X_Index,Index, y_reg);
						} */  
					}
				}
			}
							
				
				ptr_y[reg2]= y_reg;
	}
	

//------------------------------------------------------------------------------------
//   The function is a wrapper function which sets up the device memory, transfers
//   data to and from the device, and calls the MatMult kernel. 
//------------------------------------------------------------------------------------ 
 
 int SGCUDA_MatMult(PetscScalar* coeff, PetscScalar* x, PetscScalar* y, PetscInt *idx, PetscInt* idy, 
PetscInt* idz, PetscInt m, PetscInt n,PetscInt p, PetscInt nos, PetscCUSPFlag* fp,PetscInt DOF)
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
int BLOCK_SIZE;
int cons=m*DOF;
int cons1=m*n*DOF;

 //Reducing to a Single offset instead of using three offsets int the x,y and z direction.  
  
  idx[0]=0;
  idx[1]=DOF;
  idx[2]=-DOF;
  idx[3]=cons;
  idx[4]=-cons;
  if(nos==7)
    {
     idx[5]=cons1;
     idx[6]=-cons1;       
    }

	
      if(_DBGFLAG) tbegin1=getclock();
	  
	  if ((*fp == PETSC_CUSP_UNALLOCATED) ||
	  (*fp == PETSC_CUSP_CPU) )
	{
		if (*fp == PETSC_CUSP_UNALLOCATED)
		{
		size_coeff=nos*m*n*p*DOF*sizeof(PetscScalar);	
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


size_xy = m*n*p*DOF*sizeof(PetscScalar);
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

//Binding X to the texture Memory
cudaBindTexture(0, tex_x_double, d_x, size_xy);

if (_DBGFLAG){
cudaPrintfInit();
}

if(_DBGFLAG) 
	{
		tend1=getclock();
		tsetup=tend1-tbegin1;
		tbegin2=getclock();
	}

// Kernel Setup and Configurations
	
	if ((DOF%2)!=0 || (DOF==6))
		{
		BLOCK_SIZE=BLOCKWIDTH_X-BLOCKWIDTH_X%DOF;
		}
	else{
		BLOCK_SIZE=BLOCKWIDTH_X;  
		}
		
	dim3 dimBlock(BLOCK_SIZE,BLOCKWIDTH_Y);
	dim3 dimGrid((int)ceil((float)(m*p*DOF)/(float)BLOCK_SIZE),((int)ceil((float)(n)/(float)BLOCKWIDTH_Y)));
				
	if (DOF==1)
		{
		MatMul_Kernel_tex_1_DOF<<<dimGrid,dimBlock>>>(d_coeff, d_x, d_y, d_idx, m, n, p, nos,DOF);
		}
	else{
		MatMul_Kernel_tex<<<dimGrid,dimBlock>>>(d_coeff, d_x, d_y, d_idx, m, n, p, nos, DOF);
		}
   
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
 

	if(_DBGFLAG) 
	{
		cudaThreadSynchronize();
		tend2=getclock();
		tkernel=tend2-tbegin2;
	}

	if (_DBGFLAG){
	cudaPrintfDisplay(stdout, true);
	cudaPrintfEnd();
	}
	
//Read y from the Device Memory

cudaMemcpy(y, d_y, size_xy, cudaMemcpyDeviceToHost); 
 
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

for(int i=0;i<m*n;i++)
printf("Y[%d]: %lf\n",i,y[i]);
//Free Device Memory
//cudaFree(d_coeff);
cudaFree(d_x);
cudaFree(d_y);
cudaFree(d_idx);
cudaFree(d_idy);
cudaFree(d_idz);

return 0;
}




