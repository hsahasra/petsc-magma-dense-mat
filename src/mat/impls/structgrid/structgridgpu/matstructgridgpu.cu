/*  -------------------------------------------------------------------- 
     This file extends structgrid data type to make use of GPUS. The new data type
     is structgridgpu. The implementation of the new datatype emulates the seqaijcusp
     implementation which is an extension to aij matrix format. 
     Author: Chekuri S. Choudary, RNET
             Daniel Lowell, ANL-MCS
*/


#define PETSCMAT_DLL
#include "../src/mat/impls/structgrid/structgridgpu/matstructgridgpu.h"

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
//#include <omp.h>
#include "../src/mat/impls/structgrid/matstructgrid.h"
#include "private/matimpl.h"
#include "private/vecimpl.h"
#include "matstructgridgpu.h"
#include "../include/private/petscimpl.h"


#define _DBGFLAG 0
#define KERNELVERSION 2

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


//static PetscScalar* d_coeff;



// ----------------------------------------------------------
// helper function for error checking
// pops the CUDA error stack and exits on nonzero error code
// written by: dlowell ANL-MCS
// ----------------------------------------------------------
void checkCUDAError(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) ); 
    exit(EXIT_FAILURE);//<-------------------use PETScError handle
  }
}



#undef __FUNCT__
#define __FUNCT__ "MatCheckCUDAError"
PetscErrorCode MatCheckCUDAError(const char *msg) {

  PetscFunctionBegin;
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) ); 
    fflush(NULL);
    PetscFunctionReturn(PETSC_ERR_LIB);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatCheckCUDAStatus"
PetscErrorCode MatCheckCUDAStatus(cudaError_t cs,const char *msg){

  PetscFunctionBegin;
  if(cs!=cudaSuccess){
    fprintf(stderr, "Cuda error!: %s: %s.\n",msg,cudaGetErrorString(cs));
    fflush(NULL);
    PetscFunctionReturn(PETSC_ERR_LIB);
  }
  PetscFunctionReturn(0);
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
PetscErrorCode  MatDestroy_SeqSGGPU(Mat B){
  PetscFunctionBegin;
  printf("MatDestroy_SeqSGGPU(Mat B)\n");
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  Mat_SeqSG* b=(Mat_SeqSG*)B->data;

  if (b->syncState != MAT_UNALLOC){
    //if (d_coeff) cudaFree(d_coeff);
    if(b->devptr){
      cudastatus = cudaFree(b->devptr); /* if(devX) cudaFree(devX); if(devY) cudaFree(devY); */
      ierr = MatCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr); 
      b->devptr=PETSC_NULL;
    }
  }
  b->syncState = MAT_UNALLOC;
  ierr = MatDestroy_SeqSG(B);CHKERRQ(ierr);
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
PetscErrorCode  MatCreate_SeqSGGPU(Mat B){
  PetscFunctionBegin;
  printf("MatCreate_SeqSGGPU(Mat B)\n");
  PetscErrorCode ierr;
  //cudaError_t cudastatus;
  ierr             = MatCreate_SeqSG(B);CHKERRQ(ierr);
  B->ops->mult     = MatMult_SeqSGGPU;
  B->ops->destroy  = MatDestroy_SeqSGGPU;
  ierr = PetscObjectChangeTypeName((PetscObject)B,MATSTRUCTGRIDGPU);CHKERRQ(ierr);
  /* Allocate device memory for matrix A */
  Mat_SeqSG *b = (Mat_SeqSG *) B->data;
  b->syncState = MAT_UNALLOC;
  PetscFunctionReturn(0);
}
EXTERN_C_END


/* --------------------------------------------------------------------
//     This function implements matrix vector multiplication for the
//     structgridgpu datatype. It calls a CUDA kernel to do matrix
//     multiplication on the GPU.
//     Author: Daniel Lowell, ANL-MCS, Chekuri S. Choudary, RNET
//--------------------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqSGGPU"
PetscErrorCode MatMult_SeqSGGPU(Mat mat, Vec x, Vec y){
  // int i;
  PetscErrorCode ierr;
  Mat_SeqSG *a=(Mat_SeqSG *) mat->data;
  /* PetscScalar *v, *xx,*yy; */
  PetscFunctionBegin;
  if(KERNELVERSION==1){
    /* Call to Jeswin's version */
    /* v = a->a;
       ierr = VecSet(y,0.0); CHKERRQ(ierr);
       ierr = VecGetArray(x, &xx); CHKERRQ(ierr);
       ierr = VecGetArray(y, &yy); CHKERRQ(ierr);
       ierr = SGCUDA_MatMult(v,xx,yy,a->idx,a->idy,
                            a->idz,a->m,a->n,a->p,a->stpoints,
                            &(mat->valid_GPU_matrix));
       CHKERRQ(ierr);
    ierr = SGCUDA_MatMult_v2(v,xx,yy,sparams,&(mat->valid_GPU_matrix));

       ierr = VecRestoreArray(x,&xx); CHKERRQ(ierr);
       ierr = VecRestoreArray(y,&yy); CHKERRQ(ierr); */
  }else if(KERNELVERSION==2){
    /* Call to dlowell's version */
    ierr = SGCUDA_MatMult_v2(mat,x,y); CHKERRQ(ierr);
  }
  ierr = PetscLogFlops(a->nz*a->stpoints); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
EXTERN_C_END


/* ------------------------------------------------------------------------
//   This function is the matrix vector multiplication kernel
//   structgridgpu datatype. This version uses shared memory for the write
//   back vector Y. Constant memory for reused constants and indices.
//   More offloading to registers might be possible as well.
//   written by: dlowell, ANL-MCS
//------------------------------------------------------------------------- */
EXTERN_C_BEGIN
__global__ void MatMul_Kernel_v2(double* A, double* X, double* Y){
  /* indices for local accesses */
  int tbtx, tbty, tbtz;
  int ix,iy,iz,j;

  int nos  = devparams.nos; /* set to register */
  int lda1 = devparams.lda1;
  int lda2 = devparams.lda2;
  int lda3 = devparams.lda3;
  int tilex = devparams.tile_x*devparams.tsizex;
  int tiley = devparams.tile_y*devparams.tsizey;
  int tilez = devparams.tile_z*devparams.tsizez;
  int Aindex, Xindex, index, offset;
 
  __shared__ double Ys[SHDSIZE][SHDSIZE][SHDSIZE];
  __shared__ double As[SHDSIZE][SHDSIZE][SHDSIZE];

  /* ------------------------------------------------------------------------ */

  for(iz=0;iz<tilez;iz+=devparams.tsizez){/*  tiles Z loop */
    tbtz = blockDim.z*blockIdx.z+threadIdx.z + iz;

    for(iy=0;iy<tiley;iy+=devparams.tsizey){/*  tiles Y loop */
      tbty = blockDim.y*blockIdx.y+threadIdx.y + iy;

      for(ix=0;ix<tilex;ix+=devparams.tsizex){/*  tiles X loop */
        tbtx = blockDim.x*blockIdx.x+threadIdx.x + ix;

        /* initialize current return-tile */
        Ys[threadIdx.z][threadIdx.y][threadIdx.x]=0.;

        /* adjusted index for global access */
        index = tbtz*lda2 + tbty*lda3 + tbtx;

        /* ......STENCIL........................................... */
        for(j=0;j<nos;j++){/* loop over stencil pattern */
	   offset= j*lda1;
           Aindex=offset+index;/* set up Aindex and read from global A into As tile */
           if(Aindex<devparams.matsize){
             if(A[Aindex]>0.)printf("A[%d]: %e\n",Aindex,A[Aindex]);
              As[threadIdx.z][threadIdx.y][threadIdx.x]=A[Aindex];//needs to be coalesced
           }else{
              As[threadIdx.z][threadIdx.y][threadIdx.x]=0.;
           }

           __syncthreads();

           /* set up Xindex for element-wise operation using stencil pattern */
           Xindex=(devparams.idz[j]*lda2 + devparams.idy[j]*lda3 + devparams.idx[j]) + index;
           if(Xindex<devparams.vecsize_x && Xindex>=0){
              Ys[threadIdx.z][threadIdx.y][threadIdx.x]+=As[threadIdx.z][threadIdx.y][threadIdx.x]*X[Xindex];
           }

        }/* end j-for */

        if(index<devparams.vecsize_y){
           Y[index]=Ys[threadIdx.z][threadIdx.y][threadIdx.x];/* global write back */
           if(Y[index]!=0.)printf("YOUT[%d]: %e\n",index,Y[index]);
        }

   }//end ix-for
   }//end iy-for
   }//end iz-for
}//end kernel_v2
EXTERN_C_END



/* ---------------------------------------------------------------------------------
//   This function is the wrapper function which sets up the device memory, transfers
//   data to and from the device, and calls the kernel. Error checking is done at 
//   each step. Timing stats are recorded using static vars.
//   written by: Daniel Lowell, ANL-MCS
//---------------------------------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "SGCUDA_MatMult_v2"
PetscErrorCode SGCUDA_MatMult_v2(Mat mat, Vec x,Vec y){

  PetscErrorCode ierr;
  PetscFunctionBegin;
  printf("Start SGCUDA_MatMult_v2\n");
  PetscInt i;
  Mat_SeqSG *a = (Mat_SeqSG *) mat->data;
  PetscScalar *A = a->a;
  Vec_SeqGPU* xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU* yd=(Vec_SeqGPU*)y->data;
  PetscScalar *X,*Y;

  static struct Stencilparams sparams;
  static double cumktime=0.;/* cummalitive kernel time */
  static double cumtime=0.;/* cummalitive call time */
  static unsigned int kcalls=0;/* number of kernel calls */
  double cs,ce,temp;
  float elapsedtime;        /* using CUDA device timer */
  cudaEvent_t start,stop;
  static unsigned char allocflag = 1;
  static double maxshared;
  static int bx,by,bz;/* number of blocks in 3-D */
  static int tx,ty,tz;/* number of threads ber block in 3-D */
  static int maxblocks_xy;
  static int maxblocks_z;
  static int xytile;

  static dim3 dimGrid;
  static dim3 dimBlock;
  cudaError_t cudastatus;



  if(_DBGFLAG){/* create CUDA events for timer */
     cudaEventCreate(&start);
     cudaEventCreate(&stop);
     cs=getclock();
  }

  /* Allocate and Memcpy Structured Matrix A
     The matrix remains the same throughout one iteration
     of the linear solver. The following uses a flag
     defined in the base class to check the status of the
     matrix. The matrix is copied to the GPU only if
     it has been changed on the CPU side
     This feature added by Chekuri S. Choudary */
  if(a->syncState == MAT_UNALLOC){  /* Allocate device memory for matrix A */
    printf("Allocating matrix A on GPU.\n");
    cudastatus=cudaSuccess;
    cudastatus=cudaMalloc((void**)&(a->devptr),a->matsize*sizeof(PetscScalar));
    ierr = MatCheckCUDAStatus(cudastatus,"a->devptr alloc in MatCreate_SeqSGGPU");CHKERRQ(ierr); 
    a->syncState = MAT_ALLOC;
  }

  if(A && a->syncState == MAT_ALLOC){
    a->syncState = MAT_CPU;
  }else if(!A){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"Matrix unallocated.");
  }

  if(a->syncState == MAT_CPU){
    /* copy over values of A to device memory */
    printf("Copying A to gpu. Size: %d\n",a->matsize*sizeof(PetscScalar));
    cudastatus=cudaSuccess;
    cudastatus=cudaMemcpy(a->devptr,A,a->matsize*sizeof(PetscScalar),cudaMemcpyHostToDevice);
    ierr = MatCheckCUDAStatus(cudastatus,"devA copy to device in SGCUDA_MatMult_v2");CHKERRQ(ierr);
    a->syncState = MAT_SYNCHED;
  }

  /* ..................................................... */

   if(xd->syncState==VEC_UNALLOC){
     SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"VectorX unallocated.");
   }else if(xd->syncState==VEC_CPU){
     /* copy over values of X to device memory */
     printf("Copying X to gpu.\n");
     ierr = VecGetArray(x,&X); CHKERRQ(ierr);
     ierr = VecRestoreArray(x,&X); CHKERRQ(ierr);
     xd->syncState = VEC_GPU;
   }


   /* ..................................................... */
   if(yd->syncState==VEC_UNALLOC){
     SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"VectorY unallocated.");
   }

   /* memset to 0. Vector Y on device */
   printf("Memset Y to gpu.\n");
   ierr = VecSet(y,0);CHKERRQ(ierr);
   printf("Done vector setup.\n");

   /* ..................................................... */

   if(allocflag){
     /* Set up blocks and thread numbers */
     maxshared = 49152.0/(float)(2.0*sizeof(double));
     if(a->p==1){
       xytile = pow(maxshared,0.5);/* square blocks */
       maxblocks_z = 1;
     }else{
       temp=maxshared/a->p;/* lop off z */
       xytile = pow(temp,0.5);/* xyblocks */
       maxblocks_z=ceil((float)SHDSIZE/(float)a->p);
     }
     maxblocks_xy = xytile/SHDSIZE;


     /* Set up blocks and thread numbers for columns */
     if(a->m <= SHDSIZE){
       tx = a->m; bx = 1; sparams.tile_x = 1; sparams.tsizex=1;
     }else{
       tx = SHDSIZE;
       bx = ceil((float)a->m/(float)SHDSIZE);/* create enough blocks */
       if(bx>maxblocks_xy){                 /* too many blocks created */
         bx = maxblocks_xy;                 /* set to max number of blocks allowed */
         sparams.tile_x=ceil((float)a->m/(float)(bx*SHDSIZE)); /* number of tiles */
         sparams.tsizex=bx*SHDSIZE;               /* tilesize is block-thread coverage */
       }else{
         sparams.tile_x=1; sparams.tsizex=1;
       }
     }

     /* Set up blocks and thread numbers for rows */
     if(a->n <= SHDSIZE){
       ty = a->n; by = 1; sparams.tile_y = 1; sparams.tsizey=1;
     }else{
       ty = SHDSIZE;
       by = ceil((float)a->n/(float)SHDSIZE);
       if(by > maxblocks_xy){
         by = maxblocks_xy;
         sparams.tile_y=ceil((float)a->n/(float)(by*SHDSIZE));
         sparams.tsizey=by*SHDSIZE;
       }else{
         sparams.tile_y=1; sparams.tsizey=1;
       }
     }

     /* Set up blocks and thread numbers for z */
     if(a->p <= SHDSIZE){
       tz = a->p;
       bz = 1;
       sparams.tile_z = 1;
       sparams.tsizez=1;
     }else{
       tz = SHDSIZE; bz = ceil((float)a->p/(float)SHDSIZE);
       if(bz > maxblocks_z){
         bz = maxblocks_z;
         sparams.tile_z=ceil((float)a->p/(float)(bz*SHDSIZE));
         sparams.tsizez=bz*SHDSIZE;
       }else{
         sparams.tile_z=1; sparams.tsizez=1;
       }
     }

     /* set grid shape */
     dimGrid.x = bx; dimGrid.y = by; dimGrid.z = bz;
     //dimGrid.x = 1; dimGrid.y = 1; dimGrid.z = 1;// bz;

     /* set block shape */
     dimBlock.x = tx; dimBlock.y = ty; dimBlock.z = tz;
     //dimBlock.x = 1; dimBlock.y = 1; dimBlock.z = 1;// tz;


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
     sparams.lda3=a->lda3;
     sparams.lda2=a->lda2;
     sparams.lda1=a->lda1;
     sparams.matsize=a->matsize;


     /* update constant memory with structured grid parameters */
     cudastatus=cudaMemcpyToSymbol("devparams",&sparams,sizeof(Stencilparams));
     ierr = MatCheckCUDAStatus(cudastatus,"symbol copy to device in SGCUDA_MatMult_v2");CHKERRQ(ierr);

     /* toggle off allocation flag */
     allocflag = 0;

   }/* end allocflag-if */


   /* grid and block shape & device config. debugging.................................. */
   unsigned int sharebytes;
   static unsigned char dbgflag = 1;
   if(_DBGFLAG && dbgflag){
     sharebytes = 2*tx*ty*tz*bx*by*bz*sizeof(double);
     printf("VecsizeX: %d, VecsizeY: %d, Matsize: %d\n",sparams.vecsize_x,sparams.vecsize_y,sparams.matsize);
     printf("(m, n, p, nos): (%d, %d, %d, %d)\n",sparams.m,sparams.n,sparams.p,sparams.nos);
     printf("MAXBLOCKS: %d maxshared: %lf\n",maxblocks_xy+maxblocks_z,maxshared);
     printf("blocks: (%d, %d, %d), threads per block: %d\n", bx,by,bz,tx*ty*tz );
     printf("Shared elements occupied: %0.3f SharedOccupied in Bytes: %d\n",sharebytes/49152.0,sharebytes);

     printf("Blocks x,y,z: (%d, %d, %d)\n",dimGrid.x,dimGrid.y,dimGrid.z);
     printf("Threads x,y,z: (%d, %d, %d)\n",dimBlock.x,dimBlock.y,dimBlock.z);

     printf("Blocks*ThreadsPer size x,y,z: (%d, %d, %d)\n",bx*tx,by*ty,bz*tz);
     printf("Tiles (x,y,z): (%d, %d, %d)\n",sparams.tile_x,sparams.tile_y,sparams.tile_z);
     printf("Tile Size (x,y,z): (%d, %d, %d)\n",sparams.tsizex,sparams.tsizey,sparams.tsizez);
     //dbgflag=0;
   }
   /* ...End config. debug section...................................................... */

   if(_DBGFLAG){//init kernel timer and debug settings
     ce=getclock();//end setup timer
     temp=ce-cs;
     // cudaPrintfInit();//start cuda printf environ.
     cudaEventRecord(start,0);//begin recording kernel
   }

   /*   printf("Launching kernel.\n");
   if(!a->devptr){
     printf("Hey! No a->devptr allocated!\n");
   } */


   /* Launch the kernel.......................................... */
   MatMul_Kernel_v2<<<dimGrid,dimBlock>>>(a->devptr,xd->devptr,yd->devptr);
   ierr = MatCheckCUDAError("CUDA Kernel launch status"); CHKERRQ(ierr);/* check for failure */
   /* ........................................................... */

   xd->syncState = VEC_GPU;
   yd->syncState = VEC_GPU;
   cudaDeviceSynchronize();

   if(_DBGFLAG){//end kernel timer and debug settings
     // cudaPrintfDisplay(stdout, true);//choose output
     // cudaPrintfEnd();//kill cuda printf environ
     cudaEventRecord(stop,0);
     cudaEventSynchronize(stop); // event barrier
     cudaEventElapsedTime(&elapsedtime,start,stop);
     cudaEventDestroy(start);
     cudaEventDestroy(stop);
     cs=getclock();
   }
   //ierr = VecRestoreArray(x,PETSC_NULL); CHKERRQ(ierr);
   //ierr = VecRestoreArray(y,PETSC_NULL); CHKERRQ(ierr);
   // Copy back Vector Y from Kernel
   /*if(y->valid_GPU_array==PETSC_CUSP_GPU){
     ierr = VecGetArray(y, &Y); CHKERRQ(ierr);
     cudastatus7=cudaMemcpy(Y,devY,vecsize_y,cudaMemcpyDeviceToHost);
     if(cudastatus7!=cudaSuccess){
       printf("Error on copy back Y, kernel status: %s\nExiting...\n\n",cudaGetErrorString(cudastatus7));
       if(devA) cudaFree(devA);
       if(devY) cudaFree(devY);
       if(devX) cudaFree(devX);
       PetscFunctionReturn(PETSC_ERR_MEM);
       }
     y->valid_GPU_array=PETSC_CUSP_BOTH;
     ierr = VecRestoreArray(y,&Y); CHKERRQ(ierr);
     ierr = PetscLogFlops(a->nz*a->stpoints); CHKERRQ(ierr);
   }*/

   if(_DBGFLAG){//final timer and debug settings
     //for(i=0;i<P.vecsize_y;i++)printf("Y[%d]: %lf\n",i,Y[i]);//for verification
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

   printf("returning from SGCUDA_MatMult_v2()\n");
   PetscFunctionReturn(0);
}
EXTERN_C_END








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
 #define BLOCKWIDTH 8
#define BLOCKWIDTH_X 8
#define BLOCKWIDTH_Y 8
#define BLOCKWIDTH_Z 8 


/*
  
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
   
  

 
 __global__ void MatMul_Kernel(PetscScalar * ptr_coeff, PetscScalar* ptr_x, PetscScalar* ptr_y, PetscInt *idx, PetscInt* idy, PetscInt* idz, PetscInt m, PetscInt n ,PetscInt p, PetscInt nos)
{

int tx= blockDim.x * blockIdx.x + threadIdx.x;
int ty= blockDim.y * blockIdx.y + threadIdx.y;
int l,i;
int xdisp,ydisp,zdisp,offset;
int lda1=m*n*p,lda2=m*n,lda3=m;
__shared__ PetscScalar y_sm[256];
__shared__ PetscScalar x_sm[324];

//initializing to the zero

// copying a Tile from Y into the shared Memory
y_sm[threadIdx.y*BLOCKWIDTH + threadIdx.x]=0;

//Copying a tile x into the shared Memory with 2 steps.

//Copying without the Ghost Cells  
x_sm[(threadIdx.y+1)*(BLOCKWIDTH+2) + (threadIdx.x+1)]=ptr_x[ty*lda3 + tx];

//Copying the Ghost Cells
// if (tx!=0)
// {
// if (threadIdx.x==0)
// x_sm[(threadIdx.y+1)*(BLOCKWIDTH+2) + threadIdx.x]=ptr_x[ty*lda3 + tx-1];
// }

// if (ty!=0)
// {
// if (threadIdx.y==0)
// x_sm[(threadIdx.y)* (BLOCKWIDTH+2) + threadIdx.x+1]=ptr_x[(ty-1)*lda3 + tx];
// }

// if (tx != n-1) // not sure about this
// {
// if (threadIdx.x==BLOCKWIDTH-1)
// x_sm[(threadIdx.y+1)*(BLOCKWIDTH+2) + threadIdx.x + 2]=ptr_x[ty*lda3 + tx + 1];
// }

// if (ty != m-1)
// {
// if (threadIdx.y==BLOCKWIDTH-1)
// x_sm[(threadIdx.y+2)*(BLOCKWIDTH+2) + threadIdx.x +1]=ptr_x[(ty+1)*lda3 + tx];
// }

//Copying the Ghost Cells
if (tx!=0 || ty!=0 || tx != n-1 || ty != m-1 || !(tx > n-1) || !(ty > m-1))
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
//MATMUL
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
	y_sm[threadIdx.y*BLOCKWIDTH + threadIdx.x]+= (ptr_coeff[offset + i*lda2 + ty*lda3 +tx] * x_sm[(i+zdisp)*lda2 + (threadIdx.y+ydisp +1)*
	(BLOCKWIDTH+2) + (threadIdx.x+xdisp+1)]); //forgetting Z currently.. I have to Fix it.
	}
	//removing i tempararily
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



*/

