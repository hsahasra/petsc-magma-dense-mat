#include <petscconf.h>
#include <petscsys.h>
//#include <petscerror.h>
PETSC_CUDA_EXTERN_C_BEGIN
#include <string.h>

#include <stdlib.h>
#include <float.h>
#include <petsc-private/vecimpl.h>          /*I "petscvec.h" I*/
#include <petscblaslapack.h>
#include <../src/vec/vec/impls/dvecimpl.h>
#include <../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h>
PETSC_CUDA_EXTERN_C_END




EXTERN_C_BEGIN

/* Misc constant memory variables (rarely used) */
__constant__ int     integerSymbol;
__constant__ int2    integer2Symbol;
__constant__ int3    integer3Symbol;
__constant__ int     devN;/* vector length */
__constant__ double  dblScalarValue;/* utility var */
__constant__ double2 dblScalar2Value;/* utility var */
__constant__ float   fltScalarValue;/* utility var */
__constant__ float2  fltScalar2Value;/* utility var */

/* error check variables */

static cudaError_t ccs[16];
static cudaError_t cms[16];

/* timer for vector functions */
#undef __FUNCT__
#define __FUNCT__ "vec_clock"
double vec_clock(){
  struct timezone tzp;
  struct timeval tp;
  gettimeofday (&tp, &tzp);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}


__device__ void orcu_warpReduce32(int tid, volatile double* reducts) {
  reducts[tid]+=reducts[tid+16];
  reducts[tid]+=reducts[tid+8];
  reducts[tid]+=reducts[tid+4];
  reducts[tid]+=reducts[tid+2];
  reducts[tid]+=reducts[tid+1];
}

__device__ void orcu_warpReduce64(int tid, volatile double* reducts) {
  reducts[tid]+=reducts[tid+32];
  reducts[tid]+=reducts[tid+16];
  reducts[tid]+=reducts[tid+8];
  reducts[tid]+=reducts[tid+4];
  reducts[tid]+=reducts[tid+2];
  reducts[tid]+=reducts[tid+1];
}

/* Function unrolls work of warp for reductions */
__device__ void warpDotReduce(volatile double* sdata, int tid){
  if(blockDim.x>=64)sdata[tid]+= sdata[tid+32];
  if(blockDim.x>=32)sdata[tid]+= sdata[tid+16];
  if(blockDim.x>=16)sdata[tid]+= sdata[tid+8];
  if(blockDim.x>=8) sdata[tid]+= sdata[tid+4];
  if(blockDim.x>=4) sdata[tid]+= sdata[tid+2];
  if(blockDim.x>=2) sdata[tid]+= sdata[tid+1];
}

/* Function unrolls work of warp for reductions */
__device__ void warpReduce(volatile double* sdata, int tid){
  sdata[tid]+= sdata[tid+32];
  sdata[tid]+= sdata[tid+16];
  sdata[tid]+= sdata[tid+8];
  sdata[tid]+= sdata[tid+4];
  sdata[tid]+= sdata[tid+2];
  sdata[tid]+= sdata[tid+1];
}

/* Function unrolls work of warp for reductions */
__device__ void warpMaxReduce(volatile double* sdata, int tid){
  sdata[tid]= fmax(sdata[tid],sdata[tid+32]);
  sdata[tid]= fmax(sdata[tid],sdata[tid+16]);
  sdata[tid]= fmax(sdata[tid],sdata[tid+8]);
  sdata[tid]= fmax(sdata[tid],sdata[tid+4]);
  sdata[tid]= fmax(sdata[tid],sdata[tid+2]);
  sdata[tid]= fmax(sdata[tid],sdata[tid+1]);
}


/* ---------------------------------------------------------
// helper function for error checking from kernel launches
// pops the CUDA error stack and exits on nonzero error code
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCheckCUDAError"
PetscErrorCode VecCheckCUDAError(const char *msg){

  PetscFunctionBegin;
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err){
    fprintf(stderr, "Cuda kernel error: %s: %s.\n", msg,cudaGetErrorString(err));
    fflush(NULL);
    PetscFunctionReturn(PETSC_ERR_LIB);
  }
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// helper function for error checking from status codes
// exits on nonzero error code, else does nothing
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCheckCUDAStatus"
PetscErrorCode VecCheckCUDAStatus(cudaError_t cs,const char *msg){
  PetscFunctionBegin;
    if(cs!=cudaSuccess){
      SETERRQ2(PETSC_COMM_SELF,0,"Cuda error: %s: %s.\n",msg,cudaGetErrorString(cs));
    }
  PetscFunctionReturn(0);
}
/* -------------------- end error checkers ------------------- */




/* ****************************************************************************
 This code is now included in CUDA SDK 4.1+ as cuRAND, so it may be obsolete

 *****************************************************************************
 * This is a shared memory implementation that keeps the full 625 words of state
 * in shared memory. Faster for heavy random work where you can afford 
 *  the shared memory. */
/* Init by single seed - single threaded as only used once */
__device__ void mt19937si(uint seed){
    int	i;
    if(threadIdx.x == 0){
	mtNexts = 0;
	s_seeds[0] = seed;
	for(i = 1;i < NNN;i++){
	   seed = (INIT_MULT * (seed^(seed >> 30))+i);
	   s_seeds[i] = seed;
	}
    }
    __syncthreads();/* Ensure mtNexts set */
    return;
}

/* Init by array - single threaded as only used once */
__device__ void mt19937sai(uint* seeds,uint length){
    int i,j,k;
    mt19937si(ARRAY_SEED);
    if(threadIdx.x==0){
     i=1; j=0;
     for(k = NNN>length?NNN:length;k!=0;k--){
        s_seeds[i] = (s_seeds[i]^((s_seeds[i-1]^(s_seeds[i-1] >> 30))*1664525)) + seeds[j] + j;
	if(++i >= NNN){
          s_seeds[0] = s_seeds[NNN-1];
	  i = 1;
        }
        if(++j>=length)j = 0;
     }
     for(k = NNN-1; k!=0;k--){
       s_seeds[i] = (s_seeds[i] ^ ((s_seeds[i-1]^(s_seeds[i-1]>>30))*1566083941))-i;
       if(++i >= NNN){
         s_seeds[0] = s_seeds[NNN-1];
	 i=1;
       }
     }
     s_seeds[0] = 0x80000000;/* MSB is 1; assuring non-zero initial array */ 
    }
    __syncthreads();				/* Needed for mt19937w() */
    return;
}

/* Return next MT random by increasing thread ID for 1-227 threads. */
__device__ uint mt19937s(void){
    int		kk;
    uint	x;
    uint	y;
    int		tid = threadIdx.x;

    kk = (mtNexts + tid) % NNN;
    __syncthreads();				/* Finished with mtNexts */

    if (tid == blockDim.x - 1)mtNexts = kk + 1;			/* Will get modded on next call */
    x = s_seeds[kk] & UPPER_MASK;
    if(kk < NNN - MMM){
      x |= (s_seeds[kk+1]&LOWER_MASK);
      y = s_seeds[kk+MMM];
    }else if(kk < NNN-1){
      x |= (s_seeds[kk+1]&LOWER_MASK);
      y = s_seeds[kk + (MMM-NNN)];
    }else{					/* kk == N - 1 */
      x |= (s_seeds[0]&LOWER_MASK);
      y = s_seeds[MMM - 1];
    }
    y ^= x >> 1;
    if (x & 1)y ^= MATRIX_A;
    __syncthreads();				/* All done before we update */

    s_seeds[kk] = y;
    y ^= (y >> 11);				/* Tempering */
    y ^= (y <<  7) & TEMPER1;
    y ^= (y << 15) & TEMPER2;
    y ^= (y >> 18);
    return y;
}

/* General shared memory version for any number of threads.
 * Note only up to 227 threads are run at any one time,
 * the rest loop and block till all are done. */
__device__ uint mt19937sl(void){
  int jj,kk,tid;
  uint x,y;
  tid = threadIdx.x;
  kk = (mtNexts + tid) % NNN;
  __syncthreads();				/* Finished with mtNexts */

  if(tid == blockDim.x - 1)mtNexts = kk + 1;	/* Will get modded on next call */
  jj = 0;
  do{
    if(jj <= tid && tid < jj + NNN - MMM){
      x = s_seeds[kk] & UPPER_MASK;
      if(kk < NNN - MMM){
         x |= (s_seeds[kk+1]&LOWER_MASK);
	 y = s_seeds[kk + MMM];
      }else if (kk < NNN-1){
         x |= (s_seeds[kk + 1]&LOWER_MASK);
	 y = s_seeds[kk + (MMM-NNN)];
      }else{				/* kk == N - 1 */
         x |= (s_seeds[0]&LOWER_MASK);
         y = s_seeds[MMM-1];
      }

      y ^= x >> 1;
      if(x & 1) y ^= MATRIX_A;
    }
    __syncthreads();			/* All done before we update */
    if(jj <= tid && tid < jj+NNN-MMM) s_seeds[kk] = y;
    __syncthreads();

  }while ((jj += NNN-MMM) < blockDim.x);
  y ^= (y >> 11);				/* Tempering */
  y ^= (y <<  7) & TEMPER1;
  y ^= (y << 15) & TEMPER2;
  y ^= (y >> 18);
  return y;
}

#undef __FUNCT__
#define __FUNCT__ "kernRandS"
__global__ void kernRandS(uint* seeds){
  mt19937sai(seeds,gridDim.x);
}

#undef __FUNCT__
#define __FUNCT__ "kernRand"
__global__ void kernRand(double *x, int* n){
  int tid = threadIdx.x + blockDim.x*blockIdx.x;
  uint rval;
  if(tid<*n){
    rval = mt19937sl();
    x[tid] = ((double)rval/(double)UINT_MAX);
    /* printf("RAND value[%d]: %0.13f, rval: %u UINT_MAX: %u\n",
       tid,x[tid],rval,UINT_MAX); */
  }
}

#undef __FUNCT__
#define __FUNCT__ "VecSetRandom_SeqGPU"
PetscErrorCode VecSetRandom_SeqGPU(Vec x,PetscRandom r){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i;
  uint *seeds=PETSC_NULL,*devseeds=PETSC_NULL;
  PetscScalar rval;
  dim3 dimBlock,dimGrid;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  if(xd->syncState==VEC_ALLOC || xd->syncState==VEC_CPU){
    for(i=0; i<x->map->n; i++){
       ierr = PetscRandomGetValue(r,&xd->cpuptr[i]);CHKERRQ(ierr);
    }
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }else if(xd->syncState==VEC_SYNCHED || xd->syncState==VEC_GPU){
    dimGrid.x=ceil((float)x->map->n/(float)TCOUNT);
    dimBlock.x=TCOUNT;
    ierr = PetscMalloc(dimGrid.x*sizeof(PetscInt),&seeds);CHKERRQ(ierr);
    for(i=0; i<dimGrid.x; i++){
       ierr = PetscRandomGetValue(r,&rval);CHKERRQ(ierr);
       seeds[i]=(uint)(UINT_MAX*rval);
    }
    cms[0] = cudaMalloc((void**)&devseeds,dimGrid.x*sizeof(uint));
    ccs[0]=cudaMemcpy(devseeds,seeds,dimGrid.x*sizeof(uint),cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(cms[0],"error in cudaMalloc");CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(ccs[0],"on copy H2D in VecSetRandom_SeqGPU");CHKERRQ(ierr);
    #endif

    kernRandS<<<dimGrid,dimBlock>>>(devseeds);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernRandS launch");CHKERRQ(ierr);
    #endif
    kernRand<<<dimGrid,dimBlock>>>(xd->devptr,xd->length);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernRand launch");CHKERRQ(ierr);
    #endif
    ierr = PetscFree(seeds);CHKERRQ(ierr);
    cudaDeviceSynchronize();
    cms[1] = cudaFree(devseeds);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(cms[1],"in cudaFree()");CHKERRQ(ierr);
    #endif
  }
  xd->syncState = VEC_GPU;
  PetscFunctionReturn(0);
}
/*------------------------ end random generator ------------------------*/



/*------------------------------ compare ------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "VecCompare_SeqGPU"
PetscErrorCode VecCompare_SeqGPU(Vec x, Vec y, PetscBool *same, PetscInt offset, PetscInt blocksize){
  PetscFunctionBegin;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  Vec_SeqGPU* yd = (Vec_SeqGPU*)y->data;
  if(xd->syncState!=yd->syncState||xd->syncState==VEC_ALLOC||yd->syncState==VEC_ALLOC){
    *same=PETSC_FALSE;
    PetscFunctionReturn(0);
  }
  PetscErrorCode ierr;
  dim3 dimGrid, dimBlock;
  if(blocksize && !offset){
    dimGrid.x=ceil((float)blocksize/(float)TCOUNT);
  } else {
    dimGrid.x=ceil((float)x->map->n/(float)TCOUNT);
  }
  dimBlock.x=TCOUNT;
  cudaError_t cudastatus;
  int *devsame=PETSC_NULL;
  int cpusame=0;
  int2 offset_bsize;
  offset_bsize.x = offset;
  offset_bsize.y = blocksize;
  if(xd->syncState==VEC_CPU && yd->syncState==VEC_CPU){
    ierr = PetscMemcmp((void*)&xd->cpuptr[offset],(void*)&yd->cpuptr[offset],blocksize,same);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  cudastatus = cudaMalloc((void**)&devsame,sizeof(int));
  ierr = VecCheckCUDAStatus(cudastatus,"error in device malloc");CHKERRQ(ierr);

  cudastatus=cudaMemcpyToSymbol("integer2Symbol",(void*)&offset_bsize,sizeof(int2),0,cudaMemcpyHostToDevice);
  ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);

  cudastatus=cudaMemcpyToSymbol("devN",(void*)&x->map->n,sizeof(int),0,cudaMemcpyHostToDevice);
  ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);

  kernCompare<<<dimGrid,dimBlock,2*dimBlock.x*sizeof(double)>>>(xd->devptr,yd->devptr,xd->length,yd->length,devsame);
  ierr = VecCheckCUDAError("kernCompare launch");CHKERRQ(ierr);

  cudastatus=cudaMemcpy(&cpusame,devsame,sizeof(int),cudaMemcpyDeviceToHost);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H in VecCompare_SeqGPU");CHKERRQ(ierr);

  if(cpusame==1)*same=PETSC_TRUE;
  else *same=PETSC_FALSE;
  cudastatus = cudaFree(devsame);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

extern __shared__ double sharedCompare[];
#undef __FUNCT__
#define __FUNCT__ "kernCompare"
__global__ void kernCompare(double* devX, double* devY, int* lx, int* ly, int* devsame){
  __shared__ unsigned char blockflag;
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int2 localOBS = integer2Symbol;
  int localn = localOBS.x+localOBS.y;
  int index = tid+localOBS.x;
  double value=0;
  double* chunkX = sharedCompare;
  double* chunkY = sharedCompare + blockDim.x;

  if(threadIdx.x==0)blockflag=0;
  __syncthreads();
  if(index<localn){
    /* read in values to shared */
    chunkX[threadIdx.x]=devX[index];
    chunkY[threadIdx.x]=devY[index];
    value = fabs(chunkX[threadIdx.x]-chunkY[threadIdx.x]);
    if(value>1e-16){
      #if(DEBUGVEC && VVERBOSE)
         printf("In kernCompare found an element mismatch: %e\n",value);
      #endif
      blockflag=1;
    }
    if(*lx!=*ly){
      #if(DEBUGVEC && VVERBOSE)
         printf("In kernCompare found length mismatch: lx: %d vs ly: %d\n",*lx,*ly);
      #endif
      blockflag=1;
    }
  }
  __syncthreads();
  if(threadIdx.x==0){
    if(blockflag)*devsame=0;
    else *devsame=1;
  }
  return;
}

/*-------------------------- end compare ----------------------------*/

/*----------------------- Vec info functions ------------------------*/

#undef __FUNCT__
#define __FUNCT__ "VecView_SeqGPU"
PetscErrorCode VecView_SeqGPU(Vec x,PetscViewer viewer){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  if(xd->syncState==VEC_GPU){
    ierr = VecCopyOverD2H(x,xd->cpuptr); CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  cudaDeviceSynchronize();
  int i;
  ierr = PetscObjectPrintClassNamePrefixType((PetscObject)x,viewer,"Vector Object");CHKERRQ(ierr);
  for(i=0;i<x->map->n;i++){
    PetscViewerASCIIPrintf(viewer,"%G\n",xd->cpuptr[i]);
  }
  /* ierr= PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);*/
  /* ierr =VecView_Seq_ASCII(x,viewer);CHKERRQ(ierr); */
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecGetLocalSize_SeqGPU"
PetscErrorCode VecGetLocalSize_SeqGPU(Vec x, PetscInt *localsize){
  PetscFunctionBegin;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecGetLocalSize_SeqGPU\n");
  #endif
  PetscValidHeaderSpecific(x,VEC_CLASSID,1);
  PetscValidIntPointer(localsize,2);
  PetscValidType(x,1);
  *localsize=x->map->n;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecGetSize_SeqGPU"
PetscErrorCode VecGetSize_SeqGPU(Vec x, PetscInt *globalsize){
  PetscFunctionBegin;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecGetSize_SeqGPU\n");
  #endif
  PetscValidHeaderSpecific(x,VEC_CLASSID,1);
  PetscValidIntPointer(globalsize,2);
  PetscValidType(x,1);
  *globalsize=x->map->N;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecCheck_SeqGPU"
PetscErrorCode VecCheck_SeqGPU(Vec x){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  dim3 dimgrid(ceil((float)x->map->n/((float)TCOUNT)),1,1);
  dim3 dimblocks(TCOUNT,1,1);
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  printf("******************************************\n");
  kernCheck<<<dimgrid,dimblocks>>>(xd->devptr,xd->length);
  ierr = VecCheckCUDAError("Call to kernCheck. "); CHKERRQ(ierr);
  cudaDeviceSynchronize();
  printf("******************************************\n");
  fflush(NULL);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "kernCheck"
__global__ void kernCheck(double* x, int* n){
  int tid = threadIdx.x + blockDim.x*blockIdx.x;
  if(tid<*n){
    #if(DEBUGVEC && VVERBOSE)
       printf("kernCheck: x[%d]: %e, length: %d\n",tid,x[tid],*n);
    #endif
  }
}
/*------------------------------ end info -------------------------------*/

/*---------------------------- copy functions ---------------------------*/


/* ---------------------------------------------------------
// Copies a block of memory from one array to another both 
// of which are on the device
// *** Currently nonfunctional ***
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCopyBlockDevice"
PetscErrorCode VecCopyBlockDevice(Vec d, Vec s, PetscInt doffset, PetscInt soffset, PetscInt blocksize){
  PetscFunctionBegin;
  printf("Call to VecCopyBlockDevice (**** EMPTY ****)\n");
  PetscFunctionReturn(0);
}


/* ---------------------------------------------------------
// Copies all elements from one allocated array to another.
// This is done asynchronously on the device Vec's streamID.
// Both array must be allocated and be of the same size otherwise
// PETSc will return an error.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCopyOverDevice"
PetscErrorCode VecCopyOverDevice(Vec d,Vec s){
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscInfo(d,"Copying vector on device only\n"); CHKERRQ(ierr);
  Vec_SeqGPU* dd = (Vec_SeqGPU*)d->data;
  Vec_SeqGPU* sd = (Vec_SeqGPU*)s->data;
  dim3 dimGrid;  dim3 dimBlock;
  if(s->map->n!=d->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  ccs[0]=cudaMemcpyAsync(dd->devptr,sd->devptr,
                    s->map->n*sizeof(double),cudaMemcpyDeviceToDevice,dd->streamid);
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// Helper function copies an integer array on device
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernCopyLen"
__global__ void kernCopyLen(int* ly,int* lx){
  if(threadIdx.x==0)*ly=*lx;
}


/* ---------------------------------------------------------
// Copies a block of elements from one host allocated array to an
// array allocated on the device. Does not check for allocation,
// only for failure if debugging is toggled.
// Copy is done asynchronously on the destination's Vec type's
// streamID.
// Both array must be allocated and be of the same size otherwise
// PETSc will return an error
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCopyBlockH2D"
PetscErrorCode VecCopyBlockH2D(Vec v,PetscScalar *y, PetscInt offset, PetscInt blocksize){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  ierr = PetscInfo(v,"Copying vec: cpu -> gpu\n"); CHKERRQ(ierr);
  ccs[0]=cudaMemcpyAsync(&(vd->devptr[offset]),y,
               blocksize*sizeof(double),cudaMemcpyHostToDevice,vd->streamid);
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// Copies all elements from one host allocated array to an
// array allocated on the device. Does not check for allocation,
// only for failure if debugging is toggled.
// Copy is done asynchronously on the destination's Vec type's
// streamID.
// Both array must be allocated and be of the same size otherwise
// PETSc will return an error
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCopyOverH2D"
PetscErrorCode VecCopyOverH2D(Vec v,PetscScalar *y){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  ccs[0]=cudaMemcpyAsync(vd->devptr,y,
               v->map->n*sizeof(double),cudaMemcpyHostToDevice,vd->streamid);
  ierr = PetscInfo(v,"Copying vec: cpu -> gpu\n"); CHKERRQ(ierr);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(ccs[0],"on copy H2D in VecCopyOverH2D");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


/* ---------------------------------------------------------
// Copies a block of elements from one device allocated array to an
// array allocated on the host. Does not check for allocation,
// only for failure if debugging is toggled.
// Copy is done asynchronously on the destination's Vec type's
// streamID.
// Both array must be allocated and be of the same size otherwise
// PETSc will return an error
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCopyBlockD2H"
PetscErrorCode VecCopyBlockD2H(Vec v,PetscScalar *y,PetscInt offset, PetscInt blocksize){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscInfo(v,"Copying vec: gpu -> cpu\n"); CHKERRQ(ierr);
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  ccs[0]=cudaMemcpyAsync(y,&(vd->devptr[offset]),
               blocksize*sizeof(double),cudaMemcpyDeviceToHost,vd->streamid);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(ccs[0],"on copy D2H in VecCopyBlockD2H");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


/* ---------------------------------------------------------
// Copies all elements from one device allocated array to an
// array allocated on the host. Does not check for allocation,
// only for failure if debugging is toggled.
// Copy is done asynchronously on the destination's Vec type's
// streamID.
// Both array must be allocated and be of the same size otherwise
// PETSc will return an error
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecCopyOverD2H"
PetscErrorCode VecCopyOverD2H(Vec v,PetscScalar *y){
  PetscErrorCode ierr;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  PetscFunctionBegin;
  ierr = PetscInfo(v,"Copying vec: gpu -> cpu\n"); CHKERRQ(ierr);
  ccs[0]=cudaMemcpyAsync(y,vd->devptr,v->map->n*sizeof(double),cudaMemcpyDeviceToHost,vd->streamid);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(ccs[0],"on copy D2H in VecCopyOverD2H");CHKERRQ(ierr); 
  #endif
  PetscFunctionReturn(0);
}

/*---------------------------- end copy functions --------------------------*/





/*------------------------------ set functions -----------------------------*/


/* ---------------------------------------------------------
// VecSetValues - Inserts or adds values into certain locations of a vector.
// INSERT and ADD VALUES both are implemented
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecSetValues_SeqGPU"
PetscErrorCode VecSetValues_SeqGPU(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora){
  PetscErrorCode ierr;
  PetscInt i;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  PetscInt *devix;
  PetscScalar *devy;
  dim3 grid,blocks;
  PetscFunctionBegin;
  ierr = PetscInfo(x,"setting gpu values\n"); CHKERRQ(ierr);
  if(xd->syncState==VEC_CPU || xd->syncState==VEC_SYNCHED){
    if(iora==INSERT_VALUES){
      for(i=0;i<ni;i++){
         xd->cpuptr[i]=y[i];
      }
      ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
      xd->syncState=VEC_SYNCHED;
    }else{/* ADD_VALUES */
       ierr = VecCheckCUDAStatus(cudaMalloc((void**)&devix,ni*sizeof(PetscInt)),
                     "cudaMalloc ADD_VALS");CHKERRQ(ierr);
       ierr = VecCheckCUDAStatus(cudaMalloc((void**)&devy,ni*sizeof(PetscScalar)),
                     "cudaMalloc ADD_VALS");CHKERRQ(ierr);
       ierr = VecCheckCUDAStatus(cudaMemcpy(devix,ix,ni*sizeof(PetscInt),cudaMemcpyHostToDevice),
                     "cudaMemcpy ADD_VALS");CHKERRQ(ierr);
       ierr = VecCheckCUDAStatus(cudaMemcpy(devy,y,ni*sizeof(PetscScalar),cudaMemcpyHostToDevice),
                     "cudaMemcpy ADD_VALS");CHKERRQ(ierr);
       grid.x = ceil((float)ni/((float)TCOUNT));
       blocks.x = TCOUNT;
       kernAddValues<<<grid,blocks>>>(xd->devptr,x->map->n,devix,ni,devy);
       cudaDeviceSynchronize();
       ierr=VecCheckCUDAStatus(cudaFree(devy),"cudaFree devy");CHKERRQ(ierr);
       ierr=VecCheckCUDAStatus(cudaFree(devix),"cudaFree devix");CHKERRQ(ierr);
       printf("Call to VecSetValues_SeqGPU 1: ADD_VALUES\n");
    }
  }else{
      if(iora==INSERT_VALUES){/* not efficient */
        PetscScalar yval=0;
        for(i=0;i<ni;i++){
          yval=y[i];
          ierr = VecCopyBlockH2D(x,&yval,ix[i],1);CHKERRQ(ierr);
        }
      }else{/* ADD_VALUES */

        ierr = VecCheckCUDAStatus(cudaMalloc((void**)&devix,ni*sizeof(PetscInt)),
                                  "cudaMalloc ADD_VALS");CHKERRQ(ierr);
        ierr = VecCheckCUDAStatus(cudaMalloc((void**)&devy,ni*sizeof(PetscScalar)),
                                  "cudaMalloc ADD_VALS");CHKERRQ(ierr);
        ierr = VecCheckCUDAStatus(cudaMemcpy(devix,ix,ni*sizeof(PetscInt),cudaMemcpyHostToDevice),
                                  "cudaMemcpy ADD_VALS");CHKERRQ(ierr);
        ierr = VecCheckCUDAStatus(cudaMemcpy(devy,y,ni*sizeof(PetscScalar),cudaMemcpyHostToDevice),
                                  "cudaMemcpy ADD_VALS");CHKERRQ(ierr);
        grid.x = ceil((float)ni/((float)TCOUNT));
        blocks.x = TCOUNT;
        kernAddValues<<<grid,blocks>>>(xd->devptr,x->map->n,devix,ni,devy);
        ierr=VecCheckCUDAError("call to kernAddValues.");CHKERRQ(ierr);
        cudaDeviceSynchronize();
        ierr=VecCheckCUDAStatus(cudaFree(devy),"cudaFree devy");CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cudaFree(devix),"cudaFree devix");CHKERRQ(ierr);
        printf("Call to VecSetValues_SeqGPU 2: ADD_VALUES\n");
      }
      xd->syncState=VEC_GPU;
  }
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// VecSet - Sets all values of an allocated device array to
// a scalar alpha. Checks if deivce array is allocated.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecSet_SeqGPU"
PetscErrorCode VecSet_SeqGPU(Vec xin,PetscScalar alpha){
  PetscFunctionBegin;
  #if(DEBUGVEC)
    PetscErrorCode ierr;
  #endif
  dim3 dimGrid, dimBlock;
  dimGrid.x = ceil((float)xin->map->n/((float)TCOUNT));
  dimBlock.x = TCOUNT;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)xin->data;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecSet_SeqGPU, alpha: %e\n",alpha);
  #endif
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,
            PETSC_ERR_MEM,"*** In VecSet_SeqGPU, Vec not allocated. ***\n");
  }else{
    kernSet<<<dimGrid,dimBlock>>>(xd->devptr,alpha,xin->map->n);
    #if(DEBUGVEC)
      #if(VVERBOSE)
         printf("In VecSet_SeqGPU: blocks: %d, threads: %d\n",dimGrid.x, dimBlock.x);
      #endif
      ierr = VecCheckCUDAError("Call to kernSet."); CHKERRQ(ierr);
    #endif
    xd->syncState=VEC_GPU;
  }
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// Device kernel called by VecSetValues when insert type is 
// ADD VALUES.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernAddValues"
__global__ void kernAddValues(double* x, int n, int* xi, int ni,double *y){
  unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
  if(tid<ni) x[xi[tid]] += y[tid];
}

/* ---------------------------------------------------------
// Called by VecSet and set all values of an allocated array
// to a scalar alpha.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernSet"
__global__ void kernSet(double* x, double alpha, int n){
  unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
  if(tid<n) x[tid] = alpha;
}


/* ---------------------------------------------------------
// VecScale - Scales all values of an allocated device array
// by a scalar alpha. Checks if deivce array is allocated.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecScale_SeqGPU"
PetscErrorCode VecScale_SeqGPU(Vec x, PetscScalar alpha){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  dim3 dimGrid,dimBlock;
  dimGrid.x=ceil((float)x->map->n/((float)TCOUNT));
  dimBlock.x=TCOUNT;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  #if(DEBUGVEC && VVERBOSE)
     printf("VecScale_SeqGPU...alpha: %e\n",alpha);
  #endif
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,
            PETSC_ERR_MEM,
            "*** In call to VecScale_SeqGPU, arg Vec xin has not been allocated. ***\n");
  }else if(xd->syncState==VEC_CPU){
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  cudaDeviceSynchronize();
  if(alpha==0.){
    ccs[0] = cudaMemset(xd->devptr,0,x->map->n*sizeof(double));
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[0],"error in cudaMemset");CHKERRQ(ierr);
    #endif
  }else if (alpha != 1.0){
    kernScale<<<dimGrid,dimBlock>>>(xd->devptr,alpha,x->map->n);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("Call to kernScale."); CHKERRQ(ierr);
    #endif
  }
  xd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// Device array scaling kernel called by VecScale
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernScale"
__global__ void kernScale(double* x, double alpha, int n){
  unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
  if(tid<n) x[tid] *= alpha;
}
/*---------------------------- end set and scale ---------------------------*/


/*-------------------------- dot product functions -------------------------*/

/* ---------------------------------------------------------
// Does the Dot product of a transposed vector yin...
// *** Currently non-functional ***
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecTDot_SeqGPU"
PetscErrorCode VecTDot_SeqGPU(Vec xin,Vec yin,PetscScalar *z){
  PetscFunctionBegin;
  printf("VecTDot_SeqGPU (***EMPTY***)\n");
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// Computes the dot product of two vectors
// Checks for size mismatch and will synch the vectors to the
// device if needed.
// Orio tuned kernels are implemented for 3 size ranges:
// Manual tuned kernels are also included.
// Timers are available if toggled.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecDot_SeqGPU"
PetscErrorCode VecDot_SeqGPU(Vec x,Vec y,PetscScalar *z){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(x->map->n!=y->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
#if(VTIMER)
  double start,finish,elapsed;
  static double mint,maxt=0.,cumt=0.,avg=0.;
  static int ccnt=0;
  //start = vec_clock();
#endif
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  if(xd->syncState==VEC_CPU){
    #if(DEBUGVEC && VVERBOSE)
       printf("xd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  if(yd->syncState==VEC_CPU){
    #if(DEBUGVEC && VVERBOSE)
       printf("yd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  double *devScratch;
  dim3 dimGrid, dimBlock;
#if(VMANDOT)
cudaStream_t* dotstream;
  PetscInt i,chunks=0,segment,scratchsize;
  float threadscale =(DOTMPLIER*CHUNKWIDTH);
  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) /threadscale);
  dotstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));
  /* make sure the segment size for each chunk is correct */
  if(chunks>1) segment = (int) threadscale;
  else segment = x->map->n;
  dimGrid.x=ceil(((float)segment)/(float)THRDOTCNT);
  dimBlock.x = THRDOTCNT;
  /* allocate gridwide scratch array */
  scratchsize=chunks*dimGrid.x;
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize*sizeof(double));/* scratch pad */
  #if(DEBUGVEC)
    #if(VVERBOSE)
       printf("Call to VecDot_SeqGPU, chunks: %d segsize: %d grid: %d\n",chunks,segment,dimGrid.x);
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecDot_SeqGPU");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();/* make sure everyone is ready  */
  #if(VTIMER)
    start = vec_clock();
  #endif
  for(i=0;i<chunks;i++){  /* streaming async kernel calls */
    cudaStreamCreate(&(dotstream[i]));
    /* Overlapping execution */
    kernDot<<<dimGrid,dimBlock,0,dotstream[i]>>>(xd->devptr,yd->devptr,
                                                          segment,
                                                          x->map->n,
                                                          i,
                                                          devScratch+i*dimGrid.x);
  }
  dimBlock.x  = THRDOTCNT2;
  cudaDeviceSynchronize();

  while(scratchsize>1){/* begin next reduction */
    dimGrid.x = ceil((float)scratchsize/(float)THRDOTCNT2);
    if(dimGrid.x>MAXBLOCKS){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_LIB,"Grid size too large for gpu capabilities.");
    }
    kernRedDot<<<dimGrid,dimBlock>>>(scratchsize,devScratch);
    scratchsize = dimGrid.x;
  }
  for(i=0;i<chunks;i++){
     cudaStreamDestroy(dotstream[i]);
  }
  free(dotstream);
 
#else /* use Orio kernels */


  if(y->map->n>=1e7){
    int nthreads=320;
    int nstreams=4;
    /* Set per function L1 cache config */
    cudaFuncSetCacheConfig(orcu_dotkernel_1e7,cudaFuncCachePreferL1);
    dimBlock.x=nthreads;
    dimGrid.x=112;
    /*create streams*/
    int istream, soffset, boffset;
    cudaStream_t stream[nstreams+1];
    for(istream=0;istream<=nstreams;istream++) cudaStreamCreate(&stream[istream]);
    /* divide up stream work */
    int chunklen=x->map->n/nstreams;
    int chunkrem=x->map->n%nstreams;
    int blks4chunk=dimGrid.x/nstreams;
    if(dimGrid.x%nstreams!=0) blks4chunk++;
    int blks4chunks=blks4chunk*nstreams;
    /* allocate scratch pad memory */
    cudaMalloc((void**)&devScratch,(dimGrid.x+1)*sizeof(double));

    cudaDeviceSynchronize();
    for(istream=0; istream<nstreams;istream++) {
      soffset=istream*chunklen;
      boffset=istream*blks4chunk;
      orcu_dotkernel_1e7<<<blks4chunk,dimBlock,0,stream[istream]>>>
                        (chunklen,xd->devptr+soffset,yd->devptr+soffset,devScratch+boffset);
    }
    if(chunkrem!=0){/* do remaining work */
      soffset=istream*chunklen;
      boffset=istream*blks4chunk;
      orcu_dotkernel_1e7<<<blks4chunk,dimBlock,0,stream[istream]>>>
                        (chunkrem,xd->devptr+soffset,yd->devptr+soffset,devScratch+boffset);
      blks4chunks++ ;
    }
    int orcu_blks=blks4chunks;
    int orcu_n;
    while (orcu_blks>1){/* second stage reduction */
      orcu_n=orcu_blks;
      orcu_blks=(orcu_blks+319)/320;
      orcu_dotblksum_1e7<<<orcu_blks,320>>>(orcu_n,devScratch);
    }
    /* destroy streams */
    for (istream=0; istream<=nstreams; istream++)cudaStreamDestroy(stream[istream]);


  }else if(y->map->n>=1e6){
    int nthreads=160;
    int nstreams=2;
    cudaFuncSetCacheConfig(orcu_dotkernel_1e6,cudaFuncCachePreferL1);
    dimBlock.x=nthreads;
    dimGrid.x=112;
    /*create streams*/
    int istream, soffset, boffset;
    cudaStream_t stream[nstreams+1];
    for(istream=0;istream<=nstreams;istream++) cudaStreamCreate(&stream[istream]);
    /* divide up stream work */
    int chunklen=x->map->n/nstreams;
    int chunkrem=x->map->n%nstreams;
    int blks4chunk=dimGrid.x/nstreams;
    if(dimGrid.x%nstreams!=0) blks4chunk++;
    int blks4chunks=blks4chunk*nstreams;
    /* allocate scratch pad memory */
    cudaMalloc((void**)&devScratch,(dimGrid.x+1)*sizeof(double));

    cudaDeviceSynchronize();
    for(istream=0; istream<nstreams;istream++) {
      soffset=istream*chunklen;
      boffset=istream*blks4chunk;
      orcu_dotkernel_1e6<<<blks4chunk,dimBlock,0,stream[istream]>>>
                        (chunklen,xd->devptr+soffset,yd->devptr+soffset,devScratch+boffset);
    }
    if(chunkrem!=0){/* do remaining work */
      soffset=istream*chunklen;
      boffset=istream*blks4chunk;
      orcu_dotkernel_1e6<<<blks4chunk,dimBlock,0,stream[istream]>>>
                        (chunkrem,xd->devptr+soffset,yd->devptr+soffset,devScratch+boffset);
      blks4chunks++ ;
    }
    int orcu_blks=blks4chunks;
    int orcu_n;
    while (orcu_blks>1){/* second stage reduction */
      orcu_n=orcu_blks;
      orcu_blks=(orcu_blks+159)/160;
      orcu_dotblksum_1e6<<<orcu_blks,160>>>(orcu_n,devScratch);
    }
    /* kill streams */
    for (istream=0; istream<=nstreams; istream++)cudaStreamDestroy(stream[istream]);
 }else{
      /*calculate device dimensions*/
      dimBlock.x=256;
      dimGrid.x=112;
      /* allocate scratch pad memory */
      cudaMalloc((void**)&devScratch,(dimGrid.x+1)*sizeof(double));
      orcu_dotkernel_1e5<<<dimGrid,dimBlock>>>
                 (y->map->n,yd->devptr,xd->devptr,devScratch);

      int orcu_blks=dimGrid.x;
      int orcu_n;
      while (orcu_blks>1){/* second stage reduction */
        orcu_n=orcu_blks;
        orcu_blks=(orcu_blks+255)/256;
        orcu_dotblksum_1e5<<<orcu_blks,256>>>(orcu_n,devScratch);
      }
  }/*end Orio vecsize if */
#endif /* end VMANUAL */
#if(VTIMER)
  finish=vec_clock();
#endif

  ccs[4]=cudaMemcpy(z,devScratch,sizeof(double),cudaMemcpyDeviceToHost);/* copy back */
  /* delete scratch memory */
  cms[3] = cudaFree(devScratch);
#if(DEBUGVEC)
  ierr = VecCheckCUDAStatus(ccs[4],"on cudaMemcpy(devScratch)");CHKERRQ(ierr);
#endif

 #if(DEBUGVEC)
   #if(VVERBOSE)
       printf("Zdot: %e\n",*z);
   #endif
    ierr = VecCheckCUDAStatus(cms[3],"on cudaFree(devScratch)");CHKERRQ(ierr);
 #endif

 #if(VTIMER)
    // finish=vec_clock();
    elapsed=finish-start;
    cumt+=elapsed;
    if(!ccnt++){
      maxt=mint=avg=elapsed;
    }else{
      maxt=elapsed>maxt?elapsed:maxt;
      mint=elapsed<mint?elapsed:mint;
      avg=cumt/ccnt;
      if(!(ccnt%ITSHOW)){
        printf("VecDot_SeqGPU calls: %d, max: %e, min: %e, average: %e\n",
               ccnt,maxt,mint,avg);
      }
    }
  #endif
  PetscFunctionReturn(0);
}

#if(VMANDOT)
/* ---------------------------------------------------------
// Manual tuned second stage dot product parallel reduction.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernRedDot"
__global__ void kernRedDot(int n,double* scratch){/* reduction kernel */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunk[THRDOTCNT2];
  double mysum =(tid<n)?scratch[tid]:0.;
  if(threadIdx.x>127)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<128)mysum+=chunk[threadIdx.x+128];
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<64)mysum+=chunk[threadIdx.x+64];
  else return;
  chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<32) warpReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}

/* ---------------------------------------------------------
// Manual tuned first stage dot product parallel reduction.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernDot"
__global__ void kernDot(double* devX, double* devY,int segmentsize,
                        int arrsize, int offset, double* scratch){
  __shared__ double chunk[THRDOTCNT];
  unsigned int item = segmentsize*offset+blockDim.x*blockIdx.x+threadIdx.x;
  double mysum=(item<arrsize)?(devX[item]*devY[item]):0.;
  if(threadIdx.x>127)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<128)mysum+=chunk[threadIdx.x+128];
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<64)mysum+=chunk[threadIdx.x+64];
  else return;
  chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<32) warpReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}

#else

/* 1e5 tuned kernels */
#undef __FUNCT__
#define __FUNCT__ "orcu_dotkernel_1e5"
__global__ void orcu_dotkernel_1e5(const int n, double* y, double* x, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  double orcu_var8193=0;
  for (int i=tid; i<=n-1; i+=gsize) {
    orcu_var8193=orcu_var8193+x[i]*y[i];
  }
  /*reduce single-thread results within a block*/
  __shared__ double orcu_vec8194[256];
  orcu_vec8194[threadIdx.x]=orcu_var8193;
  __syncthreads();
  if (threadIdx.x<128) 
    orcu_vec8194[threadIdx.x]+=orcu_vec8194[threadIdx.x+128];
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec8194[threadIdx.x]+=orcu_vec8194[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec8194);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec8194[0];
}

#undef __FUNCT__
#define __FUNCT__ "orcu_dotblksum_1e5"
__global__ void orcu_dotblksum_1e5(int orcu_n, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double orcu_vec8194[256];
  if (tid<orcu_n) 
    orcu_vec8194[threadIdx.x]=reducts[tid];
  else 
    orcu_vec8194[threadIdx.x]=0;
  __syncthreads();
  if (threadIdx.x<128) 
    orcu_vec8194[threadIdx.x]+=orcu_vec8194[threadIdx.x+128];
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec8194[threadIdx.x]+=orcu_vec8194[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec8194);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec8194[0];
}

/* 1e6 tuned kernels */
#undef __FUNCT__
#define __FUNCT__ "orcu_dotkernel_1e6"
__global__ void orcu_dotkernel_1e6(const int n, double* y, double* x, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_y[160];
  __shared__ double shared_x[160];
  double orcu_var16389=0;
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_y[threadIdx.x]=y[i];
    shared_x[threadIdx.x]=x[i];
    orcu_var16389=orcu_var16389+shared_x[threadIdx.x]*shared_y[threadIdx.x];
  }
  /*reduce single-thread results within a block*/
  __shared__ double orcu_vec16390[160];
  orcu_vec16390[threadIdx.x]=orcu_var16389;
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec16390[threadIdx.x]+=orcu_vec16390[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec16390);
  if (threadIdx.x>=128&&threadIdx.x<144) 
    orcu_warpReduce32(threadIdx.x,orcu_vec16390);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec16390[0]+orcu_vec16390[128];
}


#undef __FUNCT__
#define __FUNCT__ "orcu_dotblksum_1e6"
__global__ void orcu_dotblksum_1e6(int orcu_n, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double orcu_vec16390[160];
  if (tid<orcu_n) 
    orcu_vec16390[threadIdx.x]=reducts[tid];
  else 
    orcu_vec16390[threadIdx.x]=0;
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec16390[threadIdx.x]+=orcu_vec16390[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec16390);
  if (threadIdx.x>=128&&threadIdx.x<144) 
    orcu_warpReduce32(threadIdx.x,orcu_vec16390);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec16390[0]+orcu_vec16390[128];
}



/* 1e7 tuned kernels */
#undef __FUNCT__
#define __FUNCT__ "orcu_dotkernel_1e7"
__global__ void orcu_dotkernel_1e7(const int n, double* y, double* x, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_y[160];
  __shared__ double shared_x[160];
  double orcu_var16389=0;
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_y[threadIdx.x]=y[i];
    shared_x[threadIdx.x]=x[i];
    orcu_var16389=orcu_var16389+shared_x[threadIdx.x]*shared_y[threadIdx.x];
  }
  /*reduce single-thread results within a block*/
  __shared__ double orcu_vec16390[160];
  orcu_vec16390[threadIdx.x]=orcu_var16389;
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec16390[threadIdx.x]+=orcu_vec16390[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec16390);
  if (threadIdx.x>=128&&threadIdx.x<144) 
    orcu_warpReduce32(threadIdx.x,orcu_vec16390);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec16390[0]+orcu_vec16390[128];
}


#undef __FUNCT__
#define __FUNCT__ "orcu_dotblksum_1e7"
__global__ void orcu_dotblksum_1e7(int orcu_n, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double orcu_vec16390[160];
  if (tid<orcu_n) 
    orcu_vec16390[threadIdx.x]=reducts[tid];
  else 
    orcu_vec16390[threadIdx.x]=0;
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec16390[threadIdx.x]+=orcu_vec16390[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec16390);
  if (threadIdx.x>=128&&threadIdx.x<144) 
    orcu_warpReduce32(threadIdx.x,orcu_vec16390);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec16390[0]+orcu_vec16390[128];
}

#endif /* end VMANDOT -#if */


/* ---------------------------------------------------------
// Computes multiple dot products over a single x and multiple y's
// val[i] = x . y[i]
// Just a for-loop wrapper call to VecDot_SeqGPU
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecMDot_SeqGPU"
PetscErrorCode  VecMDot_SeqGPU(Vec x,PetscInt nv,const Vec y[],PetscScalar val[]){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt       i;
  #if(DEBUGVEC && VVERBOSE)
    printf("Call to VecMDot_SeqGPU\n");
  #endif
  for (i=0; i<nv; i++) {
    ierr = VecDot_SeqGPU(x,y[i],&val[i]);CHKERRQ(ierr);
    if(PetscIsInfOrNanScalar(val[i])){
      SETERRQ1(((PetscObject)x)->comm,PETSC_ERR_FP,"Infinite or not-a-number generated in mdot, entry %D",i);
    }
  }
  PetscFunctionReturn(0);
}
/*----------------------------- end dot ----------------------------- */



/* ---------------------------------------------------------
// AXPBY currently not implemented
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecAXPBY_SeqGPU"
PetscErrorCode VecAXPBY_SeqGPU(Vec yin,PetscScalar beta,PetscScalar alpha,Vec xin){
  /* Y = b*Y + a*X */
  PetscFunctionBegin;
  printf("Call to VecAXPBY_SeqGPU (***EMPTY***)\n");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecAYPX_SeqGPU"
PetscErrorCode VecAYPX_SeqGPU(Vec Y, PetscScalar alpha, Vec X)
{
  /* Y = X + alpha Y */
  /* reference implementation -- not optimized */
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecScale(Y,alpha);CHKERRQ(ierr);
  ierr = VecAXPY(Y,1.0,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* ---------------------------------------------------------
// WAXPY
// Checks for vector size mismatch for each vector.
// Copies down the vector(s) to the GPU if needed
// Manual implmentation and Orio tuned kernels can be toggled
// Timers are available if toggled
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecWAXPY_SeqGPU"
PetscErrorCode VecWAXPY_SeqGPU(Vec w,PetscScalar alpha,Vec x,Vec y){
  /* w = y + alpha*x */
  PetscFunctionBegin;
  #if(VTIMER)
    double start,finish,elapsed;
    static double mint,maxt=0.,cumt=0.,avg=0.;
    static int ccnt=0;
    start = vec_clock();
  #endif
  PetscErrorCode ierr;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)w->data;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid, dimBlock;
  #if(DEBUGVEC && VVERBOSE)
     printf("VecWAXPY_SeqGPU...alpha: %e\n",alpha);
  #endif
  if(x->map->n!=y->map->n || w->map->n!=y->map->n || w->map->n!=x->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  if(yd->syncState==VEC_CPU){/* synch up y */
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  if(xd->syncState==VEC_CPU){/* synch up x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }

#if(VMANWXP)
  dimGrid.x=ceil((float)y->map->n/(float)AXPYTCOUNT);
  dimBlock.x=AXPYTCOUNT;
  cudaDeviceSynchronize();
  if(alpha==0.0){
    ierr = VecCopyOverDevice(w,y);CHKERRQ(ierr);
    cudaDeviceSynchronize();
  }else if(alpha==1.0){
    kernWXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,x->map->n,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWXPY");CHKERRQ(ierr); 
    #endif
  }else if(alpha==-1.0){
    kernWYMX<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,x->map->n,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWYMX");CHKERRQ(ierr);
    #endif
  }else{
    kernWAXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,alpha,x->map->n,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWAXPY");CHKERRQ(ierr); 
    #endif
  }
#else
  if(w->map->n>=1e7){
    int nthreads=352;
    int nstreams=5;
    cudaFuncSetCacheConfig(orcu_waxpykernel_1e7,cudaFuncCachePreferL1);
    dimBlock.x=nthreads;
    dimGrid.x=70;
    /*create streams*/
    int istream, soffset;
    cudaStream_t stream[nstreams+1];
    for (istream=0;istream<=nstreams;istream++) cudaStreamCreate(&stream[istream]);
    int chunklen=x->map->n/nstreams;
    int chunkrem=x->map->n%nstreams;


    /*invoke device kernel*/
    int blks4chunk=dimGrid.x/nstreams;
    if (dimGrid.x%nstreams!=0) blks4chunk++ ;
    for (istream=0; istream<nstreams; istream++ ) {
      soffset=istream*chunklen;
      orcu_waxpykernel_1e7<<<blks4chunk,dimBlock,0,stream[istream]>>>
                      (chunklen,alpha,xd->devptr+soffset,yd->devptr+soffset,wd->devptr+soffset);
    }
    if (chunkrem!=0) {
      soffset=istream*chunklen;
      orcu_waxpykernel_1e7<<<blks4chunk,dimBlock,0,stream[istream]>>>
                      (chunkrem,alpha,xd->devptr+soffset,yd->devptr+soffset,wd->devptr+soffset);
    }
    cudaDeviceSynchronize();
    for (istream=0; istream<=nstreams; istream++) cudaStreamDestroy(stream[istream]);

  }else if(w->map->n>=1e6){
    int nthreads=512;
    int nstreams=5;
    /*calculate device dimensions*/
    dimBlock.x=nthreads;
    dimGrid.x=70;
    /*create streams*/
    int istream, soffset;
    cudaStream_t stream[nstreams+1];
    for (istream=0; istream<=nstreams; istream++) cudaStreamCreate(&stream[istream]);
    int chunklen=x->map->n/nstreams;
    int chunkrem=x->map->n%nstreams;
    cudaFuncSetCacheConfig(orcu_waxpykernel_1e6,cudaFuncCachePreferL1);

    /*invoke device kernel*/
    int blks4chunk=dimGrid.x/nstreams;
    if (dimGrid.x%nstreams!=0) blks4chunk++ ;
    for (istream=0; istream<nstreams; istream++ ) {
      soffset=istream*chunklen;
      orcu_waxpykernel_1e6<<<blks4chunk,dimBlock,0,stream[istream]>>>
                          (chunklen,alpha,xd->devptr+soffset,yd->devptr+soffset,wd->devptr+soffset);
    }
    if (chunkrem!=0) {
      soffset=istream*chunklen;
      orcu_waxpykernel_1e6<<<blks4chunk,dimBlock,0,stream[istream]>>>
                          (chunklen,alpha,xd->devptr+soffset,yd->devptr+soffset,wd->devptr+soffset);
    }
    cudaDeviceSynchronize();
    for (istream=0; istream<=nstreams; istream++ ) cudaStreamDestroy(stream[istream]);

  }else{

    cudaFuncSetCacheConfig(orcu_waxpykernel_1e5,cudaFuncCachePreferL1);
    int nthreads=320;
    dimBlock.x=nthreads;
    dimGrid.x=112;
    /*invoke device kernel*/
    orcu_waxpykernel_1e5<<<dimGrid,dimBlock>>>(x->map->n,alpha,xd->devptr,yd->devptr,wd->devptr);

    cudaError_t err=cudaGetLastError();
    if (cudaSuccess!=err) {
      printf("CUDA runtime error: %s@",cudaGetErrorString(err));
    }
    cudaDeviceSynchronize();
 

  }
#endif /* manual if */
  wd->syncState=VEC_GPU;
#if(VTIMER)
    cudaDeviceSynchronize();
    finish=vec_clock();
    elapsed=finish-start;
    cumt+=elapsed;
    if(!ccnt++){
      maxt=mint=avg=elapsed;
    }
    maxt=elapsed>maxt?elapsed:maxt;
    mint=elapsed<mint?elapsed:mint;
    avg=cumt/ccnt;
    printf("VecWAXPY_SeqGPU calls: %d, max: %e, min: %e, average: %e\n",
             ccnt,maxt,mint,avg);
#endif
  PetscFunctionReturn(0);
}



#if(VMANWXP)
#undef __FUNCT__
#define __FUNCT__ "kernWAXPY"
__global__ void  kernWAXPY(double* devY,double* devX, double alpha, int vlen, double* devW){
  /* w <- y + alpha*x */
  unsigned int tid = blockIdx.x*blockDim.x+threadIdx.x;
  if(tid<vlen){
    devW[tid]=devY[tid]+alpha*devX[tid];
  }
}

#undef __FUNCT__
#define __FUNCT__ "kernWXPY"
__global__ void  kernWXPY(double* devY,double* devX, int vlen, double* devW){
 /* w <- y + x */
  unsigned int tid = blockIdx.x*blockDim.x+threadIdx.x;
  if(tid<vlen){
    devW[tid]=devY[tid]+devX[tid];
  }
}

#undef __FUNCT__
#define __FUNCT__ "kernWYMX"
__global__ void  kernWYMX(double* devY,double* devX, int vlen, double* devW){
 /* w <- y - x */
  unsigned int tid = blockIdx.x*blockDim.x+threadIdx.x;
  if(tid<vlen){
    devW[tid]=devY[tid]-devX[tid];
  }
}

#else

__global__ void orcu_waxpykernel_1e7(const int n, double a, double* x, double* y, double* w) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_y[352];
  __shared__ double shared_x[352];
  __shared__ double shared_w[352];
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_y[threadIdx.x]=y[i];
    shared_x[threadIdx.x]=x[i];
    shared_w[threadIdx.x]=a*shared_x[threadIdx.x]+shared_y[threadIdx.x];
    w[i]=shared_w[threadIdx.x];
  }
}

__global__ void orcu_waxpykernel_1e6(const int n, double a, double* x, double* y, double* w) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_y[512];
  __shared__ double shared_x[512];
  __shared__ double shared_w[512];
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_y[threadIdx.x]=y[i];
    shared_x[threadIdx.x]=x[i];
    shared_w[threadIdx.x]=a*shared_x[threadIdx.x]+shared_y[threadIdx.x];
    w[i]=shared_w[threadIdx.x];
  }
}

__global__ void orcu_waxpykernel_1e5(const int n, double a, double* x, double* y, double* w) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_y[320];
  __shared__ double shared_x[320];
  __shared__ double shared_w[320];
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_y[threadIdx.x]=y[i];
    shared_x[threadIdx.x]=x[i];
    shared_w[threadIdx.x]=a*shared_x[threadIdx.x]+shared_y[threadIdx.x];
    w[i]=shared_w[threadIdx.x];
  }
}
#endif







/* ---------------------------------------------------------
// MAXPY
// Accumilation of AXPY on to x, currently only the manual
// version is available. This is a simple loop over y[i] vectors.
// The function checks each pairing for size mismatch and will
// synch the vectors to the GPU is necessary.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecMAXPY_SeqGPU"
PetscErrorCode VecMAXPY_SeqGPU(Vec x,PetscInt nv,const PetscScalar* alpha,Vec *y){
  /* x = x + sum(a[i]*y[i]) */
  PetscFunctionBegin;
  if(DEBUGVEC && VVERBOSE)printf("VecMAXPY_SeqGPU: alpha: %e\n",*alpha);
  PetscErrorCode ierr;
  PetscInt i;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=PETSC_NULL;
  dim3 dimGrid;  dim3 dimBlock;
  dimGrid.x=ceil((float)x->map->n/(float)AXPYTCOUNT);
  dimBlock.x=AXPYTCOUNT;
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Number of vectors in MAXPY: %d, blocks: %d, threads: %d\n",nv,dimGrid.x,dimBlock.x);
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"error in device malloc VecMAXPY_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"error in device memset VecMAXPY_SeqGPU");CHKERRQ(ierr);
  #endif
  if(xd->syncState==VEC_CPU){/* synch x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }

  for(i=0;i<nv;i++){
     if(y[i]->map->n!=x->map->n){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
    }
    yd=(Vec_SeqGPU*)y[i]->data;
    if(yd->syncState==VEC_CPU){/* synch x */
      ierr = VecCopyOverH2D(y[i],yd->cpuptr);CHKERRQ(ierr);
      xd->syncState=VEC_SYNCHED;
    }
    cudaDeviceSynchronize();
    if(alpha[i]==0){
      /* printf("no-op, continuing...\n"); */
      continue;
    }else{
      kernAXPY<<<dimGrid,dimBlock,0,xd->streamid>>>(xd->devptr,yd->devptr,alpha[i],y[i]->map->n);
      #if(DEBUGVEC)
        #if(VVERBOSE)
           if(!xd->devptr){
             printf("xd points to nothing.\n");
           }else if(!yd->devptr){
             printf("y[%d] points to nothing.\n",i);
           }
           printf("nv: %d, ylen[%d]: %d, alpha[%d]: %e, xlen: %d\n",nv,i,y[i]->map->n,i,alpha[i],x->map->n);
        #endif
        ierr = VecCheckCUDAError("kernel call to kernAXPY in VecMAXPY_SeqGPU");CHKERRQ(ierr);
      #endif
    }
  }
  xd->syncState=VEC_GPU;
  cudaDeviceSynchronize();
  PetscFunctionReturn(0);
}


/* ---------------------------------------------------------
// Device function XPY: x = x + y
// Labels are off, but this function is merely AXPY is alpha
// is equal to 1.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernXPY"
__global__ void  kernXPY(double* devY,double* devX, int vlen){
 /* y <- y + x */
  unsigned int tid = blockIdx.x*blockDim.x+threadIdx.x;
  if(tid<vlen){
    devY[tid]+=devX[tid];
  }
}



/* ---------------------------------------------------------
// MAXPY
// Accumilation of AXPY on to x, currently only the manual
// version is available. This is a simple loop over y[i] vectors.
// The function checks each pairing for size mismatch and will
// synch the vectors to the GPU is necessary.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecAXPY_SeqGPU"
PetscErrorCode VecAXPY_SeqGPU(Vec y,PetscScalar alpha,Vec x){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  if(x->map->n!=y->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  if(yd->syncState==VEC_CPU){/* synch y */
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  if(xd->syncState==VEC_CPU){/* synch x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  #if(DEBUGVEC && VVERBOSE)
     printf("VecAXPY_SeqGPU\n");
  #endif
  dim3 dimGrid,dimBlock;

#if(VMANXPY)
  dimGrid.x=ceil((float)x->map->n/(float)AXPYTCOUNT);
  dimBlock.x=AXPYTCOUNT;
  cudaDeviceSynchronize();
  kernAXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,alpha,y->map->n);
  #if(DEBUGVEC)
   ierr = VecCheckCUDAError("kernel call in VecAXPY_SeqGPU");CHKERRQ(ierr);
  #endif

#else
   if(x->map->n>=1e7){
     int nthreads=352;
     dimBlock.x=nthreads;
     dimGrid.x=14;
     cudaFuncSetCacheConfig(orcu_axpykernel_1e7,cudaFuncCachePreferL1);
     /*copy data from host to device*/
     orcu_axpykernel_1e7<<<dimGrid,dimBlock>>>(x->map->n,alpha,yd->devptr,xd->devptr);
     cudaError_t err=cudaGetLastError();
     if (cudaSuccess!=err) {
       printf("orcuda axpy_1e7 kernel error: %s\n",cudaGetErrorString(err));
       PetscFunctionReturn(PETSC_ERR_LIB);
     }
   }else if(x->map->n>=1e6){
       int nthreads=288;
       int nstreams=2;
       /*calculate device dimensions*/
       dimBlock.x=nthreads;
       dimGrid.x=28;
       /*create streams*/
       int istream, soffset;
       cudaStream_t stream[nstreams+1];
       for (istream=0; istream<=nstreams;istream++)cudaStreamCreate(&stream[istream]);
       int chunklen=x->map->n/nstreams;
       int chunkrem=x->map->n%nstreams;

       cudaFuncSetCacheConfig(orcu_axpykernel_1e6,cudaFuncCachePreferL1);
       /*invoke device kernel*/
       int blks4chunk=dimGrid.x/nstreams;
       if(dimGrid.x%nstreams!=0) blks4chunk++ ;
       for(istream=0; istream<nstreams; istream++ ) {
         soffset=istream*chunklen;
         orcu_axpykernel_1e6<<<blks4chunk,dimBlock,0,stream[istream]>>>
                             (chunklen,alpha,yd->devptr+soffset,xd->devptr+soffset);
       }
       if (chunkrem!=0) {
         soffset=istream*chunklen;
         orcu_axpykernel_1e6<<<blks4chunk,dimBlock,0,stream[istream]>>>
                             (chunkrem,alpha,yd->devptr+soffset,xd->devptr+soffset);
       }
       cudaDeviceSynchronize();
       for (istream=0; istream<=nstreams; istream++ )  cudaStreamDestroy(stream[istream]);

   }else{
     /*calculate device dimensions*/
     dimBlock.x=512;
     dimGrid.x=112;
     /*invoke device kernel*/
     orcu_axpykernel_1e5<<<dimGrid,dimBlock>>>(x->map->n,alpha,yd->devptr,xd->devptr);
   }

#endif

  yd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}



#if(VMANXPY)
__global__ void orcu_axpykernel_1e5(const int n, double a, double* y, double* x) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  for (int i=tid; i<=n-1; i+=gsize) {
    y[i]=y[i]+a*x[i];
  }
}

__global__ void orcu_axpykernel_1e6(const int n, double a, double* y, double* x) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_y[288];
  __shared__ double shared_x[288];
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_y[threadIdx.x]=y[i];
    shared_x[threadIdx.x]=x[i];
    shared_y[threadIdx.x]=shared_y[threadIdx.x]+a*shared_x[threadIdx.x];
    y[i]=shared_y[threadIdx.x];
  }
}

__global__ void orcu_axpykernel_1e7(const int n, double a, double* y, double* x) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_y[352];
  __shared__ double shared_x[352];
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_y[threadIdx.x]=y[i];
    shared_x[threadIdx.x]=x[i];
    shared_y[threadIdx.x]=shared_y[threadIdx.x]+a*shared_x[threadIdx.x];
    y[i]=shared_y[threadIdx.x];
  }
}
#endif




/* ---------------------------------------------------------
// Device kernel for AXPY
// Manual version
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernAXPY"
__global__ void  kernAXPY(double* devY,double* devX,double alpha, int vlen){
 /* y <- y + alpha*x */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  if(tid<vlen){
    devY[tid]+=alpha*devX[tid];
  }
}

#undef __FUNCT__
#define __FUNCT__ "kernSwap"
__global__ void kernSwap(int n, double *devx, double *devy, double *devscratch)
{
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  if (tid<n) {
    devscratch[tid] = devx[tid];
    devx[tid] = devy[tid];
    devy[tid] = devscratch[tid];
  }
}


/* ---------------------------------------------------------
// AXPBYPCZ: x = a*x + b*y + c*z
// Implemented, but currently only used by bicgs
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecAXPBYPCZ_SeqGPU"
PetscErrorCode VecAXPBYPCZ_SeqGPU(Vec x,PetscScalar alpha,PetscScalar beta,PetscScalar gamma,Vec y,Vec z){
  PetscFunctionBegin;
  #if(DEBUGVEC)
     PetscErrorCode ierr;
     #if(VVERBOSE)
        printf("Call to VecAXPBYPCZ_SeqGPU\n");
     #endif
  #endif
  Vec_SeqGPU* devX = (Vec_SeqGPU*)x->data;
  Vec_SeqGPU* devY = (Vec_SeqGPU*)y->data;
  Vec_SeqGPU* devZ = (Vec_SeqGPU*)z->data;
  double2 alphabeta;  alphabeta.x = alpha;  alphabeta.y = beta;
  dim3 dimGrid, dimBlock;
  dimGrid.x=ceil((float)x->map->n/(float)AXPBYPCZTCOUNT);
  dimBlock.x=AXPBYPCZTCOUNT;
  ccs[0]=cudaMemcpyToSymbol("dblScalar2Value",(void*)&alphabeta,sizeof(double2),0,cudaMemcpyHostToDevice);
  ccs[1]=cudaMemcpyToSymbol("dblScalarValue",(void*)&gamma,sizeof(double),0,cudaMemcpyHostToDevice);
  #if(DEBUGVEC)
   ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
   ierr = VecCheckCUDAStatus(ccs[1],"error in symbol copy to device");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();
  kernAXPBYPCZ<<<dimGrid,dimBlock>>>(devX->devptr,devY->devptr,devZ->devptr,devX->length);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAError("launch kernAXPBYPCZ");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


/* ---------------------------------------------------------
// Device kernel for AXPBYPCZ: x = a*x + b*y + c*z
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
extern __shared__ double sharedAXPBYPCZ[];
#undef __FUNCT__
#define __FUNCT__ "kernAXPBYPCZ"
__global__ void kernAXPBYPCZ(double* devX, double* devY, double* devZ, int* len){
  /* x <- alpha*x + beta*y + gamma*z */
  int localn = *len;
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  double work;
  if(tid<localn){
    /* do flops */
    if(dblScalarValue){
      work=dblScalarValue*devZ[tid];
    }else{
      work=0.;
    }

    if(dblScalar2Value.y){
      work+=dblScalar2Value.y*devY[tid];
    }
    if(dblScalar2Value.x){
      work+=dblScalar2Value.x*devX[tid];
    }
    /* write back */
    devX[tid]=work;
  }
  return;
}
/*---------------------------- end level 2 ------------------------------ */



/*------------------------- pointwise functions ------------------------- */

/* ---------------------------------------------------------
// Function which multiplies elementwise two vectors X and Y
// storing the result into a third vector W.
// Checks for size mismatch and synchs to the device if needed
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecPointwiseMult_SeqGPU"
PetscErrorCode VecPointwiseMult_SeqGPU(Vec w,Vec x,Vec y){
  PetscFunctionBegin;
  #if(DEBUGVEC && VERBOSE)
     printf("VecPointwiseMult_SeqGPU\n");
  #endif
  PetscErrorCode ierr;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)w->data;
  dim3 dimGrid, dimBlock;
  if(x->map->n!=y->map->n || w->map->n!=y->map->n || w->map->n!=x->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  if(yd->syncState==VEC_CPU){/* synch up y */
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  if(xd->syncState==VEC_CPU){/* synch up x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  dimGrid.x=ceil((float)y->map->n/(float)PMULTCOUNT);
  dimBlock.x=PMULTCOUNT;
  cudaDeviceSynchronize();
  kernPMULT<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAError("kernel call to kernPMULT");CHKERRQ(ierr);
  #endif
  wd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// Device kernel for pointwise multiply
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "kernPMULT"
__global__ void  kernPMULT(double* devY,double* devX, int* vlen, double* devW){
 /* w <- x./y */
  unsigned int tid = blockIdx.x*blockDim.x+threadIdx.x;
  if(tid<*vlen){
    devW[tid]=devX[tid]*devY[tid];
  }
}

/* ---------------------------------------------------------
// VecMaxPointwiseDivide_SeqGPU
// Function which calculates the elementwise division of vector
// X/Y, if one element is zero then X is just returned, the
// maximum value of all the resulting elements is returned to
// the host (max).
// The function implements a two stage reduction.
// Currently only the manual tuned version is implemented.
// Checks for size mismatch of arrays and synchs to device if
// needed.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecMaxPointwiseDivide_SeqGPU"
PetscErrorCode VecMaxPointwiseDivide_SeqGPU(Vec x,Vec y,PetscReal *max){
  PetscFunctionBegin;
  if(x->map->n!=y->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  #if(VTIMER)
   double start,finish,elapsed;
   static double mint,maxt=0.,cumt=0.,avg=0.;
   static int ccnt=0;
   start = vec_clock();
  #endif
  PetscErrorCode ierr;
  PetscScalar *devScratch;
  PetscInt i,chunks=0,segment,scratchsize;
  cudaStream_t* pwdstream;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid, dimBlock;
  /* Size of workload for the kernels */
  float threadscale = MAXMPLIER*CHUNKWIDTH;

  if(yd->syncState==VEC_CPU){/* synch up y */
    #if(DEBUGVEC && VVERBOSE)
       printf("yd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  if(xd->syncState==VEC_CPU){/* synch up x */
    #if(DEBUGVEC && VERBOSE)
       printf("xd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) / threadscale);
  pwdstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));
  /* make sure the segment size for each chunk is correct */
  if(chunks>1)segment = (int) threadscale;
  else segment = x->map->n;
  dimGrid.x=ceil(((float)segment)/(float)PDIVTCOUNT);
  dimBlock.x  = PDIVTCOUNT;
  #if(DEBUGVEC && VVERBOSE)
     printf("VecMaxPointwiseDivide_SeqGPU, chunks: %d segsize: %d\n",chunks,segment);
  #endif
  /* Divide up workload among streams and allocate scratch memory */
  scratchsize = chunks*dimGrid.x*sizeof(double);
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize);
  ccs[0] = cudaMemsetAsync(devScratch,0,scratchsize,yd->streamid);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecMPWD_SeqGPU");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();
  for(i=0;i<chunks;i++){/* first kernel */
    cudaStreamCreate(&(pwdstream[i]));
    /* Overlapping execution */
    kernMAXPDIV<<<dimGrid,dimBlock,0,pwdstream[i]>>>(xd->devptr,yd->devptr,
                                                     segment,
                                                     x->map->n,
                                                     i,
                                                     devScratch+i*dimGrid.x);
  }/* end for-loop */
  dimBlock.x  = PDIVTCOUNT2;
  scratchsize = chunks*dimGrid.x;
  cudaDeviceSynchronize();
  while(scratchsize>1){/* begin next reduction */
    dimGrid.x = ceil((float)scratchsize/(float)PDIVTCOUNT2);
    kernMAX<<<dimGrid,dimBlock>>>(scratchsize,devScratch);
    scratchsize = dimGrid.x;
  }

  /* copy back result */
  ccs[4]=cudaMemcpy(max,devScratch,sizeof(double),cudaMemcpyDeviceToHost);/* copy back */
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(ccs[4],"on cudaMemcpy(devScratch)");CHKERRQ(ierr);
  #endif

  /* Free temp resources */
  cms[3] = cudaFree(devScratch);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAStatus(cms[3],"on cudaFree(devScratch)");CHKERRQ(ierr);
  #endif
  for(i=0;i<chunks;i++) cudaStreamDestroy(pwdstream[i]);

  #if(DEBUGVEC && VVERBOSE)
    printf("max: %e\n",*max);
  #endif
  #if(VTIMER)
    finish=vec_clock();
    elapsed=finish-start;
    cumt+=elapsed;
    if(!ccnt++){
      maxt=mint=avg=elapsed;
    }else{
      maxt=elapsed>maxt?elapsed:maxt;
      mint=elapsed<mint?elapsed:mint;
      avg=cumt/ccnt;
      printf("VecMAXPWD_SeqGPU calls: %d, max: %e, min: %e, average: %e\n",
               ccnt,maxt,mint,avg);
    }
  #endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernMAX"
__global__ void  kernMAX(int n,double* scratch){
  unsigned int tid = blockDim.x*blockIdx.x+threadIdx.x;
  __shared__ double chunk[PDIVTCOUNT2];
  double mymax;
  mymax=(tid<n)?scratch[tid]:0.;
  if(threadIdx.x>127)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<128)mymax=fmax(mymax,chunk[threadIdx.x+128]);
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<64)mymax=fmax(mymax,chunk[threadIdx.x+64]);
  else return;
  chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<32)warpMaxReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}

#undef __FUNCT__
#define __FUNCT__ "kernMAXPDIV"
__global__ void  kernMAXPDIV(double* devX,double* devY, int segmentsize,
                             int n,int offset,double* scratch){
 /* w <- max(abs(x./y)) */
  __shared__ double chunk[PDIVTCOUNT];
  double mymax;
  unsigned int item = segmentsize*offset+blockDim.x*blockIdx.x+threadIdx.x;
  if(item<n){
    mymax=devY[item];
    if(mymax!=0.)mymax=fabs(devX[item]/mymax);//reusing register
    else mymax=fabs(devX[item]);
  }else{
    mymax=0.;
  }
  if(threadIdx.x>127)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<128)mymax=fmax(mymax,chunk[threadIdx.x+128]);
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<64)mymax=fmax(mymax,chunk[threadIdx.x+64]);
  else return;
  chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<32) warpReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}




/* ---------------------------------------------------------
// VecPointwiseDivide_SeqGPU
// Function which calculates the elementwise division of vector
// X/Y and stores the results in a third array. If one element 
// is zero then X is just returned. Currently only the manual
// tuned version is implemented.
// Checks for size mismatch of arrays and synchs to device if
// needed.
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecPointwiseDivide_SeqGPU"
PetscErrorCode VecPointwiseDivide_SeqGPU(Vec w,Vec x,Vec y){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)w->data;
  dim3 dimGrid, dimBlock;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecPointwiseDivide_SeqGPU\n");
  #endif
  if(x->map->n!=y->map->n || w->map->n!=y->map->n || w->map->n!=x->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  if(yd->syncState==VEC_CPU){/* synch up y */
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  if(xd->syncState==VEC_CPU){/* synch up x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  dimGrid.x=ceil((float)y->map->n/(float)PDIVTCOUNT);
  dimBlock.x=PDIVTCOUNT;
  cudaDeviceSynchronize();
  kernPDIV<<<dimGrid,dimBlock,2*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAError("kernel call to kernPDIV");CHKERRQ(ierr);
  #endif
  wd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}
extern __shared__ double sharedPDIV[];
#undef __FUNCT__
#define __FUNCT__ "kernPDIV"
__global__ void  kernPDIV(double* devY,double* devX, int* vlen, double* devW){
  /* w <- x./y */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  double* chunkX = sharedPDIV;
  double* chunkY = sharedPDIV + blockDim.x;
  double work;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    if(chunkX[threadIdx.x]*chunkY[threadIdx.x]!=0){
      work=chunkX[threadIdx.x]/chunkY[threadIdx.x];
    }else{
      work=0;
    }
    devW[tid]=work;
  }
}

/*--------------------------- end pointwise ---------------------------- */

/*--------------------------- norm functions --------------------------- */

/* ---------------------------------------------------------
// VecDotNorm2_SeqGPU
// Simple wrapper function for two calls
// Never seen this function called
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecDotNorm2_SeqGPU"
PetscErrorCode VecDotNorm2_SeqGPU(Vec s, Vec t, PetscScalar *dp, PetscScalar *nm){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecDotNorm2_SeqGPU\n");
  #endif
  ierr = VecDot(s,t,dp); CHKERRQ(ierr);
  ierr = VecNorm(t,NORM_2,nm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------
// VecNorm_SeqGPU
// Function which computes the norm of a vector
// Currently implements two norm types: infinity norm and norm2
// Norm2 is implemented as a manually tuned kernel as well a
// Orio tuned kernels for three vector size ranges
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecNorm_SeqGPU"
PetscErrorCode VecNorm_SeqGPU(Vec x,NormType type,PetscReal* z){
  /* NormType: NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4 */
  /* dealing with NORM_2 and NORM_INFINITY for now... */
  PetscFunctionBegin;
  #if(VTIMER)
    double start,finish,elapsed;
    static double mint,maxt=0.,cumt=0.,avg=0.;
    static int ccnt=0;
    start = vec_clock();
  #endif
  PetscErrorCode ierr;
  double *devScratch,zhost;
  PetscInt i,chunks=0,segment,scratchsize;
  cudaStream_t* nrmstream;
  dim3 dimGrid, dimBlock;
  /* defining per-stream work load */
  float threadscale = NRMMPLIER*CHUNKWIDTH;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  if(xd->syncState==VEC_CPU){
    #if(DEBUGVEC && VVERBOSE)
       printf("xd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }

  if(type==NORM_INFINITY){
    printf("Infinity NORM.\n");
    /* figure out how many chunks will be needed */
    chunks = ceil( ((float)x->map->n) / threadscale);
    nrmstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));
    /* make sure the segment size for each chunk is correct */
    if(chunks>1) segment = (int) threadscale;
    else segment = x->map->n;
    dimGrid.x=ceil(((float)segment)/(float)THRNRMCNT);
    dimBlock.x  = THRNRMCNT;
    #if(DEBUGVEC && VVERBOSE)
      printf("Call to VecNorm_SeqGPU chunks: %d segsize: %d\n",chunks,segment);
    #endif
    /* allocate gridwide scratch array */
    scratchsize = chunks*dimGrid.x;
    cms[0] = cudaMalloc((void**)&devScratch,scratchsize*sizeof(double));
    #if(DEBUGVEC)
      #if(VVERBOSE)
         printf("NORM2: chunks: %d, seg: %d, blks: %d\n",chunks,segment,dimGrid.x);
      #endif
      ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecNorm_SeqGPU"); CHKERRQ(ierr);
    //ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();
    for(i=0;i<chunks;i++){/* streaming async kernel calls */
      cudaStreamCreate(&(nrmstream[i]));
      /* Overlapping execution */
      kernInfNorm<<<dimGrid,dimBlock,0,nrmstream[i]>>>(xd->devptr,segment,x->map->n,i,
                                                       devScratch+i*dimGrid.x);
    }/* end for-loop */
    dimBlock.x  = THRNRMCNT2;
    cudaDeviceSynchronize();
    while(scratchsize>1){/* begin next reduction */
      /* printf("Seconds stage reduction.\n"); */
      dimGrid.x = ceil((float)scratchsize/(float)THRNRMCNT2);
      kernRedInfNorm<<<dimGrid,dimBlock>>>(scratchsize,devScratch);
      scratchsize = dimGrid.x;
    }
  }else{/* NORM2 etc... */

#if(VMANNRM)
  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) / threadscale);
  nrmstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));
  /* make sure the segment size for each chunk is correct */
  if(chunks>1) segment = (int) threadscale;
  else segment = x->map->n;
  dimGrid.x=ceil(((float)segment)/(float)THRNRMCNT);
  dimBlock.x  = THRNRMCNT;
#if(DEBUGVEC && VVERBOSE)
  printf("Call to VecNorm_SeqGPU chunks: %d segsize: %d\n",chunks,segment);
#endif
  /* allocate gridwide scratch array */
  scratchsize = chunks*dimGrid.x;
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize*sizeof(double));
  #if(DEBUGVEC)
    #if(VVERBOSE)
       printf("NORM2: chunks: %d, seg: %d, blks: %d\n",chunks,segment,dimGrid.x);
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecNorm_SeqGPU"); CHKERRQ(ierr);
    //ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecNorm_SeqGPU");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();
    for(i=0;i<chunks;i++){/* streaming async kernel calls */
      cudaStreamCreate(&(nrmstream[i]));
      /* Overlapping execution */
      kernNorm2<<<dimGrid,dimBlock,0,nrmstream[i]>>>(xd->devptr,segment,x->map->n,i,
                                                     devScratch+i*dimGrid.x);
    }/* end for-loop */
    dimBlock.x  = THRNRMCNT2;
    cudaDeviceSynchronize();
    while(scratchsize>1){/* begin next reduction */
      /* printf("Seconds stage reduction.\n"); */
      dimGrid.x = ceil((float)scratchsize/(float)THRNRMCNT2);
      kernRedNorm<<<dimGrid,dimBlock>>>(scratchsize,devScratch);
      scratchsize = dimGrid.x;
    }

    cudaDeviceSynchronize();
    for(i=0;i<chunks;i++) cudaStreamDestroy(nrmstream[i]);
    free(nrmstream);
#else

    if(x->map->n>=1e7){
      cudaFuncSetCacheConfig(orcu_norm2kernel_1e7,cudaFuncCachePreferL1);
      int nthreads=512;
      int nstreams=2;
      /*calculate device dimensions*/
      dim3 dimGrid, dimBlock;
      dimBlock.x=nthreads;
      dimGrid.x=112;
      cudaMalloc((void**)&devScratch,(dimGrid.x+1)*sizeof(double));
      /*create streams*/
      int istream, soffset, boffset;
      cudaStream_t stream[nstreams+1];
      for (istream=0; istream<=nstreams;istream++)cudaStreamCreate(&stream[istream]);
      int chunklen=x->map->n/nstreams;
      int chunkrem=x->map->n%nstreams;
   
      /*invoke device kernel*/
      int blks4chunk=dimGrid.x/nstreams;
      if(dimGrid.x%nstreams!=0)blks4chunk++ ;
      int blks4chunks=blks4chunk*nstreams;
      for(istream=0; istream<nstreams; istream++){
        soffset=istream*chunklen;
        boffset=istream*blks4chunk;
        orcu_norm2kernel_1e7<<<blks4chunk,dimBlock,0,stream[istream]>>>
                           (chunklen,xd->devptr+soffset,devScratch+boffset);
      }
      if (chunkrem!=0) {
        soffset=istream*chunklen;
        boffset=istream*blks4chunk;
        orcu_norm2kernel_1e7<<<blks4chunk,dimBlock,0,stream[istream]>>>
                           (chunkrem,xd->devptr+soffset,devScratch+boffset);
        blks4chunks++ ;
      }
      int orcu_blks=blks4chunks;
      int orcu_n;
      while (orcu_blks>1) {
        orcu_n=orcu_blks;
        orcu_blks=(orcu_blks+511)/512;
        orcu_norm2blksum_1e7<<<orcu_blks,512>>>(orcu_n,devScratch);
      }
      for (istream=0; istream<=nstreams; istream++)cudaStreamDestroy(stream[istream]);
    }else if(x->map->n>=1e6){
      /*calculate device dimensions*/
      dimBlock.x=228;
      dimGrid.x=56;
      cudaMalloc((void**)&devScratch,(dimGrid.x+1)*sizeof(double));
      orcu_norm2kernel_1e6<<<dimGrid,dimBlock>>>(x->map->n,xd->devptr,devScratch);
      int orcu_blks=dimGrid.x;
      int orcu_n;
      while (orcu_blks>1) {
        orcu_n=orcu_blks;
        orcu_blks=(orcu_blks+227)/228;
        orcu_norm2blksum_1e6<<<orcu_blks,228>>>(orcu_n,devScratch);
      }
    }else{
      /*calculate device dimensions*/
      dimBlock.x=128;
      dimGrid.x=112;
      cudaMalloc((void**)&devScratch,(dimGrid.x+1)*sizeof(double));
      orcu_norm2kernel_1e5<<<dimGrid,dimBlock>>>(x->map->n,xd->devptr,devScratch);
      int orcu_blks=dimGrid.x;
      int orcu_n;
      while (orcu_blks>1) {
        orcu_n=orcu_blks;
        orcu_blks=(orcu_blks+127)/128;
        orcu_norm2blksum_1e5<<<orcu_blks,128>>>(orcu_n,devScratch);
      }
    }
#endif /* end VMANNRM norm2 */

    ccs[4]=cudaMemcpy(&zhost,devScratch,sizeof(double),cudaMemcpyDeviceToHost);/* copy back */
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(ccs[4],"on cudaMemcpy(devScratch)");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    *z = PetscSqrtScalar(zhost);
  }/* end NORMTYPE if */


  /* clean up resources */

  cms[3] = cudaFree(devScratch);
  #if(DEBUGVEC)
   #if(VVERBOSE)
        printf("Znorm: %e\n",*z);
   #endif
     ierr = VecCheckCUDAStatus(cms[3],"on cudaFree(devScratch)");CHKERRQ(ierr);
  #endif
  #if(VTIMER)
    finish=vec_clock();
    elapsed=finish-start;
    cumt+=elapsed;
    if(!ccnt++){
      maxt=mint=avg=elapsed;
    }else{
      maxt=elapsed>maxt?elapsed:maxt;
      mint=elapsed<mint?elapsed:mint;
      avg=cumt/ccnt;
      if(!(ccnt%(ITSHOW/2))){
        printf("VecNorm_SeqGPU calls: %d, max: %e, min: %e, average: %e\n",
               ccnt,maxt,mint,avg);
      }
    }
  #endif
  PetscFunctionReturn(0);
}

/*-------------- Device kernels for infinite norm ----------------*/
#undef __FUNCT__
#define __FUNCT__ "kernRedInfNorm"
__global__ void kernRedInfNorm(int n,double* scratch){/* reduction kernel */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunk[THRNRMCNT2];
  double mymax;
  mymax=(tid<n)?scratch[tid]:0.;
  if(threadIdx.x>127)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<128)mymax=fmax(mymax,chunk[threadIdx.x+128]);
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<64)mymax=fmax(mymax,chunk[threadIdx.x+64]);
  else return;
  chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<32)warpMaxReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}

#undef __FUNCT__
#define __FUNCT__ "kernInfNorm"
__global__ void kernInfNorm(double* devX,int segmentsize,
                        int arrsize, int offset, double* scratch){
  __shared__ double chunk[THRNRMCNT];
  unsigned int item = segmentsize*offset+blockDim.x*blockIdx.x+threadIdx.x;
  double mymax=0.;
  mymax=(item<arrsize)?fabs(devX[item]):0.;

  if(threadIdx.x>127)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<128)mymax=fmax(mymax,chunk[threadIdx.x+128]);
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<64)mymax=fmax(mymax,chunk[threadIdx.x+64]);
  else return;
  chunk[threadIdx.x]=mymax;
  __syncthreads();
  if(threadIdx.x<32) warpMaxReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}
/*---------------------------------------------------------*/

/*-------------- Device kernels for norm2 ----------------*/
#undef __FUNCT__
#define __FUNCT__ "kernRedNorm"
__global__ void kernRedNorm(int n,double* scratch){/* reduction kernel */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunk[THRNRMCNT2];
  double mysum;
  mysum=(tid<n)?scratch[tid]:0.;
  if(threadIdx.x>127)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<128)mysum+=chunk[threadIdx.x+128];
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<64)mysum+=chunk[threadIdx.x+64];
  else return;
  chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<32)warpReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}

#undef __FUNCT__
#define __FUNCT__ "kernNorm2"
__global__ void kernNorm2(double* devX,int segmentsize,
                        int arrsize, int offset, double* scratch){
  __shared__ double chunk[THRNRMCNT];
  unsigned int item = segmentsize*offset+blockDim.x*blockIdx.x+threadIdx.x;
  double mysum=0.;
  mysum=(item<arrsize)?devX[item]:0.;
  mysum*=mysum;
  if(threadIdx.x>127)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<128)mysum+=chunk[threadIdx.x+128];
  else return;
  if(threadIdx.x>63)chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<64)mysum+=chunk[threadIdx.x+64];
  else return;
  chunk[threadIdx.x]=mysum;
  __syncthreads();
  if(threadIdx.x<32) warpReduce(chunk,threadIdx.x);
  else return;
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunk[0];
  }else return;
}
/*---------------------------------------------------------*/


/*------------------- Orio Norm2 kernels ------------------*/
__global__ void orcu_norm2kernel_1e5(const int n, double* x, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  double orcu_var10241=0;
  for (int i=tid; i<=n-1; i+=gsize) {
    orcu_var10241=orcu_var10241+x[i]*x[i];
  }
  /*reduce single-thread results within a block*/
  __shared__ double orcu_vec10242[128];
  orcu_vec10242[threadIdx.x]=orcu_var10241;
  __syncthreads();
  if (threadIdx.x<64) orcu_vec10242[threadIdx.x]+=orcu_vec10242[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) orcu_warpReduce64(threadIdx.x,orcu_vec10242);
  __syncthreads();
  if (threadIdx.x==0) reducts[blockIdx.x]=orcu_vec10242[0];
}
__global__ void orcu_norm2blksum_1e5(int orcu_n, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double orcu_vec10242[128];

  if (tid<orcu_n) orcu_vec10242[threadIdx.x]=reducts[tid];
  else  orcu_vec10242[threadIdx.x]=0;
  __syncthreads();
  if (threadIdx.x<64) orcu_vec10242[threadIdx.x]+=orcu_vec10242[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) orcu_warpReduce64(threadIdx.x,orcu_vec10242);
  __syncthreads();
  if (threadIdx.x==0) reducts[blockIdx.x]=orcu_vec10242[0];
}

__global__ void orcu_norm2kernel_1e6(const int n, double* x, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  double orcu_var20485=0;
  for (int i=tid; i<=n-1; i+=gsize) {
    orcu_var20485=orcu_var20485+x[i]*x[i];
  }
  /*reduce single-thread results within a block*/
  __shared__ double orcu_vec20486[288];
  orcu_vec20486[threadIdx.x]=orcu_var20485;
  __syncthreads();
  if (threadIdx.x<128) 
    orcu_vec20486[threadIdx.x]+=orcu_vec20486[threadIdx.x+128];
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec20486[threadIdx.x]+=orcu_vec20486[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec20486);
  if (threadIdx.x>=256&&threadIdx.x<272) 
    orcu_warpReduce32(threadIdx.x,orcu_vec20486);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec20486[0]+orcu_vec20486[256];
}
__global__ void orcu_norm2blksum_1e6(int orcu_n, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double orcu_vec20486[288];
  if (tid<orcu_n) 
    orcu_vec20486[threadIdx.x]=reducts[tid];
  else 
    orcu_vec20486[threadIdx.x]=0;
  __syncthreads();
  if (threadIdx.x<128) 
    orcu_vec20486[threadIdx.x]+=orcu_vec20486[threadIdx.x+128];
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec20486[threadIdx.x]+=orcu_vec20486[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec20486);
  if (threadIdx.x>=256&&threadIdx.x<272) 
    orcu_warpReduce32(threadIdx.x,orcu_vec20486);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec20486[0]+orcu_vec20486[256];
}

__global__ void orcu_norm2kernel_1e7(const int n, double* x, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  const int gsize=gridDim.x*blockDim.x;
  __shared__ double shared_x[512];
  double orcu_var30729=0;
  for (int i=tid; i<=n-1; i+=gsize) {
    shared_x[threadIdx.x]=x[i];
    orcu_var30729=orcu_var30729+shared_x[threadIdx.x]*shared_x[threadIdx.x];
  }
  /*reduce single-thread results within a block*/
  __shared__ double orcu_vec30730[512];
  orcu_vec30730[threadIdx.x]=orcu_var30729;
  __syncthreads();
  if (threadIdx.x<256) 
    orcu_vec30730[threadIdx.x]+=orcu_vec30730[threadIdx.x+256];
  __syncthreads();
  if (threadIdx.x<128) 
    orcu_vec30730[threadIdx.x]+=orcu_vec30730[threadIdx.x+128];
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec30730[threadIdx.x]+=orcu_vec30730[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec30730);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec30730[0];
}
__global__ void orcu_norm2blksum_1e7(int orcu_n, double* reducts) {
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double orcu_vec30730[512];
  if (tid<orcu_n) 
    orcu_vec30730[threadIdx.x]=reducts[tid];
  else 
    orcu_vec30730[threadIdx.x]=0;
  __syncthreads();
  if (threadIdx.x<256) 
    orcu_vec30730[threadIdx.x]+=orcu_vec30730[threadIdx.x+256];
  __syncthreads();
  if (threadIdx.x<128) 
    orcu_vec30730[threadIdx.x]+=orcu_vec30730[threadIdx.x+128];
  __syncthreads();
  if (threadIdx.x<64) 
    orcu_vec30730[threadIdx.x]+=orcu_vec30730[threadIdx.x+64];
  __syncthreads();
  if (threadIdx.x<32) 
    orcu_warpReduce64(threadIdx.x,orcu_vec30730);
  __syncthreads();
  if (threadIdx.x==0) 
    reducts[blockIdx.x]=orcu_vec30730[0];
}
/*---------------------------------------------------------*/
/* --------------------- end norms ----------------------- */


/* ---------------------------------------------------------
// VecGetArray_SeqGPU
// Grabs the pointer to the cpu memory, if necessary copies
// that array up from the device
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecGetArray_SeqGPU"
PetscErrorCode VecGetArray_SeqGPU(Vec v,PetscScalar **a){
#ifdef PETSC_USE_DEBUG
  PetscInt flg1=0,flg2=0,flg3=0,flg4=0;
#endif
  PetscErrorCode ierr;
  Vec_SeqGPU *vd=(Vec_SeqGPU*)v->data;

  PetscFunctionBegin;

  if(vd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"in VecGetArray_SeqGPU the vector has not been created.");
  }
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecGetArray_SeqGPU\n");
  #endif
#ifdef PETSC_USE_DEBUG
  /* PETSc in debug mode uses a macro to VecValidValues
     to test values before trying to use the vector. In order to
     prevent these checks (which all require cudamemcpy), 
     the stack is checked to make sure it's a real need for the values
  */
   MyPetscStackCheckByName("DMDAVecGetArray",flg1);
   MyPetscStackCheckByName("DMGlobalToLocalBegin",flg2);
   MyPetscStackCheckByName("SNESDefaultComputeJacobian",flg3);
   MyPetscStackCheckByName("DMComputeJacobianDefault",flg4);
  if(flg1 || flg2 || flg3 || flg4 ){
    if(vd->syncState==VEC_GPU){
      ierr = VecCopyOverD2H(v,vd->cpuptr); CHKERRQ(ierr);
    }
    vd->syncState = VEC_CPU;
  }
#else
  if(vd->syncState==VEC_GPU){
    ierr = VecCopyOverD2H(v,vd->cpuptr); CHKERRQ(ierr);
  }
  vd->syncState = VEC_CPU;
#endif
  cudaDeviceSynchronize();
  *a=vd->cpuptr;
  PetscFunctionReturn(0);
}


/* ---------------------------------------------------------
// VecRestoreArray_SeqGPU
// Returns data back to the vector type and copying back
// the memory to device only if necessary
// written by: dlowell ANL-MCS
// --------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecRestoreArray_SeqGPU"
PetscErrorCode VecRestoreArray_SeqGPU(Vec v,PetscScalar **a){
  PetscErrorCode ierr;
  Vec_SeqGPU *vd=(Vec_SeqGPU*)v->data;
#ifdef PETSC_USE_DEBUG
  PetscInt flg1=0,flg2=0,flg3=0;
#endif
  PetscFunctionBegin;

#ifdef PETSC_USE_DEBUG
  /* PETSc in debug mode uses a macro to VecValidValues
     to test values before trying to use the vector. In order to
     prevent these checks (which all require cudamemcpy), 
     the stack is checked to make sure it's a real need for the values
  */
  MyPetscStackCheckByName("VecRestoreArrayRead",flg1);
  MyPetscStackCheckByName("DMDAVecRestoreArray",flg2);
  MyPetscStackCheckByName("DMGlobalToLocalBegin",flg3);
  if(vd->syncState==VEC_CPU||(!flg1||flg2||flg3)){
    if(a){
      ierr = VecCopyOverH2D(v,*a);CHKERRQ(ierr);
      vd->syncState=VEC_GPU;
    }else{
      ierr = VecCopyOverH2D(v,vd->cpuptr);CHKERRQ(ierr);
      vd->syncState=VEC_SYNCHED;
    }
  }
#else
  if(a){
    ierr = VecCopyOverH2D(v,*a);CHKERRQ(ierr);
    vd->syncState=VEC_GPU;
  }else{
    ierr = VecCopyOverH2D(v,vd->cpuptr);CHKERRQ(ierr);
    vd->syncState=VEC_SYNCHED;
  }
#endif
  cudaDeviceSynchronize();
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCreateSeqGPU"
PetscErrorCode  VecCreateSeqGPU(MPI_Comm comm,PetscInt n,Vec *v){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = VecCreate(comm,v);CHKERRQ(ierr);
  ierr = VecSetSizes(*v,n,n);CHKERRQ(ierr);
  ierr = VecSetType(*v,VECSEQGPU);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCopy_SeqGPU"
PetscErrorCode VecCopy_SeqGPU(Vec s,Vec d){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU *sd=(Vec_SeqGPU*)s->data;
  Vec_SeqGPU *dd=(Vec_SeqGPU*)d->data;
  if(d->map->n!=s->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  if(dd->syncState==VEC_UNALLOC){
     SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Destination vector unalloced.");
  }
  if(sd->syncState==VEC_ALLOC){
      PetscFunctionReturn(0);/* nothing to do */
  }
  if(sd->syncState==VEC_CPU){
    ierr = PetscMemcpy((void*)dd->cpuptr,(void*)sd->cpuptr,s->map->n*sizeof(double));CHKERRQ(ierr);
    dd->syncState = VEC_CPU;
    PetscFunctionReturn(0);
  }
  ierr = VecCopyOverDevice(d,s); CHKERRQ(ierr);
  ierr = cudaDeviceSynchronize();
  dd->syncState=VEC_GPU;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecCopy_SeqGPU\n");
  #endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecSwap_SeqGPU"
PetscErrorCode VecSwap_SeqGPU(Vec xin,Vec yin){
  PetscErrorCode ierr;
  Vec_SeqGPU *x=(Vec_SeqGPU*)xin->data;
  Vec_SeqGPU *y=(Vec_SeqGPU*)yin->data;
  PetscBLASInt one=1,bn=PetscBLASIntCast(xin->map->n);
  PetscScalar *devScratch;
  dim3 dimGrid,dimBlock;
  PetscFunctionBegin;
  dimGrid.x=ceil((float)xin->map->n/(float)AXPYTCOUNT);
  dimBlock.x=AXPYTCOUNT;
  if (xin != yin) {
    if ((x->syncState == VEC_GPU && y->syncState == VEC_GPU) ||
        (x->syncState == VEC_SYNCHED && y->syncState == VEC_GPU) ||
        (x->syncState == VEC_GPU && y->syncState == VEC_SYNCHED)) {

      /* If both vectors current on GPU */
      cms[0] = cudaMalloc((void**)&devScratch,xin->map->n*sizeof(PetscScalar));
      ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecSwap_SeqGPU");CHKERRQ(ierr);
      kernSwap<<<dimGrid,dimBlock>>>(xin->map->n,x->devptr,y->devptr,devScratch);
      x->syncState = VEC_GPU;
      y->syncState = VEC_GPU;
      cms[1] = cudaFree(devScratch);
      ierr = VecCheckCUDAStatus(cms[1],"on VecSwap(devScratch)");CHKERRQ(ierr);
      cudaDeviceSynchronize();

    } else {
      if (y->syncState != VEC_CPU) {
        ierr = VecCopyOverD2H(yin,y->cpuptr); CHKERRQ(ierr);
        cudaDeviceSynchronize();
      } else if (x->syncState != VEC_CPU) {
        ierr = VecCopyOverD2H(xin,x->cpuptr); CHKERRQ(ierr);
        cudaDeviceSynchronize();
      }
      BLASswap_(&bn,x->cpuptr,&one,y->cpuptr,&one);
      x->syncState = VEC_CPU;
      y->syncState = VEC_CPU;

    }
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecDuplicate_SeqGPU"
PetscErrorCode VecDuplicate_SeqGPU(Vec win,Vec *V){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecDuplicate_SeqGPU\n");
  #endif
  ierr = VecCreate(((PetscObject)win)->comm,V);CHKERRQ(ierr);
  ierr = VecSetType(*V,VECSEQGPU);CHKERRQ(ierr);
  ierr = PetscObjectSetPrecision((PetscObject)*V,((PetscObject)win)->precision);CHKERRQ(ierr);
  ierr = VecSetSizes(*V,win->map->n,win->map->N);CHKERRQ(ierr);
  ierr = PetscLayoutReference(win->map,&(*V)->map);CHKERRQ(ierr);
  ierr = PetscOListDuplicate(((PetscObject)win)->olist,&((PetscObject)(*V))->olist);CHKERRQ(ierr);
  ierr = PetscFListDuplicate(((PetscObject)win)->qlist,&((PetscObject)(*V))->qlist);CHKERRQ(ierr);
  (*V)->stash.ignorenegidx = win->stash.ignorenegidx;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecDuplicateVecs_SeqGPU"
PetscErrorCode VecDuplicateVecs_SeqGPU(Vec vin, PetscInt m, Vec **Vlist){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i=0;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecDuplicateVecs_SeqGPU\n"); 
  #endif
  PetscValidHeaderSpecific(vin,VEC_CLASSID,1);
  PetscValidPointer(Vlist,3);
  if (m <= 0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"m must be > 0: m = %D",m);
  ierr = PetscMalloc(m*sizeof(Vec),Vlist);CHKERRQ(ierr);
  for(i=0;i<m;i++){
    ierr = VecDuplicate_SeqGPU(vin,*Vlist+i);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecPlaceArray_SeqGPU"
PetscErrorCode  VecPlaceArray_SeqGPU(Vec x,const PetscScalar* array){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecPlaceArray_SeqGPU\n"); 
  #endif
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"in VecPlaceArray_SeqGPU the vector has not been created.");
  }
  if(xd->unplacedarray){
     SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,
       "VecPlaceArray() was already called on this vector, without a call to VecResetArray()");
  }
  if(xd->syncState==VEC_GPU){/* assuming there is a logical reason for this copy up */
    ierr = VecCopyOverD2H(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  xd->unplacedarray=xd->cpuptr;
  xd->cpuptr=(PetscScalar*)array;
  ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
  xd->syncState=VEC_SYNCHED;
  cudaDeviceSynchronize();
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecResetArray_SeqGPU"
PetscErrorCode  VecResetArray_SeqGPU(Vec x){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecResetArray_SeqGPU\n"); 
  #endif
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"in VecResetArray_SeqGPU the vector has not been created.");
  }
  if(xd->cpuptr){
    ierr = PetscFree(xd->cpuptr);CHKERRQ(ierr);
  }
  xd->cpuptr=xd->unplacedarray;
  xd->unplacedarray=PETSC_NULL;
  ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
  xd->syncState=VEC_SYNCHED;
  cudaDeviceSynchronize();
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecReplaceArray_SeqGPU"
PetscErrorCode  VecReplaceArray_SeqGPU(Vec x,const PetscScalar* array){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecReplaceArray_SeqGPU\n"); 
  #endif
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"in VecResetArray_SeqGPU the vector has not been created.");
  }
  if(xd->cpuptr){
    ierr = PetscFree(xd->cpuptr);CHKERRQ(ierr);
  }
  xd->cpuptr=(PetscScalar*)array;
  ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
  xd->syncState=VEC_SYNCHED;
  cudaDeviceSynchronize();
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecMax_SeqGPU"
PetscErrorCode VecMax_SeqGPU(Vec xin,PetscInt* idx,PetscReal * z)
{
  PetscInt          i,j=0,n = xin->map->n;
  PetscReal         max,tmp;
  Vec_SeqGPU*       xd = (Vec_SeqGPU*)xin->data;
  PetscErrorCode    ierr;

  PetscFunctionBegin;

  if(xd->syncState==VEC_GPU){
    ierr = VecCopyOverD2H(xin,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
    cudaDeviceSynchronize();
  }

  if (!n) {
    max = PETSC_MIN_REAL;
    j   = -1;
  } else {
#if defined(PETSC_USE_COMPLEX)
      max = PetscRealPart(xd->cpuptr[0]); j = 0;
#else
      max = xd->cpuptr[0]; j = 0;
#endif
    for (i=1; i<n; i++) {
#if defined(PETSC_USE_COMPLEX)
      if ((tmp = PetscRealPart(xd->cpuptr[i])) > max) { j = i; max = tmp;}
#else
      if ((tmp = xd->cpuptr[i]) > max) { j = i; max = tmp; }
#endif
    }
  }
  *z   = max;
  if (idx) *idx = j;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecMin_SeqGPU"
PetscErrorCode VecMin_SeqGPU(Vec xin,PetscInt* idx,PetscReal * z)
{
  PetscInt          i,j=0,n = xin->map->n;
  PetscReal         min,tmp;
  Vec_SeqGPU*       xd = (Vec_SeqGPU*)xin->data;
  PetscErrorCode    ierr;

  PetscFunctionBegin;

  if(xd->syncState==VEC_GPU){
    ierr = VecCopyOverD2H(xin,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }

  if (!n) {
    min = PETSC_MAX_REAL;
    j   = -1;
  } else {
#if defined(PETSC_USE_COMPLEX)
      min = PetscRealPart(xd->cpuptr[0]); j = 0;
#else
      min = xd->cpuptr[0]; j = 0;
#endif
    for (i=1; i<n; i++) {
#if defined(PETSC_USE_COMPLEX)
      if ((tmp = PetscRealPart(xd->cpuptr[i])) < min) { j = i; min = tmp;}
#else
      if ((tmp = xd->cpuptr[i]) < min) { j = i; min = tmp; }
#endif
    }
  }
  *z   = min;
  if (idx) *idx = j;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PinnedMalloc"
PetscErrorCode  PinnedMalloc(PetscScalar** x,PetscInt n){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscInfo1(0,"Allocating %d bytes on GPU\n",n); CHKERRQ(ierr);
  ierr=VecCheckCUDAStatus(cms[0],"before PinnedMalloc");CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF,"Allocating %d bytes on GPU\n",n); CHKERRQ(ierr);
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to PinnedMalloc\n"); 
  #endif
  cms[0]=cudaHostAlloc((void**)x,n,0);

  ierr=VecCheckCUDAStatus(cms[0],"in PinnedMalloc");CHKERRQ(ierr);

  //SETERRQ1(PETSC_COMM_SELF,0,"Error Allocating Memory -- %d bytes requested",n);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PinnedFree"
PetscErrorCode  PinnedFree(PetscScalar* x){
  PetscFunctionBegin;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to PinnedFree\n"); 
  #endif
  cms[0]=cudaFreeHost(x);
  #if(DEBUGVEC)
    PetscErrorCode ierr;
    ierr=VecCheckCUDAStatus(cms[0],"in PinnedFree");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCreate_SeqGPU"
PetscErrorCode  VecCreate_SeqGPU(Vec V){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscMPIInt    size;
  Vec_SeqGPU* seqgpu=PETSC_NULL;
  ierr = PetscMalloc(sizeof(Vec_SeqGPU),&seqgpu);
  V->data=(void*)seqgpu;
  ierr = MPI_Comm_size(((PetscObject)V)->comm,&size);CHKERRQ(ierr);
  if  (size > 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Cannot create VECSEQGPU on more than one process");

  if (V->map->bs == -1) V->map->bs = 1;
  ierr = PetscLayoutSetUp(V->map);CHKERRQ(ierr);
  ierr = PetscObjectChangeTypeName((PetscObject)V,VECSEQGPU);CHKERRQ(ierr);

  V->ops->dot             = VecDot_SeqGPU;
  V->ops->norm            = VecNorm_SeqGPU;
  V->ops->tdot            = VecTDot_SeqGPU;
  V->ops->scale           = VecScale_SeqGPU;
  V->ops->copy            = VecCopy_SeqGPU;
  V->ops->set             = VecSet_SeqGPU;
  V->ops->setvalues       = VecSetValues_SeqGPU;
  V->ops->swap            = VecSwap_SeqGPU;
  V->ops->axpy            = VecAXPY_SeqGPU;
  V->ops->axpby           = VecAXPBY_SeqGPU;
  V->ops->axpbypcz        = VecAXPBYPCZ_SeqGPU;
  V->ops->pointwisemult   = VecPointwiseMult_SeqGPU;
  V->ops->pointwisedivide = VecPointwiseDivide_SeqGPU;
  V->ops->maxpointwisedivide = VecMaxPointwiseDivide_SeqGPU;
  V->ops->setrandom       = VecSetRandom_SeqGPU;
  V->ops->dot_local       = VecDot_SeqGPU;
  V->ops->tdot_local      = VecTDot_SeqGPU;
  V->ops->norm_local      = VecNorm_SeqGPU;
  V->ops->maxpy           = VecMAXPY_SeqGPU;
  V->ops->mdot            = VecMDot_SeqGPU;
  V->ops->aypx            = VecAYPX_SeqGPU; 
  V->ops->waxpy           = VecWAXPY_SeqGPU;
  V->ops->dotnorm2        = VecDotNorm2_SeqGPU;
  V->ops->placearray      = VecPlaceArray_SeqGPU;
  V->ops->replacearray    = VecReplaceArray_SeqGPU;
  V->ops->resetarray      = VecResetArray_SeqGPU;
  V->ops->destroy         = VecDestroy_SeqGPU;
  V->ops->destroyvecs     = VecDestroyVecs_SeqGPU;
  V->ops->duplicate       = VecDuplicate_SeqGPU;
  V->ops->duplicatevecs   = VecDuplicateVecs_SeqGPU;
  V->ops->getarray        = VecGetArray_SeqGPU;
  V->ops->restorearray    = VecRestoreArray_SeqGPU;
  V->ops->getlocalsize    = VecGetLocalSize_SeqGPU;
  V->ops->getsize         = VecGetSize_SeqGPU;
  V->ops->view            = VecView_SeqGPU;
  V->ops->max             = VecMax_SeqGPU;
  V->ops->min             = VecMin_SeqGPU;
  V->petscnative=PETSC_FALSE;
  seqgpu->syncState      = VEC_UNALLOC;
  seqgpu->unplacedarray=PETSC_NULL;
  seqgpu->array_allocated=PETSC_NULL;
  seqgpu->array=PETSC_NULL;
  /* create an associated stream */
  cms[0] = cudaStreamCreate(&(seqgpu->streamid));
  /* allocate the variable for vector size */
  cms[1]=cudaMalloc((void**)&(seqgpu->length),sizeof(int));
  /* send vec length size to device */
  ccs[0]=cudaMemcpyAsync((void*)seqgpu->length,
               (void*)&(V->map->n),sizeof(int),cudaMemcpyHostToDevice,seqgpu->streamid);
  /* allocate the vector on device */
  cms[2]=cudaMalloc((void**)&(seqgpu->devptr),V->map->n*sizeof(double));
  ccs[1]=cudaMemsetAsync((void*)seqgpu->devptr,0,V->map->n*sizeof(double),seqgpu->streamid);
  /* allocate the variable for vector offsets */
  cms[3]=cudaMalloc((void**)&(seqgpu->offset),sizeof(int));
  /* allocate the variable for vector segment length */
  cms[4]=cudaMalloc((void**)&(seqgpu->segment),sizeof(int));
  /* allocate the variable for vector single value result */
  cms[5]=cudaMalloc((void**)&(seqgpu->zval),sizeof(double));
  cms[6]=cudaMalloc((void**)&(seqgpu->scalar),sizeof(double));
  /* using pinned memory (could be a resource hog with very large arrays) */
  ierr = PinnedMalloc(&(seqgpu->cpuptr),V->map->n*sizeof(double));CHKERRQ(ierr);

  /* ierr = PetscMalloc(V->map->n*sizeof(PetscScalar),&(seqgpu->cpuptr)); */
  ierr = PetscMemzero(seqgpu->cpuptr,V->map->n*sizeof(double));CHKERRQ(ierr);
  seqgpu->syncState=VEC_ALLOC;
  /* printf("VmapN: %d\n",V->map->n);*/
  #if(DEBUGVEC)
    #if(VVERBOSE)
       printf("Call to VecCreate_SeqGPU\n");
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"on cudaStreamCreate VecCreate_SeqGPU");  CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[1],"Alloc devlength in VecCreate_SeqGPU");   CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"Copy H2D devlength in VecCreate_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[2],"Alloc of devptr in VecCreate_SeqGPU");   CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[1],"on device cudaMemSet VecCreate_SeqGPU"); CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[3],"Alloc devoffset in VecCreate_SeqGPU");   CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[4],"Alloc dev segment in VecCreate_SeqGPU"); CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[5],"Alloc dev zval in VecCreate_SeqGPU");    CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[6],"Alloc dev scalar in VecCreate_SeqGPU");    CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VecDestroy_SeqGPU"
PetscErrorCode VecDestroy_SeqGPU(Vec v){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
#if(DEBUGVEC && VVERBOSE)
     printf("Call to VecDestroyArray_SeqGPU\n"); 
  #endif
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  if(vd && vd->syncState != VEC_UNALLOC){
      cms[0]=cudaFree(vd->devptr);  vd->devptr=PETSC_NULL;
      cms[1]=cudaFree(vd->length);  vd->length=PETSC_NULL;
      cms[2]=cudaFree(vd->segment); vd->segment=PETSC_NULL;
      cms[3]=cudaFree(vd->zval);    vd->zval=PETSC_NULL;
      cms[4]=cudaFree(vd->scalar);  vd->scalar=PETSC_NULL;
      cms[5] = cudaStreamDestroy(vd->streamid);
      ierr = PinnedFree(vd->cpuptr); CHKERRQ(ierr);
      /* ierr = PetscFree(vd->cpuptr);CHKERRQ(ierr); */
      #if(DEBUGVEC)
        ierr=VecCheckCUDAStatus(cms[0],"destroying devptr in VecDestroy_SeqGPU"); CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[1],"destroying length in VecDestroy_SeqGPU"); CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[2],"destroying segment in VecDestroy_SeqGPU");CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[3],"destroying zval in VecDestroy_SeqGPU");   CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[4],"destroying scalar in VecDestroy_SeqGPU"); CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[5],"destroying stream in VecDestroy_SeqGPU"); CHKERRQ(ierr);
      #endif
      vd->syncState = VEC_UNALLOC;
  }
  ierr = PetscObjectDepublish(v);CHKERRQ(ierr);
#if defined(PETSC_USE_LOG)
  PetscLogObjectState((PetscObject)v,"Length=%D",v->map->n);
#endif
  ierr = PetscFree(v->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecDestroyVecs_SeqGPU"
PetscErrorCode  VecDestroyVecs_SeqGPU(PetscInt m,Vec *vv){
  PetscFunctionBegin;
  #if(DEBUGVEC && VVERBOSE)
     printf("Call to VecDestroyVecs_SeqGPU\n");
  #endif
  PetscErrorCode ierr;
  PetscInt i;
   /* destroy the internal part */
  for(i=0;i<m;i++){
    ierr = VecDestroy(&vv[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(vv); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecView_SeqGPU_ASCII"
PetscErrorCode VecView_SeqGPU_ASCII(Vec xin,PetscViewer viewer){
  PetscFunctionBegin;
  printf("VecView_Seq_ASCII() (***EMPTY***)\n");
  PetscFunctionReturn(0);
}

EXTERN_C_END
