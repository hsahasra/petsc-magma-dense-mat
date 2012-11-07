#include <petscconf.h>
#include <petscsys.h>
#include <petscerror.h>
PETSC_CUDA_EXTERN_C_BEGIN
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include <float.h>
#include <petsc-private/vecimpl.h>          /*I "petscvec.h" I*/
#include <../src/vec/vec/impls/dvecimpl.h>
#include <../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h>
PETSC_CUDA_EXTERN_C_END



EXTERN_C_BEGIN
__constant__ int     integerSymbol;
__constant__ int2    integer2Symbol;
__constant__ int3    integer3Symbol;
__constant__ int     devN;/* vector length */
__constant__ double  dblScalarValue;/* utility var */
__constant__ double2 dblScalar2Value;/* utility var */
__constant__ float   fltScalarValue;/* utility var */
__constant__ float2  fltScalar2Value;/* utility var */

static cudaError_t ccs[16];
static cudaError_t cms[16];

/* Valid pointer check function (probably doesn't work) */
PetscBool valid(void *p){
  extern char _etext;
  if((p != PETSC_NULL) && ((char*) p > &_etext)){
    return PETSC_TRUE;
  }else{
    return PETSC_FALSE;
  }
}

/* ---------------------------------------------------------
// helper function for error checking
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


#undef __FUNCT__
#define __FUNCT__ "VecCheckCUDAStatus"
PetscErrorCode VecCheckCUDAStatus(cudaError_t cs,const char *msg){
  PetscFunctionBegin;
    if(cs!=cudaSuccess){
      fprintf(stderr, "Cuda error: %s: %s.\n",msg,cudaGetErrorString(cs));
      fflush(NULL);
      PetscFunctionReturn(PETSC_ERR_LIB);
    }
  PetscFunctionReturn(0);
}

/* -------------------- end error checkers ------------------- */




/* ****************************************************************************
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
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecSetRandom_SeqGPU\n");
  #endif
  if(xd->syncState==VEC_ALLOC || xd->syncState==VEC_CPU){
    for(i=0; i<x->map->n; i++){
       ierr = PetscRandomGetValue(r,&xd->cpuptr[i]);CHKERRQ(ierr);
    }
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }else if(xd->syncState==VEC_SYNCHED || xd->syncState==VEC_GPU){
    dimGrid.x=ceil((float)x->map->n/(float)TCOUNT);
    dimBlock.x=TCOUNT;
    while(dimGrid.x>MAXBLOCKS){
      dimGrid.x/=2;
      dimBlock.x*=2;
    }
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

/*------------------------end random generator ------------------------*/



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
  while(dimGrid.x>MAXBLOCKS){
      dimGrid.x/=2;
      dimBlock.x*=2;
  }
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
      #if(DEBUGVEC && VERBOSE)
      printf("In kernCompare found an element mismatch: %e\n",value);
      #endif
      blockflag=1;
    }
    if(*lx!=*ly){
      #if(DEBUGVEC && VERBOSE)
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
  for(i=0;i<x->map->n;i++){
    if(xd->cpuptr[i]!=0)printf("cpu[%d]: %e\n",i,xd->cpuptr[i]);
  }
  /* ierr= PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);*/
  /* ierr =VecView_Seq_ASCII(x,viewer);CHKERRQ(ierr); */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecGetLocalSize_SeqGPU"
PetscErrorCode VecGetLocalSize_SeqGPU(Vec x, PetscInt *localsize){
  PetscFunctionBegin;
  #if(DEBUGVEC && VERBOSE)
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
  #if(DEBUGVEC && VERBOSE)
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
    #if(DEBUGVEC && VERBOSE)
    printf("kernCheck: x[%d]: %e, length: %d\n",tid,x[tid],*n);
    #endif
  }
}

/*------------------------------ end info -------------------------------*/


/*---------------------------- copy functions ---------------------------*/
#undef __FUNCT__
#define __FUNCT__ "VecCopyBlockDevice"
PetscErrorCode VecCopyBlockDevice(Vec d, Vec s, PetscInt doffset, PetscInt soffset, PetscInt blocksize){
  PetscFunctionBegin;
  printf("Call to VecCopyBlockDevice (**** EMPTY ****)\n");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecCopyOverDevice"
PetscErrorCode VecCopyOverDevice(Vec d,Vec s){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU* dd = (Vec_SeqGPU*)d->data;
  Vec_SeqGPU* sd = (Vec_SeqGPU*)s->data;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecCopyOverDevice\n");
  #endif
  dim3 dimGrid;
  dim3 dimBlock;

  if(s->map->n!=d->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  if(sd->syncState==VEC_CPU){/* synch y */
    ierr = VecCopyOverH2D(s,sd->cpuptr);CHKERRQ(ierr);
    sd->syncState=VEC_SYNCHED;
    cudaStreamSynchronize(sd->stream);
  }
  ccs[0]=cudaMemcpyAsync(dd->devptr,sd->devptr,
               s->map->n*sizeof(PetscScalar),cudaMemcpyDeviceToDevice,dd->stream);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(ccs[0],"on copy D2D in VecCopyOverDevice");CHKERRQ(ierr);
    /*PetscBool same;
    ierr = VecCompare_SeqGPU(d, s, &same,0,s->map->n);CHKERRQ(ierr);
    printf("**** compare**** s and d the same?: %d\n",same); */
  #endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "kernCopyLen"
__global__ void kernCopyLen(int* ly,int* lx){
  if(threadIdx.x==0)*ly=*lx;
}

#undef __FUNCT__
#define __FUNCT__ "VecCopyBlockH2D"
PetscErrorCode VecCopyBlockH2D(Vec v,PetscScalar *y, PetscInt offset, PetscInt blocksize){
  PetscFunctionBegin;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  ccs[0]=cudaMemcpyAsync(&(vd->devptr[offset]),y,
               blocksize*sizeof(PetscScalar),cudaMemcpyHostToDevice,vd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Call to VecCopyBlockH2D\n");
    #endif
    PetscErrorCode ierr;
    ierr = VecCheckCUDAStatus(ccs[0],"on copy H2D in VecCopyBlockH2D");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCopyOverH2D"
PetscErrorCode VecCopyOverH2D(Vec v,PetscScalar *y){
  PetscFunctionBegin;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  /*
  int i;
  for(i=0;i<v->map->n;i++){
    if(y[i]!=0)printf("y[%d]: %e\n",i,y[i]);
  }
  */
  ccs[0]=cudaMemcpyAsync(vd->devptr,y,
               v->map->n*sizeof(PetscScalar),cudaMemcpyHostToDevice,vd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Call to VecCopyOverH2D\n");
    #endif
    PetscErrorCode ierr;
    ierr = VecCheckCUDAStatus(ccs[0],"on copy H2D in VecCopyOverH2D");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCopyBlockD2H"
PetscErrorCode VecCopyBlockD2H(Vec v,PetscScalar *y,PetscInt offset, PetscInt blocksize){
  PetscFunctionBegin;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  ccs[0]=cudaMemcpyAsync(y,&(vd->devptr[offset]),
               blocksize*sizeof(PetscScalar),cudaMemcpyDeviceToHost,vd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Call to VecCopyBlockD2H\n");
    #endif
    PetscErrorCode ierr;
    ierr = VecCheckCUDAStatus(ccs[0],"on copy D2H in VecCopyBlockD2H");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCopyOverD2H"
PetscErrorCode VecCopyOverD2H(Vec v,PetscScalar *y){
  PetscFunctionBegin;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  ccs[0]=cudaMemcpyAsync(y,vd->devptr,
               v->map->n*sizeof(PetscScalar),cudaMemcpyDeviceToHost,vd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
      printf("Call to VecCopyOverD2H\n");
    #endif
    PetscErrorCode ierr;
    ierr = VecCheckCUDAStatus(ccs[0],"on copy D2H in VecCopyOverD2H");CHKERRQ(ierr); 
  #endif
  PetscFunctionReturn(0);
}

/*---------------------------- end copy functions --------------------------*/

/*------------------------------ set functions -----------------------------*/
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "VecSetValues_SeqGPU"
/*
   VecSetValues - Inserts or adds values into certain locations of a vector.
*/
PetscErrorCode VecSetValues_SeqGPU(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecSetValues_SeqGPU\n");
  #endif
  if(xd->syncState==VEC_CPU || xd->syncState==VEC_SYNCHED){
    if(iora==INSERT_VALUES){
      for(i=0;i<ni;i++){
         xd->cpuptr[i]=y[i];
      }
      ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
      xd->syncState=VEC_SYNCHED;
    }else{
      /* ADD_VALUES not supported now */
      printf("Call to VecSetValues_SeqGPU: ADD_VALUES (*** EMPTY ***)\n");
    }
  }else{
      if(iora==INSERT_VALUES){/* not efficient */
        PetscScalar yval=0;
        for(i=0;i<ni;i++){
          yval=y[i];
          ierr = VecCopyBlockH2D(x,&yval,ix[i],1);CHKERRQ(ierr);
        }
      }
      xd->syncState=VEC_GPU;
  }
  PetscFunctionReturn(0);
}
EXTERN_C_END



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
  while(dimGrid.x>MAXBLOCKS){
     dimGrid.x/=2;
     dimBlock.x*=2;
  }
  Vec_SeqGPU* xd = (Vec_SeqGPU*)xin->data;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecSet_SeqGPU, alpha: %e\n",alpha);
  #endif
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,
            PETSC_ERR_MEM,"*** In VecSet_SeqGPU, Vec not allocated. ***\n");
  }else{
    ccs[0]=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
    #endif
       kernSet<<<dimGrid,dimBlock,dimBlock.x*sizeof(double)>>>(xd->devptr,xd->length);
    #if(DEBUGVEC)
       #if(VERBOSE)
          printf("In VecSet_SeqGPU: blocks: %d, threads: %d\n",dimGrid.x, dimBlock.x);
       #endif
       ierr = VecCheckCUDAError("Call to kernSet. "); CHKERRQ(ierr);
    #endif
    xd->syncState=VEC_GPU;
  }
  PetscFunctionReturn(0);
}

extern __shared__ double sharedSet[];
#undef __FUNCT__
#define __FUNCT__ "kernSet"
__global__ void kernSet(double* x, int* n){
  int tid = threadIdx.x + blockDim.x*blockIdx.x;
  double* setptr = sharedSet;
  setptr[threadIdx.x] = dblScalarValue;
  if(tid<*n) x[tid] = setptr[threadIdx.x];
}




#undef __FUNCT__
#define __FUNCT__ "VecScale_SeqGPU"
PetscErrorCode VecScale_SeqGPU(Vec x, PetscScalar alpha){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  dim3 dimGrid,dimBlock;
  dimGrid.x=ceil((float)x->map->n/((float)TCOUNT));
  dimBlock.x=TCOUNT;
  while(dimGrid.x>MAXBLOCKS){
    dimGrid.x/=2;
    dimBlock.x*=2;
  }
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  #if(DEBUGVEC && VERBOSE)
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
    ccs[0] = cudaMemsetAsync(xd->devptr,0,x->map->n*sizeof(double),xd->stream);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[0],"error in cudaMemset");CHKERRQ(ierr);
    #endif
  }else if (alpha != 1.0){
    ccs[0]=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
    #endif
       kernScale<<<dimGrid,dimBlock,dimBlock.x*sizeof(double),xd->stream>>>(xd->devptr,xd->length);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("Call to kernScale."); CHKERRQ(ierr);
    #endif
  }
  xd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}

extern __shared__ double sharedScale[];
#undef __FUNCT__
#define __FUNCT__ "kernScale"
__global__ void kernScale(double* x, int* n){
  int tid = threadIdx.x + blockDim.x*blockIdx.x;
  double* scaleptr = sharedScale;
  __shared__ double scalar;
  scalar=dblScalarValue;
  if(tid<*n){
    scaleptr[threadIdx.x] = x[tid];
    scaleptr[threadIdx.x]*= scalar;
    x[tid] = scaleptr[threadIdx.x];
  }
}

/*---------------------------- end set and scale ---------------------------*/


/*-------------------------- dot product functions -------------------------*/

#undef __FUNCT__
#define __FUNCT__ "VecTDot_SeqGPU"
PetscErrorCode VecTDot_SeqGPU(Vec xin,Vec yin,PetscScalar *z){
  PetscFunctionBegin;
  printf("VecTDot_SeqGPU (***EMPTY***)\n");
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecDot_SeqGPU"
PetscErrorCode VecDot_SeqGPU(Vec x,Vec y,PetscScalar *z){
  PetscFunctionBegin;
  if(x->map->n!=y->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  PetscErrorCode ierr;
  double *devScratch,*devPartial,*hostPartial;
  PetscInt i,chunks=0,secondPhase,segment,partialsize,scratchsize;
  cudaStream_t* dotstream;
  dim3 dimGrid, dimBlock;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  if(xd->syncState==VEC_CPU){
    #if(DEBUGVEC && VERBOSE)
       printf("xd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  if(yd->syncState==VEC_CPU){
    #if(DEBUGVEC && VERBOSE)
       printf("yd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) /(float)(CHUNKWIDTH));
  dotstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));
  /* make sure the segment size for each chunk is correct */
  if(chunks>1) segment = (int) CHUNKWIDTH;
  else segment = x->map->n;
  dimGrid.x=ceil(((float)segment)/(float)THRDOTCNT);
  dimBlock.x = THRDOTCNT;
  /* allocate gridwide scratch array */
  scratchsize=chunks*dimGrid.x*sizeof(double);
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize);/* scratch pad */
  ccs[0] = cudaMemsetAsync(devScratch,0,scratchsize,xd->stream);
  ccs[1]=cudaMemcpyAsync(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice,yd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
      printf("Call to VecDot_SeqGPU\n");
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecDot_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecDot_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[1],"on copy segment size H2D in VecDot_SeqGPU");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();/* make sure everyone is ready */
  for(i=0;i<chunks;i++){  /* streaming async kernel calls */
    cudaStreamCreate(&(dotstream[i]));
    cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,dotstream[i]);
    /* Overlapping execution */
    kernDot<<<dimGrid,dimBlock,0,dotstream[i]>>>(xd->devptr,yd->devptr,
                                                          xd->segment,
                                                          xd->length,
                                                          xd->offset,
                                                          (devScratch+i*dimGrid.x));
  }
  secondPhase = scratchsize/sizeof(double);
  if(secondPhase>1){/* begin next reduction */
    dimGrid.x = ceil((float)secondPhase/(float)THRDOTCNT2);
    dimBlock.x  = THRDOTCNT2;
    /* allocate last reduction array */
    partialsize=dimGrid.x*sizeof(double);
    cms[1]=cudaMalloc((void**)&devPartial,partialsize);/* partial results to be combined */
    ccs[2]=cudaMemsetAsync(devPartial,0,partialsize,yd->stream);
    ccs[3] = cudaMemcpyAsync(xd->segment,&secondPhase,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
    #if(DEBUGVEC)
       #if(VERBOSE)
         printf("DOT phase2: blks: %d, partial: %d\n",dimGrid.x,partialsize/sizeof(double));
       #endif
       ierr = VecCheckCUDAStatus(cms[1],"devPartial alloc in VecDot_SeqGPU");CHKERRQ(ierr);
       ierr = VecCheckCUDAStatus(ccs[2],"devPartial memset in VecDot_SeqGPU");CHKERRQ(ierr);
       ierr = VecCheckCUDAStatus(ccs[3],"on cudaMemcpy(xd->segment)");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    kernRedDot<<<dimGrid,dimBlock,dimBlock.x*sizeof(double)>>>(xd->segment,devScratch,devPartial);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernRedDot launch in VecDot_SeqGPU");CHKERRQ(ierr); 
    #endif
    /* setup copy back array while waiting for kernel to finish */
    ierr = PetscMalloc(partialsize,&hostPartial);CHKERRQ(ierr);
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    ccs[4]=cudaMemcpy(hostPartial,devPartial,partialsize,cudaMemcpyDeviceToHost);/* copy back */
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[4],"on cudaMemcpy(devPartial)");CHKERRQ(ierr);
    #endif
    if(dimGrid.x>1){/* final reduction */
      *z=0.;
      for(i=0;i<dimGrid.x;i++)*z+=hostPartial[i];
    }else{
      *z=hostPartial[0];
    }
    ierr = PetscFree(hostPartial); CHKERRQ(ierr);
    cms[2] = cudaFree(devPartial);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(cms[2],"on cudaFree(devPartial)");CHKERRQ(ierr);
    #endif
  }else{
    ccs[4]=cudaMemcpy(z,devScratch,dimGrid.x*sizeof(double),cudaMemcpyDeviceToHost);/* copy back */
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[4],"on cudaMemcpy(devScratch)");CHKERRQ(ierr);
    #endif
  }

  /* clean up resources */
  for(i=0;i<chunks;i++){
     cudaStreamDestroy(dotstream[i]);
  }
  free(dotstream);
  cms[3] = cudaFree(devScratch);
  #if(DEBUGVEC)
     #if(VERBOSE)
       printf("Zdot: %e\n",*z);
    #endif
    ierr = VecCheckCUDAStatus(cms[3],"on cudaFree(devScratch)");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}

extern __shared__ double sharedRedDot[];
#undef __FUNCT__
#define __FUNCT__ "kernRedDot"
__global__ void kernRedDot(int* size,double* scratch, double* z){/* reduction kernel */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  double* zDot = sharedRedDot;
  zDot[threadIdx.x]=(tid<*size)?scratch[tid]:0.;
  while(i>0){
    if(threadIdx.x<i){
      zDot[threadIdx.x]+=zDot[threadIdx.x+i];
    }
    __syncthreads();
    i/=2;
  }
  if(threadIdx.x==0){
    z[blockIdx.x]=zDot[0];
    //printf("ZDOT block[%d]: %e\n",blockIdx.x,z[blockIdx.x]);
  }
}

#undef __FUNCT__
#define __FUNCT__ "kernDot"
__global__ void kernDot(double* devX, double* devY,
                        int* segmentsize, int* arrsize,
                        int* offset, double* scratch){
  __shared__ double chunkX[THRDOTCNT];
  __shared__ double chunkY[THRDOTCNT];
  __shared__ int n;    n   = *arrsize;

  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  int item = *segmentsize**offset+tid;

  if(item<n){/* read in values to shared mem */
    chunkX[threadIdx.x]=devX[item]; /* offset values */
    chunkY[threadIdx.x]=devY[item]; /* offset values */
  }else{
    chunkX[threadIdx.x]=0.;
    chunkY[threadIdx.x]=0.;
  }
  chunkX[threadIdx.x]*=chunkY[threadIdx.x];
  __syncthreads();
  while(i>0){/* block level reduction */
     if(threadIdx.x<i){
       chunkX[threadIdx.x]+=chunkX[threadIdx.x+i];
     }
     __syncthreads();
     i/=2;
  }/* end while */

  if(threadIdx.x==0)  scratch[blockIdx.x]=chunkX[0];
}


#undef __FUNCT__
#define __FUNCT__ "VecMDot_SeqGPU"
PetscErrorCode  VecMDot_SeqGPU(Vec x,PetscInt nv,const Vec y[],PetscScalar val[]){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt       i;
  for (i=0; i<nv; i++) {
    ierr = VecDot_SeqGPU(x,y[i],&val[i]);CHKERRQ(ierr);
    if(PetscIsInfOrNanScalar(val[i])){
      SETERRQ1(((PetscObject)x)->comm,PETSC_ERR_FP,"Infinite or not-a-number generated in mdot, entry %D",i);
    }
  }
  PetscFunctionReturn(0);
}

/*----------------------------- end dot ----------------------------- */





#undef __FUNCT__
#define __FUNCT__ "VecAXPBY_SeqGPU"
PetscErrorCode VecAXPBY_SeqGPU(Vec yin,PetscScalar beta,PetscScalar alpha,Vec xin){
  /* Y = b*Y + a*X */
  PetscFunctionBegin;
  printf("Call to VecAXPBY_SeqGPU (***EMPTY***)\n");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWAXPY_SeqGPU"
PetscErrorCode VecWAXPY_SeqGPU(Vec w,PetscScalar alpha,Vec x,Vec y){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)w->data;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid, dimBlock;
  #if(DEBUGVEC && VERBOSE)
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
  dimGrid.x=ceil((float)y->map->n/(float)AXPYTCOUNT);
  dimBlock.x=AXPYTCOUNT;
  while(dimGrid.x>MAXBLOCKS){
    dimGrid.x/=2;
    dimBlock.x*=2;
  }
  cudaDeviceSynchronize();
  if(alpha==0.0){
    ierr = VecCopyOverDevice(w,y);CHKERRQ(ierr);
  }else if(alpha==1.0){
    kernWXPY<<<dimGrid,dimBlock,3*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWXPY");CHKERRQ(ierr); 
    #endif
  }else if(alpha==-1.0){
    kernWXMY<<<dimGrid,dimBlock,3*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWXMY");CHKERRQ(ierr);
    #endif
  }else{
    ccs[0]=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
    #endif
    kernWAXPY<<<dimGrid,dimBlock,3*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWAXPY");CHKERRQ(ierr); 
    #endif
  }
  wd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}


extern __shared__ double sharedWAXPY[];
#undef __FUNCT__
#define __FUNCT__ "kernWAXPY"
__global__ void  kernWAXPY(double* devY,double* devX, int* vlen, double* devW){
 /* w <- y + alpha*x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double alphaShared;
  double* chunkX = sharedWAXPY;
  double* chunkY = sharedWAXPY + blockDim.x;
  double* chunkW = sharedWAXPY + 2*blockDim.x;
  alphaShared = dblScalarValue;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]+(chunkX[threadIdx.x]*alphaShared);
    devW[tid]=chunkW[threadIdx.x];
  }
}

extern __shared__ double sharedWXPY[];
#undef __FUNCT__
#define __FUNCT__ "kernWXPY"
__global__ void  kernWXPY(double* devY,double* devX, int* vlen, double* devW){
 /* w <- y + x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  double* chunkX = sharedWXPY;
  double* chunkY = sharedWXPY + blockDim.x;
  double* chunkW = sharedWXPY + 2*blockDim.x;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]+chunkX[threadIdx.x];
    devW[tid]=chunkW[threadIdx.x];
  }
}

extern __shared__ double sharedWXMY[];
#undef __FUNCT__
#define __FUNCT__ "kernWXMY"
__global__ void  kernWXMY(double* devY,double* devX, int* vlen, double* devW){
 /* w <- y + alpha*x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  double* chunkX = sharedWXMY;
  double* chunkY = sharedWXMY + blockDim.x;
  double* chunkW = sharedWXMY + 2*blockDim.x;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]-chunkX[threadIdx.x];
    devW[tid]=chunkW[threadIdx.x];
  }
}

#undef __FUNCT__
#define __FUNCT__ "VecMAXPY_SeqGPU"
PetscErrorCode VecMAXPY_SeqGPU(Vec x,PetscInt nv,const PetscScalar* alpha,Vec *y){
  /* y = y + sum(a[i]*x[i]) */
  PetscFunctionBegin;
  if(DEBUGVEC && VERBOSE)printf("VecMAXPY_SeqGPU: alpha: %e\n",*alpha);
  PetscErrorCode ierr;
  PetscInt i;
  PetscScalar *devW;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=PETSC_NULL;
  cms[0] = cudaMalloc((void**)&devW,x->map->n*sizeof(double));
  ccs[0] = cudaMemsetAsync(devW,0,x->map->n*sizeof(double),xd->stream);
  dim3 dimGrid;  dim3 dimBlock;
  dimGrid.x=ceil((float)x->map->n/(float)AXPYTCOUNT);
  dimBlock.x=AXPYTCOUNT;
  while(dimGrid.x>MAXBLOCKS){
    dimGrid.x/=2;
    dimBlock.x*=2;
  }
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Number of vectors in MAXPY: %d, blocks: %d, threads: %d\n",nv,dimGrid.x,dimBlock.x);
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"error in device malloc VecMAXPY_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"error in device memset VecMAXPY_SeqGPU");CHKERRQ(ierr);
  #endif

  for(i=0;i<nv;i++){
     if(y[i]->map->n!=x->map->n){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
    }
    yd=(Vec_SeqGPU*)y[i]->data;
    if(yd->syncState==VEC_CPU){/* synch x */
      ierr = VecCopyOverH2D(y[i],yd->cpuptr);CHKERRQ(ierr);
      yd->syncState=VEC_SYNCHED;
    }
    ccs[1]=cudaMemcpy(yd->scalar,&alpha[i],sizeof(double),cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[1],"error in symbol copy to device");CHKERRQ(ierr); 
    #endif
    cudaDeviceSynchronize();
    if(alpha[i]==0){
      continue;
    }else if(alpha[i]==1.){
      kernXPY<<<dimGrid,dimBlock,2*dimBlock.x*sizeof(double)>>>(devW,yd->devptr,yd->length);
    }else{
      kernAXPY<<<dimGrid,dimBlock,2*dimBlock.x*sizeof(double)>>>(devW,yd->devptr,yd->length,yd->scalar);
    }
    #if(DEBUGVEC)
         ierr = VecCheckCUDAError("kernel call to kernAXPY or kernXPY in VecMAXPY_SeqGPU");CHKERRQ(ierr);
    #endif
  }/* end for */
  if(xd->syncState==VEC_CPU){/* synch x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  cudaDeviceSynchronize();
  kernXPY<<<dimGrid,dimBlock,2*dimBlock.x*sizeof(double)>>>(xd->devptr,devW,xd->length);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAError("kernel call to kernXPY");CHKERRQ(ierr);
  #endif
  cms[1] = cudaFree(devW);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAStatus(cms[1],"on cudaFree");CHKERRQ(ierr);
  #endif
  xd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}

extern __shared__ double sharedXPY[];
#undef __FUNCT__
#define __FUNCT__ "kernXPY"
__global__ void  kernXPY(double* devY,double* devX, int* vlen){
 /* y <- y + x */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  double* chunkX = sharedXPY;
  double* chunkY = sharedXPY + blockDim.x;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkY[threadIdx.x]+=chunkX[threadIdx.x];
    devY[tid]=chunkY[threadIdx.x];
  }
}

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
  #if(DEBUGVEC && VERBOSE)
      printf("VecAXPY_SeqGPU\n");
  #endif
  dim3 dimGrid, dimBlock;
  dimGrid.x=ceil((float)x->map->n/(float)AXPYTCOUNT);
  dimBlock.x=AXPYTCOUNT;
  while(dimGrid.x>MAXBLOCKS){
    dimGrid.x/=2;
    dimBlock.x*=2;
  }
  cudaDeviceSynchronize();
  if(alpha==1.){
    kernXPY<<<dimGrid,dimBlock,2*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,yd->length);
  }else if(alpha!=0.){
    ccs[0]=cudaMemcpy(yd->scalar,&alpha,sizeof(double),cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
    #endif
    kernAXPY<<<dimGrid,dimBlock,2*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,yd->length,yd->scalar);
  }
  #if(DEBUGVEC)
   ierr = VecCheckCUDAError("kernel call in VecAXPY_SeqGPU");CHKERRQ(ierr);
  #endif
  yd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}


extern __shared__ double sharedAXPY[];
#undef __FUNCT__
#define __FUNCT__ "kernAXPY"
__global__ void  kernAXPY(double* devY,double* devX, int* vlen,double *scalar){
 /* y <- y + alpha*x */
  __shared__ double alphaShared;
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  double* chunkX = sharedAXPY;
  double* chunkY = sharedAXPY + blockDim.x;
  alphaShared = *scalar;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkY[threadIdx.x]+=chunkX[threadIdx.x]*alphaShared;
    devY[tid]=chunkY[threadIdx.x];
  }
}

#undef __FUNCT__
#define __FUNCT__ "VecAXPBYPCZ_SeqGPU"
PetscErrorCode VecAXPBYPCZ_SeqGPU(Vec x,PetscScalar alpha,PetscScalar beta,PetscScalar gamma,Vec y,Vec z){
  PetscFunctionBegin;
  #if(DEBUGVEC)
     PetscErrorCode ierr;
     #if(VERBOSE)
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
  while(dimGrid.x>MAXBLOCKS){
    dimGrid.x/=2;
    dimBlock.x*=2;
  }
  ccs[0]=cudaMemcpyToSymbol("dblScalar2Value",(void*)&alphabeta,sizeof(double2),0,cudaMemcpyHostToDevice);
  ccs[1]=cudaMemcpyToSymbol("dblScalarValue",(void*)&gamma,sizeof(double),0,cudaMemcpyHostToDevice);
  #if(DEBUGVEC)
   ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
   ierr = VecCheckCUDAStatus(ccs[1],"error in symbol copy to device");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();
  kernAXPBYPCZ<<<dimGrid,dimBlock,4*dimBlock.x*sizeof(double)>>>(devX->devptr,devY->devptr,devZ->devptr,devX->length);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAError("launch kernAXPBYPCZ");CHKERRQ(ierr); 
  #endif
  PetscFunctionReturn(0);
}

extern __shared__ double sharedAXPBYPCZ[];
#undef __FUNCT__
#define __FUNCT__ "kernAXPBYPCZ"
__global__ void kernAXPBYPCZ(double* devX, double* devY, double* devZ, int* len){
  /* x <- alpha*x + beta*y + gamma*z */
  __shared__ int localn;
  localn = *len;
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  double* work = sharedAXPBYPCZ;
  double* chunkX = sharedAXPBYPCZ + blockDim.x;
  double* chunkY = sharedAXPBYPCZ + 2*blockDim.x;
  double* chunkZ = sharedAXPBYPCZ + 3*blockDim.x;
  if(tid<localn){
    /* read in values to shared */
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkZ[threadIdx.x]=devZ[tid];

    /* do flops */
    if(dblScalarValue){
      work[threadIdx.x]=dblScalarValue*chunkZ[threadIdx.x];
    }else{
      work[threadIdx.x]=0.;
    }

    if(dblScalar2Value.y){
      work[threadIdx.x]+=dblScalar2Value.y*chunkY[threadIdx.x];
    }
    if(dblScalar2Value.x){
      work[threadIdx.x]+=dblScalar2Value.x*chunkX[threadIdx.x];
    }

    /* write back */
    devX[tid]=work[threadIdx.x];
  }
  return;
}

/*---------------------------- end level 2 ------------------------------ */

/*------------------------- pointwise functions ------------------------- */
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
  Vec_SeqGPU *wd=(Vec_SeqGPU*)y->data;
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
  while(dimGrid.x>MAXBLOCKS){
    dimGrid.x/=2;
    dimBlock.x*=2;
  }
  cudaDeviceSynchronize();
  kernPMULT<<<dimGrid,dimBlock,3*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAError("kernel call to kernPMULT");CHKERRQ(ierr);
  #endif

  PetscFunctionReturn(0);
}

extern __shared__ double sharedPMULT[];
#undef __FUNCT__
#define __FUNCT__ "kernPMULT"
__global__ void  kernPMULT(double* devY,double* devX, int* vlen, double* devW){
 /* w <- x./y */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  double* chunkX = sharedPMULT;
  double* chunkY = sharedPMULT + blockDim.x;
  double* chunkW = sharedPMULT + 2*blockDim.x;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkX[threadIdx.x]*chunkY[threadIdx.x];
    devW[tid]=chunkW[threadIdx.x];
  }
}


#undef __FUNCT__
#define __FUNCT__ "VecMaxPointwiseDivide_SeqGPU"
PetscErrorCode VecMaxPointwiseDivide_SeqGPU(Vec x,Vec y,PetscReal *max){
  PetscFunctionBegin;
  #if(DEBUGVEC && VERBOSE)
     printf("VecMaxPointwiseDivide_SeqGPU...");
  #endif
  if(x->map->n!=y->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }
  PetscErrorCode ierr;
  PetscScalar *devScratch,*devPartial,*hostPartial;
  PetscInt i,chunks=0,segment,partialsize,scratchsize,secondPhase;
  cudaStream_t* pwdstream;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid;  dim3 dimBlock;
  if(yd->syncState==VEC_CPU){/* synch up y */
    #if(DEBUGVEC && VERBOSE)
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
  chunks = ceil( ((float)x->map->n) /(float)(CHUNKWIDTH));
  pwdstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));
  /* make sure the segment size for each chunk is correct */
  if(chunks>1)segment = (int) (CHUNKWIDTH);
  else segment = x->map->n;
  dimGrid.x=ceil(((float)segment)/(float)PDIVTCOUNT);
  dimBlock.x  = PDIVTCOUNT;

  scratchsize = chunks*dimGrid.x*sizeof(double);
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize);
  ccs[0] = cudaMemsetAsync(devScratch,0,scratchsize,yd->stream);
  ccs[1]=cudaMemcpyAsync(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[1],"on copy segment size H2D in VecMPWD_SeqGPU");CHKERRQ(ierr);
  #endif

  cudaDeviceSynchronize();
  for(i=0;i<chunks;i++){
    cudaStreamCreate(&(pwdstream[i]));
    cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,pwdstream[i]);
    /* Overlapping execution */
    kernMAXPDIV<<<dimGrid,dimBlock,0,pwdstream[i]>>>(xd->devptr,yd->devptr,
                                                     xd->segment,
                                                     xd->length,
                                                     xd->offset,
                                                     (devScratch+i*dimGrid.x));
  }/* end for-loop */

  secondPhase = scratchsize/sizeof(double);
  if(secondPhase>1){/* begin next reduction */
    dimGrid.x = ceil((float)secondPhase/(float)PDIVTCOUNT2);
    dimBlock.x  = PDIVTCOUNT2;
    /* allocate last reduction array */
    partialsize = dimGrid.x*sizeof(double);
    cms[1] = cudaMalloc((void**)&devPartial,partialsize);
    ccs[2] = cudaMemsetAsync(devPartial,0,partialsize,yd->stream);
    ccs[3] = cudaMemcpyAsync(xd->segment,&secondPhase,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
    #if(DEBUGVEC)
       #if(VERBOSE)
         printf("MAXDIV phase2: blks: %d, partial: %d\n",dimGrid.x,partialsize/sizeof(double));
       #endif
       ierr = VecCheckCUDAStatus(cms[1],"devPartial alloc in VecMPWD_SeqGPU"); CHKERRQ(ierr);
       ierr = VecCheckCUDAStatus(ccs[2],"devPartial memset in VecMPWD_SeqGPU");CHKERRQ(ierr);
       ierr = VecCheckCUDAStatus(ccs[3],"on copy chunks H2D in VecMPWD_SeqGPU");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    kernMAX<<<dimGrid,dimBlock,dimBlock.x*sizeof(double)>>>(xd->segment,devScratch,devPartial);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernRedNorm_double launch in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
    ierr = PetscMalloc(partialsize,&hostPartial);CHKERRQ(ierr);
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    ccs[4]=cudaMemcpy(hostPartial,devPartial,partialsize,cudaMemcpyDeviceToHost);/* copy back */
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[4],"on devPartial copy D2H");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    /* final reduction */
    if(dimGrid.x>1){
      *max=0.;
      for(i=0;i<dimGrid.x;i++){
        *max=PetscMax(hostPartial[i],*max);
      }
    }else{
      *max=hostPartial[0];
    }
    ierr = PetscFree(hostPartial); CHKERRQ(ierr);
    cms[2] = cudaFree(devPartial);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(cms[2],"on cudaFree(devPartial)");CHKERRQ(ierr);
    #endif

  }else{
    ccs[4]=cudaMemcpy(max,devScratch,sizeof(double),cudaMemcpyDeviceToHost);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[4],"on max copy D2H");CHKERRQ(ierr);
    #endif
  }
  cms[3] = cudaFree(devScratch);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAStatus(cms[3],"on cudaFree(devScratch)");CHKERRQ(ierr);
  #endif
  for(i=0;i<chunks;i++) cudaStreamDestroy(pwdstream[i]);

  #if(DEBUGVEC && VERBOSE)
     printf("max: %e\n",*max);
  #endif
  PetscFunctionReturn(0);
}


extern __shared__ double sharedMAX[];
#undef __FUNCT__
#define __FUNCT__ "kernMAX"
__global__ void  kernMAX(int* size, double* maxlist,double* max){
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  double* mlist = sharedMAX;
  mlist[threadIdx.x]=(tid<*size)?maxlist[tid]:0.;
  /* printf("mlist[%d]: %e\n",threadIdx.x,mlist[threadIdx.x]); */
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i){
      mlist[threadIdx.x]=fmax(mlist[threadIdx.x],mlist[threadIdx.x+i]);
    }
    __syncthreads();
    i/=2;
  }
  if(threadIdx.x==0){
    max[blockIdx.x] = mlist[0];
    //printf("MAXDIV block[%d]: %e\n",blockIdx.x,max[blockIdx.x]);
  }
}


#undef __FUNCT__
#define __FUNCT__ "kernMAXPDIV"
__global__ void  kernMAXPDIV(double* devX,double* devY, int* segmentsize,
                             int* arrsize,int* offset,double* scratch){
 /* w <- max(abs(x./y)) */
  __shared__ double chunkY[PDIVTCOUNT];
  __shared__ double chunkX[PDIVTCOUNT];
  __shared__ double chunkW[PDIVTCOUNT];
  __shared__ int n;    n   = *arrsize;
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  int item = *segmentsize**offset+tid;
  if(item<n){
    chunkX[threadIdx.x]=devX[item];
    chunkY[threadIdx.x]=devY[item];
  }else{
    chunkX[threadIdx.x]=0.;
    chunkY[threadIdx.x]=0.;
  }
  if(chunkY[threadIdx.x]!=0.)chunkW[threadIdx.x]=fabs(chunkX[threadIdx.x]/chunkY[threadIdx.x]);
  else chunkW[threadIdx.x]=fabs(chunkX[threadIdx.x]);
  /* block Level reduction */
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i){
      chunkW[threadIdx.x]=fmax(chunkW[threadIdx.x],chunkW[threadIdx.x+i]);
    }
    __syncthreads();
    i/=2;
  }
  if(threadIdx.x==0) scratch[blockIdx.x]=chunkW[0];
}


#undef __FUNCT__
#define __FUNCT__ "VecPointwiseDivide_SeqGPU"
PetscErrorCode VecPointwiseDivide_SeqGPU(Vec w,Vec x,Vec y){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid, dimBlock;
  #if(DEBUGVEC && VERBOSE)
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
  while(dimGrid.x>MAXBLOCKS){
    dimGrid.x/=2;
    dimBlock.x*=2;
  }
  cudaDeviceSynchronize();
  kernPDIV<<<dimGrid,dimBlock,3*dimBlock.x*sizeof(double)>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  #if(DEBUGVEC) 
     ierr = VecCheckCUDAError("kernel call to kernPDIV");CHKERRQ(ierr); 
  #endif
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
  double* chunkW = sharedPDIV + 2*blockDim.x;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    if(chunkX[threadIdx.x]*chunkY[threadIdx.x]!=0){/* using intrinsic div op */
      chunkW[threadIdx.x]=chunkX[threadIdx.x]/chunkY[threadIdx.x];
    }else{
      chunkW[threadIdx.x]=0;
    }
    devW[tid]=chunkW[threadIdx.x];
  }
}

/*--------------------------- end pointwise ---------------------------- */


/*-------------------------- norm functions ---------------------------- */
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


#undef __FUNCT__
#define __FUNCT__ "VecNorm_SeqGPU"
PetscErrorCode VecNorm_SeqGPU(Vec x,NormType type,PetscReal* z){
  /* NormType: NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4 */
  /* dealing with NORM_2 for now... */
  PetscFunctionBegin;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecNorm_SeqGPU\n");
  #endif
  PetscErrorCode ierr;
  double *devScratch,*devPartial,*hostPartial,zhost;
  PetscInt i,chunks=0,segment,partialsize,secondPhase,scratchsize;
  cudaStream_t* nrmstream;
  dim3 dimGrid, dimBlock;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  if(xd->syncState==VEC_CPU){
    #if(DEBUGVEC && VERBOSE)
       printf("xd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }

  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) /(float)(CHUNKWIDTH));
  nrmstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));
  /* make sure the segment size for each chunk is correct */
  if(chunks>1) segment = (int) (CHUNKWIDTH);
  else segment = x->map->n;
  dimGrid.x=ceil(((float)segment)/(float)THRNRMCNT);
  dimBlock.x  = THRNRMCNT;
  /* allocate gridwide scratch array */
  scratchsize = chunks*dimGrid.x*sizeof(double);
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize);
  ccs[0] = cudaMemsetAsync(devScratch,0,scratchsize,xd->stream);
  ccs[1]=cudaMemcpy(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice);
  #if(DEBUGVEC)
    #if(VERBOSE)
      printf("NORM: chunks: %d, seg: %d, blks: %d, scr: %d, chksize: %d, N: %d\n",
            chunks,segment,dimGrid.x,scratchsize/8,(int)(CHUNKWIDTH),x->map->n);
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecNorm_SeqGPU"); CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecNorm_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[1],"on copy segment length H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();
  for(i=0;i<chunks;i++){/* streaming async kernel calls */
    cudaStreamCreate(&(nrmstream[i]));
    cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,nrmstream[i]);
    /* Overlapping execution */
    kernNorm2<<<dimGrid,dimBlock,0,nrmstream[i]>>>(xd->devptr,
                                                   xd->segment,
                                                   xd->length,
                                                   xd->offset,
                                                   devScratch+i*dimGrid.x);
  }/* end for-loop */

  secondPhase = scratchsize/sizeof(double);
  if(secondPhase>1){ /* begin next reduction */
    dimGrid.x = ceil((float)secondPhase/(float)THRNRMCNT2);
    dimBlock.x  = THRNRMCNT2;
    /* allocate last reduction array */
    partialsize = dimGrid.x*sizeof(double);
    cms[1] = cudaMalloc((void**)&devPartial,partialsize);
    ccs[2] = cudaMemsetAsync(devPartial,0,partialsize,xd->stream);
    ccs[3] = cudaMemcpy(xd->segment,&secondPhase,sizeof(int),cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
      #if(VERBOSE)
         printf("NORM phase2: blks: %d, partial: %d\n",dimGrid.x,partialsize/sizeof(double));
      #endif
      ierr = VecCheckCUDAStatus(cms[1],"devPartial alloc in VecNorm_SeqGPU"); CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(ccs[2],"devPartial memset in VecNorm_SeqGPU");CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(ccs[3],"on copy chunks H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    kernRedNorm<<<dimGrid,dimBlock,dimBlock.x*sizeof(double)>>>(xd->segment,devScratch,devPartial);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAError("kernRedNorm_double launch in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
    ierr = PetscMalloc(partialsize,&hostPartial);CHKERRQ(ierr);

    /* Copy back norm z */
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    ccs[4]=cudaMemcpy(hostPartial,devPartial,partialsize,cudaMemcpyDeviceToHost);/* copy back */
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(ccs[4],"on devPartial copy D2H");CHKERRQ(ierr);
    #endif
    /* final reduction */
    if(dimGrid.x>1){
     zhost=0.;
     for(i=0;i<dimGrid.x;i++)zhost+=hostPartial[i];
    }else{
     zhost=hostPartial[0];
    }
    ierr = PetscFree(hostPartial); CHKERRQ(ierr);
    cms[2] = cudaFree(devPartial);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(cms[2],"on cudaFree(devPartial)");CHKERRQ(ierr);
    #endif
  }else{/* only copy back necessary */
     ccs[4]=cudaMemcpy(&zhost,devScratch,sizeof(double),cudaMemcpyDeviceToHost);/* copy back */
     #if(DEBUGVEC)
        ierr = VecCheckCUDAStatus(ccs[4],"on zhost copy D2H");CHKERRQ(ierr);
     #endif
  }
  *z = PetscSqrtScalar(zhost);


  /* clean up resources */
  for(i=0;i<chunks;i++){
    cudaStreamDestroy(nrmstream[i]);
  }
  free(nrmstream);
  cms[3] = cudaFree(devScratch);
  #if(DEBUGVEC)
     #if(VERBOSE)
        printf("Znorm: %e\n",*z);
     #endif
     ierr = VecCheckCUDAStatus(cms[3],"on cudaFree(devScratch)");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}

extern __shared__ double sharedRedNorm[];
#undef __FUNCT__
#define __FUNCT__ "kernRedNorm"
__global__ void kernRedNorm(int* size,double* scratch,double* z){/* reduction kernel */
  int i = (blockDim.x+1)/2;
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  double* zptr = sharedRedNorm;
  zptr[threadIdx.x]=(tid<*size)?scratch[tid]:0.;
  //printf("zptr[%d]: %e\n",threadIdx.x,zptr[threadIdx.x]);
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i) zptr[threadIdx.x]+=zptr[threadIdx.x+i];
    __syncthreads();
    i/=2;
  }/* end while */
  if(threadIdx.x==0){
    z[blockIdx.x]=zptr[0];
    //printf("ZNorm block[%d]: %e\n",blockIdx.x,z[blockIdx.x]);
  }
}


#undef __FUNCT__
#define __FUNCT__ "kernNorm2"
__global__ void kernNorm2(double* devX,int* segmentsize,int* arrsize,
                          int* offset,double *scratch){
  __shared__ double chunkX[THRNRMCNT];
  __shared__ int n;    n   = *arrsize;

  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  int item = *segmentsize**offset+tid;

  /* read in values to shared */
  chunkX[threadIdx.x]=(item<n)?devX[item]:0.;
  chunkX[threadIdx.x]*=chunkX[threadIdx.x];
  __syncthreads();

  /* block level reduction */
  while(i>0){
     if(threadIdx.x<i){
       chunkX[threadIdx.x]+=chunkX[threadIdx.x+i];
     }
     __syncthreads();
     i/=2;
  }/* end while */
  if(threadIdx.x==0) scratch[blockIdx.x]=chunkX[0];
}


/*
#undef __FUNCT__
#define __FUNCT__ "VecNorm1_SeqGPU"
PetscErrorCode VecNorm1_SeqGPU(Vec xin,NormType type,PetscReal* z)
{*/

/* NormType: NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4 */
/* dealing with NORM_2 for now... */
/* z has 2 elements */

/*
  PetscErrorCode ierr;
  PetscFunctionBegin;
  printf("Call to VecNorm_SeqGPU\n");
  ierr = VecDot_SeqGPU(xin,xin,&z[0]);CHKERRQ(ierr);
  z[0]=PetscSqrtScalar(z[0]);
  printf("ZNORM: %f\n\n",*z);
  PetscFunctionReturn(0);
}*/


/*
#undef __FUNCT__
#define __FUNCT__ "kernReduceAbsSum"
PetscErrorCode kernReduceAbsSum(double * x, PetscReal* z){

}
*/
/* ------------------------------ end norms -------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "VecGetArray_SeqGPU"
PetscErrorCode VecGetArray_SeqGPU(Vec v,PetscScalar **a){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU *vd=(Vec_SeqGPU*)v->data;
  if(vd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"in VecGetArray_SeqGPU the vector has not been created.");
  }
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecGetArray_SeqGPU\n");
  #endif
  PetscInt flg1=0,flg2=0;
  PetscStackCheckByName(4,"DMDAVecGetArray",flg1);
  PetscStackCheckByName(6,"DMGlobalToLocalBegin",flg2);

  if((flg1 || flg2) && vd->syncState==VEC_GPU){
    ierr = VecCopyOverD2H(v,vd->cpuptr); CHKERRQ(ierr);
    vd->syncState = VEC_CPU;
  }
  cudaDeviceSynchronize();
  *a=vd->cpuptr;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecRestoreArray_SeqGPU"
PetscErrorCode VecRestoreArray_SeqGPU(Vec v,PetscScalar **a){
  PetscFunctionBegin;
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecRestoreArray_SeqGPU\n");
  #endif
  PetscErrorCode ierr;
  Vec_SeqGPU *vd=(Vec_SeqGPU*)v->data;
  PetscInt flg1=0;
  PetscStackCheckByName(1,"VecRestoreArrayRead",flg1);
  if(!flg1){
    if(a){
      ierr = VecCopyOverH2D(v,*a);CHKERRQ(ierr);
      vd->syncState=VEC_GPU;
    }else{
      ierr = VecCopyOverH2D(v,vd->cpuptr);CHKERRQ(ierr);
      vd->syncState=VEC_SYNCHED;
    }
  }
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
    ierr = PetscMemcpy((void*)dd->cpuptr,(void*)sd->cpuptr,s->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = VecCopyOverDevice(d,s); CHKERRQ(ierr);
  dd->syncState=sd->syncState;/* synch signal copy */
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecCopy_SeqGPU\n");
  #endif
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecSwap_SeqGPU"
PetscErrorCode VecSwap_SeqGPU(Vec xin,Vec yin){
  /* PetscErrorCode ierr; */
  PetscFunctionBegin;
  printf("VecSwap_SeqGPU (***EMPTY***)\n");
  if (xin != yin) {
#if defined(PETSC_USE_REAL_SINGLE)
    //////// cublasSswap(bn,VecCUSPCastToRawPtr(*xarray),one,VecCUSPCastToRawPtr(*yarray),one);
#else
    //////   cublasDswap(bn,VecCUSPCastToRawPtr(*xarray),one,VecCUSPCastToRawPtr(*yarray),one);
#endif

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
#if(DEBUGVEC && VERBOSE)
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
  #if(DEBUGVEC && VERBOSE)
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
#if(DEBUGVEC && VERBOSE)
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
#define __FUNCT__ "PinnedMalloc"
static PetscErrorCode  PinnedMalloc(PetscScalar** x,PetscInt n){
  PetscFunctionBegin;
#if(DEBUGVEC && VERBOSE)
     printf("Call to PinnedMalloc\n"); 
  #endif
  cms[0]=cudaHostAlloc((void**)x,n,0);
  #if(DEBUGVEC)
     PetscErrorCode ierr;
     ierr=VecCheckCUDAStatus(cms[0],"in PinnedMalloc");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PinnedFree"
static PetscErrorCode  PinnedFree(PetscScalar* x){
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
  /* V->ops->aypx            = VecAYPX_SeqGPU; */
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
  V->petscnative=PETSC_FALSE;
  seqgpu->syncState      = VEC_UNALLOC;
  seqgpu->unplacedarray=PETSC_NULL;
  seqgpu->array_allocated=PETSC_NULL;
  seqgpu->array=PETSC_NULL;
  /* create an associated stream */
  cms[0] = cudaStreamCreate(&(seqgpu->stream));
  /* allocate the variable for vector size */
  cms[1]=cudaMalloc((void**)&(seqgpu->length),sizeof(int));
  /* send vec length size to device */
  ccs[0]=cudaMemcpyAsync((void*)seqgpu->length,
               (void*)&(V->map->n),sizeof(int),cudaMemcpyHostToDevice,seqgpu->stream);
  /* allocate the vector on device */
  cms[2]=cudaMalloc((void**)&(seqgpu->devptr),V->map->n*sizeof(double));
  ccs[1]=cudaMemsetAsync((void*)seqgpu->devptr,0,V->map->n*sizeof(double),seqgpu->stream);
  /* allocate the variable for vector offsets */
  cms[3]=cudaMalloc((void**)&(seqgpu->offset),sizeof(int));
  /* allocate the variable for vector segment length */
  cms[4]=cudaMalloc((void**)&(seqgpu->segment),sizeof(int));
  /* allocate the variable for vector single value result */
  cms[5]=cudaMalloc((void**)&(seqgpu->zval),sizeof(double));
  cms[6]=cudaMalloc((void**)&(seqgpu->scalar),sizeof(double));
  /* using pinned memory */
  ierr = PinnedMalloc(&(seqgpu->cpuptr),V->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
  //ierr = PetscMalloc(V->map->n*sizeof(PetscScalar),&(seqgpu->cpuptr));
  ierr = PetscMemzero(seqgpu->cpuptr,V->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
  seqgpu->syncState=VEC_ALLOC;


  #if(DEBUGVEC)
    #if(VERBOSE)
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
#if(DEBUGVEC && VERBOSE)
     printf("Call to VecDestroyArray_SeqGPU\n"); 
  #endif
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  if(vd && vd->syncState != VEC_UNALLOC){
      cms[0]=cudaFree(vd->devptr);  vd->devptr=PETSC_NULL;
      cms[1]=cudaFree(vd->length);  vd->length=PETSC_NULL;
      cms[2]=cudaFree(vd->segment); vd->segment=PETSC_NULL;
      cms[3]=cudaFree(vd->zval);    vd->zval=PETSC_NULL;
      cms[4]=cudaFree(vd->scalar);  vd->scalar=PETSC_NULL;
      cms[5] = cudaStreamDestroy(vd->stream);
      ierr = PinnedFree(vd->cpuptr); CHKERRQ(ierr);
      //ierr = PetscFree(vd->cpuptr);CHKERRQ(ierr);
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
  #if(DEBUGVEC && VERBOSE)
     printf("Call to VecDestroyVecs_SeqGPU\n");
  #endif
  PetscErrorCode ierr;
  PetscInt i;
   /* destroy the internal part */
  for(i=0;i<m;i++){
    ierr = VecDestroy(&vv[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecView_Seq_ASCII"
static PetscErrorCode VecView_Seq_ASCII(Vec xin,PetscViewer viewer){
  PetscFunctionBegin;
  printf("VecView_Seq_ASCII() (***EMPTY***)\n");
  PetscFunctionReturn(0);
}

EXTERN_C_END

