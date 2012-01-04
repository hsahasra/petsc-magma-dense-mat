#include <petscconf.h>
#include <petscsys.h>
PETSC_CUDA_EXTERN_C_BEGIN
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include <float.h>
#include <private/vecimpl.h>          /*I "petscvec.h" I*/
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
  static PetscBool seed_flag=PETSC_TRUE;
  PetscErrorCode ierr;
  PetscInt i,bx,tx;
  uint *seeds=PETSC_NULL,*devseeds=PETSC_NULL;
  PetscScalar rval;
  dim3 dimBlock,dimGrid;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  /* assuming width mem load isn't going to be an issue */
  /* printf("Call to VecSetRandom_SeqGPU\n");*/
  if(xd->syncState==VEC_ALLOC || xd->syncState==VEC_CPU){
    for(i=0; i<x->map->n; i++){
       ierr = PetscRandomGetValue(r,&xd->cpuptr[i]);CHKERRQ(ierr);
    }
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }else if(xd->syncState==VEC_SYNCHED || xd->syncState==VEC_GPU){
    bx=ceil((float)x->map->n/(float)TCOUNT);
    ierr = PetscMalloc(bx*sizeof(PetscInt),&seeds);CHKERRQ(ierr);
    tx=TCOUNT;
    dimGrid.x=bx; dimGrid.y=1;
    dimBlock.x=tx; dimBlock.y=1;
    if(seed_flag){
      for(i=0; i<bx; i++){
         ierr = PetscRandomGetValue(r,&rval);CHKERRQ(ierr);
         seeds[i]=(uint)(UINT_MAX*rval);
      }

      cms[0] = cudaMalloc((void**)&devseeds,bx*sizeof(uint));
      ccs[0]=cudaMemcpy(devseeds,seeds,bx*sizeof(uint),cudaMemcpyHostToDevice);
      #if(DEBUGVEC)
        ierr = VecCheckCUDAStatus(cms[0],"error in cudaMalloc");CHKERRQ(ierr);
        ierr = VecCheckCUDAStatus(ccs[0],"on copy H2D in VecSetRandom_SeqGPU");CHKERRQ(ierr);
      #endif

      kernRandS<<<dimGrid,dimBlock>>>(devseeds);
      #if(DEBUGVEC)
        ierr = VecCheckCUDAError("kernRandS launch");CHKERRQ(ierr);
      #endif

      ierr = PetscFree(seeds);CHKERRQ(ierr);
      cudaDeviceSynchronize();
      cms[1] = cudaFree(devseeds);
      #if(DEBUGVEC)
         ierr = VecCheckCUDAStatus(cms[1],"in cudaFree()");CHKERRQ(ierr);
      #endif

      seed_flag=PETSC_FALSE;
    }
    kernRand<<<dimGrid,dimBlock>>>(xd->devptr,xd->length);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernRand launch");CHKERRQ(ierr);
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
  int blocks,threads;/* assuming shared memory size is not an issue */
  if(blocksize && !offset){
    blocks=ceil((float)blocksize/(float)TCOUNT);
  } else {
    blocks=ceil((float)x->map->n/(float)TCOUNT);
    blocksize = x->map->n;
  }
  threads=TCOUNT;
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

  dim3 dimGrid; dimGrid.x=blocks; dimGrid.y=1;
  dim3 dimBlock; dimBlock.x=threads; dimBlock.y=1;
  kernCompare<<<dimGrid,dimBlock>>>(xd->devptr,yd->devptr,xd->length,yd->length,devsame);
  ierr = VecCheckCUDAError("kernCompare launch");CHKERRQ(ierr);

  cudastatus=cudaMemcpy(&cpusame,devsame,sizeof(int),cudaMemcpyDeviceToHost);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H in VecCompare_SeqGPU");CHKERRQ(ierr);

  if(cpusame==1)*same=PETSC_TRUE;
  else *same=PETSC_FALSE;
  cudastatus = cudaFree(devsame);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernCompare"
__global__ void kernCompare(double* devX, double* devY, int* lx, int* ly, int* devsame){

  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int2 localOBS = integer2Symbol;
  int localn = localOBS.x+localOBS.y;
  int index = tid+localOBS.x;
  double value=0;
  __shared__ unsigned char blockflag;
  __shared__ double chunkX[TCOUNT];
  __shared__ double chunkY[TCOUNT];

  if(threadIdx.x==0)blockflag=0;
  __syncthreads();
  if(index<localn){
    /* read in values to shared */
    chunkX[threadIdx.x]=devX[index];
    chunkY[threadIdx.x]=devY[index];
    value = fabs(chunkX[threadIdx.x]-chunkY[threadIdx.x]);
    if(value>1e-16){
      //printf("In kernCompare found an element mismatch: %e\n",value);
      blockflag=1;
    }
    if(*lx!=*ly){
      //printf("In kernCompare found length mismatch: lx: %d vs ly: %d\n",*lx,*ly);
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

/*------------------------------- end compare --------------------------*/


/*---------------------------- Vec info functions ----------------------*/

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
  ierr= PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);/* forced ASCII for now */
  ierr =VecView_Seq_ASCII(x,viewer);CHKERRQ(ierr);
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
    printf("kernCheck: x[%d]: %e, length: %d\n",tid,x[tid],*n);
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
  ccs[0]=cudaMemcpyAsync(vd->devptr,y,
               v->map->n*sizeof(PetscScalar),cudaMemcpyHostToDevice,vd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Call to VecCopyBlockH2D\n");
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
#undef __FUNCT__
#define __FUNCT__ "VecSetValues_SeqGPU"
/*@
   VecSetValues - Inserts or adds values into certain locations of a vector.
@*/
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
      #pragma omp parallel for
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



#undef __FUNCT__
#define __FUNCT__ "VecSet_SeqGPU"
PetscErrorCode VecSet_SeqGPU(Vec xin,PetscScalar alpha){
  PetscFunctionBegin;
  #if(DEBUGVEC)
    PetscErrorCode ierr;
  #endif
  dim3 dimgrid(ceil((float)xin->map->n/((float)TCOUNT)),1,1);
  dim3 dimblocks(TCOUNT,1,1);
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
    kernSet<<<dimgrid,dimblocks>>>(xd->devptr,xd->length);
    #if(DEBUGVEC) 
       ierr = VecCheckCUDAError("Call to kernSet. "); CHKERRQ(ierr);
    #endif
    xd->syncState=VEC_GPU;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernSet"
__global__ void kernSet(double* x, int* n){
  int tid = threadIdx.x + blockDim.x*blockIdx.x;
  __shared__ double chunkX[TCOUNT];
  chunkX[threadIdx.x] = dblScalarValue;
  if(tid<*n){
    x[tid] = chunkX[threadIdx.x]; /* arr[threadIdx.x]; */
  }
}




#undef __FUNCT__
#define __FUNCT__ "VecScale_SeqGPU"
PetscErrorCode VecScale_SeqGPU(Vec x, PetscScalar alpha){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  dim3 dimgrid(ceil((float)x->map->n/((float)TCOUNT)),1,1);
  dim3 dimblocks(TCOUNT,1,1);
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
    kernScale<<<dimgrid,dimblocks,0,xd->stream>>>(xd->devptr,xd->length);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("Call to kernScale. "); CHKERRQ(ierr); 
    #endif
  }
  xd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernScale"
__global__ void kernScale(double* x, int* n){
  int tid = threadIdx.x + blockDim.x*blockIdx.x;
  __shared__ double arr[TCOUNT];
  double localdbl=dblScalarValue;
  if(tid<*n){
    arr[threadIdx.x] = x[tid];
    arr[threadIdx.x] *= localdbl;
    x[tid] = arr[threadIdx.x];
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
  PetscErrorCode ierr;
  double *devScratch,*devPartial,zhost;
  PetscInt i,chunks=0,*devChunks,segment,partialsize,scratchsize;
  cudaStream_t* dotstream;
  dim3 dimGrid, dimBlock;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
 
  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) /(float)(CHUNKWIDTH));
  dotstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));

  if(chunks>1){
    segment = (int) CHUNKWIDTH;
    dimGrid.x=ceil((CHUNKWIDTH)/(float)THRDOTCNT);
  }else{
    segment = x->map->n;
    dimGrid.x=ceil(((float)segment)/(float)THRDOTCNT);
  }
  if(dimGrid.x&1)dimGrid.x++;/* make sure even number of blocks */
  dimBlock.x = THRDOTCNT;


  /* set up on x stream */
  if(xd->syncState==VEC_CPU){
    if(DEBUGVEC && VERBOSE)printf("xd state VEC_CPU: copying to device.\n");
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  scratchsize=chunks*dimGrid.x*sizeof(double);
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize);/* scratch pad */
  ccs[0] = cudaMemsetAsync(devScratch,0,scratchsize,xd->stream); 

  /* set up on y stream */
  if(yd->syncState==VEC_CPU){
    #if(DEBUGVEC && VERBOSE)
       printf("yd state VEC_CPU: copying to device.\n");
    #endif
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  partialsize=chunks*sizeof(double);
  cms[1]=cudaMalloc((void**)&devPartial,partialsize);/* partial results to be combined */
  ccs[1]=cudaMemsetAsync(devPartial,0,partialsize,yd->stream);
  ccs[4]=cudaMemcpyAsync(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice,yd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
      printf("Call to VecDot_SeqGPU\n");
      printf("DOT chunks: %d, seg: %d, blks: %d, par: %d, scr: %d, chkwidth: %d\n",chunks,xd->segment,dimGrid.x,partialsize,scratchsize,(int)(CHUNKWIDTH));
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecDot_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecDot_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[1],"devPartial alloc in VecDot_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[1],"devPartial memset in VecDot_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[4],"on copy segment size H2D in VecDot_SeqGPU");CHKERRQ(ierr);
  #endif

  cudaDeviceSynchronize();/* make sure everyone is ready */
  for(i=0;i<chunks;i++){  /* streaming async kernel calls */
    ccs[2]=cudaStreamCreate(&(dotstream[i]));
    ccs[3]=cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,dotstream[i]);

    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(ccs[2],"on cudaStreamCreate");CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(ccs[3],"on copy array length H2D in VecDot_SeqGPU");CHKERRQ(ierr);
    #endif

    /* Overlapping execution */
    kernDot<<<dimGrid,dimBlock,0,dotstream[i]>>>(xd->devptr,yd->devptr,
                                                          xd->segment,
                                                          xd->length,
                                                          xd->offset,
                                                          (devScratch+i*dimGrid.x),
                                                          (devPartial+i));
    #if(DEBUGVEC)
      ierr = VecCheckCUDAError("kernDot launch in VecDot_SeqGPU");CHKERRQ(ierr);
    #endif
  }
  /* dot product block reduction */
  if(chunks>1){
    dimGrid.x  = 1;
    dimBlock.x = chunks;
    if(dimBlock.x&1)dimBlock.x++;
    cms[2] = cudaMalloc((void**)&devChunks,sizeof(int));
    ccs[5] = cudaMemcpyAsync(devChunks,&chunks,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
    #if(DEBUGVEC)
        ierr = VecCheckCUDAStatus(cms[2],"on cudamalloc()");CHKERRQ(ierr);
        ierr = VecCheckCUDAStatus(ccs[5],"on copy chunks H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    kernRedDot<<<dimGrid,dimBlock,chunks*sizeof(double)>>>(devPartial,devChunks,xd->zval);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernRedDot launch in VecDot_SeqGPU");CHKERRQ(ierr); 
    #endif
    ccs[6]=cudaMemcpy(&zhost,xd->zval,sizeof(double),cudaMemcpyDeviceToHost);/* copy back z */
    cms[3]=cudaFree(devChunks);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(ccs[6],"on copy zval D2H in VecDot_SeqGPU");CHKERRQ(ierr); 
      ierr = VecCheckCUDAStatus(cms[3],"on cudaFree()");CHKERRQ(ierr);
    #endif
  }else{
    ccs[7]=cudaMemcpy(&zhost,devPartial,sizeof(double),cudaMemcpyDeviceToHost);/* copy back z */
    #if(DEBUGVEC) 
       ierr = VecCheckCUDAStatus(ccs[7],"on copy devP D2H in VecDot_SeqGPU");CHKERRQ(ierr);
    #endif
  }
  *z=zhost;

  /* clean up resources */
  for(i=0;i<chunks;i++){
     cudaStreamDestroy(dotstream[i]);
  }
  free(dotstream);
  cms[4] = cudaFree(devPartial);
  cms[5] = cudaFree(devScratch);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Zdot: %e\n",*z);
    #endif
    ierr = VecCheckCUDAStatus(cms[4],"on cudaFree()");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[5],"on cudaFree()");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


extern __shared__ double zDot[];
#undef __FUNCT__
#define __FUNCT__ "kernRedDot"
__global__ void kernRedDot(double* arr,int* chunks, double* z){/* reduction kernel */

  int i = (blockDim.x+1)/2;
  zDot[threadIdx.x]=0.;
  if(threadIdx.x<*chunks)zDot[threadIdx.x]=arr[threadIdx.x];
  while(i>0){
    if(threadIdx.x<i){
      zDot[threadIdx.x]+=zDot[threadIdx.x+i];
    }
    __syncthreads();
    i/=2;
  }
  if(threadIdx.x==0){
    *z=zDot[0];
  }
}



#undef __FUNCT__
#define __FUNCT__ "kernDot"
__global__ void kernDot(double* devX, double* devY,
                        int* segmentsize, int* arrsize,
                        int* offset, double* scratch, double* z){
  __shared__ double chunkX[THRDOTCNT];
  __shared__ double chunkY[THRDOTCNT];
  __shared__ int n;    n   = *arrsize;
  __shared__ int seg;  seg = *segmentsize;
  __shared__ int off;  off = *offset;

  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  int j = (gridDim.x+1)/2;
  int item = seg*off+tid;

  if(item<n){
    /* read in values to shared */
    chunkX[threadIdx.x]=devX[item]; /* offset values */
    chunkY[threadIdx.x]=devY[item]; /* offset values */
  }else{
    chunkX[threadIdx.x]=0.;
    chunkY[threadIdx.x]=0.;
  }

  chunkX[threadIdx.x]*=chunkY[threadIdx.x];
  __syncthreads();

  /* block level reduction */
  while(i>0){
     if(threadIdx.x<i){
       chunkX[threadIdx.x]+=chunkX[threadIdx.x+i];
     }
     __syncthreads();
     i/=2;
  }/* end while */

  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunkX[0];
  }
  __syncthreads();


  /* grid level reduction */
  while(j>0){
    if(threadIdx.x==0 && blockIdx.x<j){
      scratch[blockIdx.x]+=scratch[blockIdx.x+j];
    }
    __syncthreads();
    j/=2;
  }
  if(tid==0)*z=scratch[blockIdx.x];
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
  PetscInt bx,tx;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)w->data;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid;
  dim3 dimBlock;
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
  cudaDeviceSynchronize();
  /* assuming width mem load isn't going to be an issue */
  bx=ceil((float)y->map->n/(float)AXPYTCOUNT);
  tx=AXPYTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;

  if(alpha==0.0){
    ierr = VecCopyOverDevice(w,y);CHKERRQ(ierr);
  }else if(alpha==1.0){
    kernWXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWXPY");CHKERRQ(ierr); 
    #endif
  }else if(alpha==-1.0){
    kernWXMY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWXMY");CHKERRQ(ierr);
    #endif
  }else{
    ccs[0]=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
    #endif
    kernWAXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernWAXPY");CHKERRQ(ierr); 
    #endif
  }
  wd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernWAXPY"
__global__ void  kernWAXPY(double* devY,double* devX, int* vlen, double* devW){

 /* w <- y + alpha*x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double alphaShared;
  __shared__ double chunkY[AXPYTCOUNT];
  __shared__ double chunkX[AXPYTCOUNT];
  __shared__ double chunkW[AXPYTCOUNT];

  alphaShared = dblScalarValue;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]+(chunkX[threadIdx.x]*alphaShared);
    devW[tid]=chunkW[threadIdx.x];
  }
}

#undef __FUNCT__
#define __FUNCT__ "kernWXPY"
__global__ void  kernWXPY(double* devY,double* devX, int* vlen, double* devW){

 /* w <- y + x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunkY[AXPYTCOUNT];
  __shared__ double chunkX[AXPYTCOUNT];
  __shared__ double chunkW[AXPYTCOUNT];
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]+chunkX[threadIdx.x];
    devW[tid]=chunkW[threadIdx.x];
  }
}

#undef __FUNCT__
#define __FUNCT__ "kernWXMY"
__global__ void  kernWXMY(double* devY,double* devX, int* vlen, double* devW){

 /* w <- y + alpha*x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunkY[AXPYTCOUNT];
  __shared__ double chunkX[AXPYTCOUNT];
  __shared__ double chunkW[AXPYTCOUNT];
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
  PetscInt i;  PetscInt bx,tx;
  dim3 dimGrid;
  dim3 dimBlock;
  PetscScalar *devW;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=PETSC_NULL;

  cms[0] = cudaMalloc((void**)&devW,x->map->n*sizeof(double));
  ccs[0] = cudaMemset(devW,0,x->map->n*sizeof(double));

  /* assuming xwidth mem load isn't going to be an issue */
  bx=ceil((float)x->map->n/(float)AXPYTCOUNT);
  tx=AXPYTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;

  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Number of vectors in MAXPY: %d, bx: %d, tx: %d\n",nv,bx,tx);
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
    ccs[1]=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha[i],sizeof(double),0,cudaMemcpyHostToDevice);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[1],"error in symbol copy to device");CHKERRQ(ierr); 
    #endif
    cudaDeviceSynchronize();
    if(alpha[i]==0){
      continue;
    }else if(alpha[i]==1.){
      /* assuming width mem load isn't going to be an issue */
      kernXPY<<<dimGrid,dimBlock>>>(devW,yd->devptr,yd->length);
      #if(DEBUGVEC)
        ierr = VecCheckCUDAError("kernel call to kernXPY");CHKERRQ(ierr);
      #endif
    }else{
      /* assuming width mem load isn't going to be an issue */
      kernAXPY<<<dimGrid,dimBlock>>>(devW,yd->devptr,yd->length);
      #if(DEBUGVEC)
         ierr = VecCheckCUDAError("kernel call to kernAXPY");CHKERRQ(ierr);
      #endif
    }
  }
  if(xd->syncState==VEC_CPU){/* synch x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  cudaDeviceSynchronize();
  kernXPY<<<dimGrid,dimBlock>>>(xd->devptr,devW,xd->length);
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

#undef __FUNCT__
#define __FUNCT__ "kernXPY"
__global__ void  kernXPY(double* devY,double* devX, int* vlen){

 /* y <- y + x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;

  __shared__ double chunkY[AXPYTCOUNT];
  __shared__ double chunkX[AXPYTCOUNT];

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
  PetscInt bx,tx;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid;
  dim3 dimBlock;
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
  ccs[0]=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
  #if(DEBUGVEC)
   #if(VERBOSE)
      printf("VecAXPY_SeqGPU\n");
   #endif
   ierr = VecCheckCUDAStatus(ccs[0],"error in symbol copy to device");CHKERRQ(ierr);
  #endif
  tx=AXPYTCOUNT; bx=ceil((float)x->map->n/(float)AXPYTCOUNT);
  dimGrid.x=bx;  dimBlock.x=tx;
  cudaDeviceSynchronize();
  if(alpha==1.){
    kernXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,yd->length);
  }else if(alpha!=0){
    kernAXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,yd->length);
  }
  #if(DEBUGVEC)
   ierr = VecCheckCUDAError("kernel call in VecAXPY_SeqGPU");CHKERRQ(ierr);
  #endif
  yd->syncState=VEC_GPU;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernAXPY"
__global__ void  kernAXPY(double* devY,double* devX, int* vlen){

 /* y <- y + alpha*x */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double alphaShared;
  __shared__ double chunkY[AXPYTCOUNT];
  __shared__ double chunkX[AXPYTCOUNT];
  alphaShared = dblScalarValue;
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkY[threadIdx.x]+=chunkX[threadIdx.x]*alphaShared;
    devY[tid]=chunkY[threadIdx.x];
  }
}

#undef __FUNCT__
#define __FUNCT__ "VecAXPBYPCZ_SeqGPU"
PetscErrorCode VecAXPBYPCZ_SeqGPU(Vec x, PetscScalar alpha, PetscScalar beta,\
                           PetscScalar gamma, Vec y, Vec z){

  PetscFunctionBegin;
  #if(DEBUGVEC)
     PetscErrorCode ierr;
  #endif
  int blocks=ceil((float)x->map->n/(float)AXPBYPCZTCOUNT);/* assuming shared memory size is not an issue */
  int threads=AXPBYPCZTCOUNT;
  Vec_SeqGPU* devX = (Vec_SeqGPU*)x->data;
  Vec_SeqGPU* devY = (Vec_SeqGPU*)y->data;
  Vec_SeqGPU* devZ = (Vec_SeqGPU*)z->data;

  double2 alphabeta;
  alphabeta.x = alpha;
  alphabeta.y = beta;
  dim3 dimGrid; dimGrid.x=blocks; dimGrid.y=1;
  dim3 dimBlock; dimBlock.x=threads; dimBlock.y=1;

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

#undef __FUNCT__
#define __FUNCT__ "kernAXPBYPCZ"
__global__ void kernAXPBYPCZ(double* devX, double* devY, double* devZ, int* len){
  /* x <- alpha*x + beta*y + gamma*z */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int localn = *len;

  __shared__ double work[AXPBYPCZTCOUNT];
  __shared__ double chunkX[AXPBYPCZTCOUNT];
  __shared__ double chunkY[AXPBYPCZTCOUNT];
  __shared__ double chunkZ[AXPBYPCZTCOUNT];

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
  if(DEBUGVEC && VERBOSE)printf("VecPointwiseMult_SeqGPU\n");
  PetscErrorCode ierr;
  PetscInt bx,tx;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid;
  dim3 dimBlock;
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
  bx=ceil((float)y->map->n/(float)PMULTCOUNT);
  tx=PMULTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;

  cudaDeviceSynchronize();
  kernPMULT<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  #if(DEBUGVEC)
     ierr = VecCheckCUDAError("kernel call to kernPMULT");CHKERRQ(ierr);
  #endif

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernPMULT"
__global__ void  kernPMULT(double* devY,double* devX, int* vlen, double* devW){

 /* w <- x./y */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunkY[PMULTCOUNT];
  __shared__ double chunkX[PMULTCOUNT];
  __shared__ double chunkW[PMULTCOUNT];
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
  //printf("VecMaxPointwiseDivide_SeqGPU...");
  PetscErrorCode ierr;
  PetscScalar *devMax,*devScratch,*devPartial;
  PetscInt i,chunks=0,*devChunks,segment,partialsize,scratchsize;
  cudaStream_t* pwdstream;
  dim3 dimGrid, dimBlock;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;

  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) /(float)(CHUNKWIDTH));
  pwdstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));

  if(chunks>1){
    segment = (int) (CHUNKWIDTH);
    dimGrid.x=ceil((CHUNKWIDTH)/(float)PDIVTCOUNT);
  }else{
    segment = x->map->n;
    dimGrid.x=ceil(((float)segment)/(float)PDIVTCOUNT);
  }
  dimBlock.x  = PDIVTCOUNT;
  if(dimGrid.x&1)dimGrid.x++;

  if(x->map->n!=y->map->n){
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

  cms[0]=cudaMalloc((void**)&devMax,sizeof(PetscScalar));
  scratchsize = chunks*dimGrid.x*sizeof(double);
  cms[1] = cudaMalloc((void**)&devScratch,scratchsize);
  ccs[0] = cudaMemsetAsync(devScratch,0,scratchsize,yd->stream);

  partialsize = chunks*sizeof(PetscScalar);
  cms[2]=cudaMalloc((void**)&devPartial,partialsize);
  ccs[1] = cudaMemsetAsync(devPartial,0,partialsize,xd->stream);
  ccs[3]=cudaMemcpyAsync(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
  #if(DEBUGVEC)
    ierr = VecCheckCUDAStatus(cms[0],"devMax alloc in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[1],"devScratch alloc in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[2],"devPartial alloc in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[1],"devPartial memset in VecMPWD_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[3],"on copy segment size H2D in VecMPWD_SeqGPU");CHKERRQ(ierr);
  #endif

  cudaDeviceSynchronize();
  for(i=0;i<chunks;i++){
    cms[3]=cudaStreamCreate(&(pwdstream[i]));
    ccs[2]=cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,pwdstream[i]);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(cms[3],"on cudaStreamCreate");CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(ccs[2],"on copy array length H2D in VecMPWD_SeqGPU");CHKERRQ(ierr);
    #endif
    /* Overlapping execution */
    kernMAXPDIV<<<dimGrid,dimBlock,0,pwdstream[i]>>>(xd->devptr,yd->devptr,
                                                     xd->segment,
                                                     xd->length,
                                                     xd->offset,
                                                     (devScratch+i*dimGrid.x),
                                                     (devPartial+i));
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernMAXPDIV launch in VecMPWD_SeqGPU");CHKERRQ(ierr);
    #endif
  }

  if(chunks>1){
    dimGrid.x=1;
    dimBlock.x = chunks;
    if(dimBlock.x&1)dimBlock.x++;
    cms[4] = cudaMalloc((void**)&devChunks,sizeof(int));
    ccs[4] = cudaMemcpyAsync(devChunks,&chunks,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
    #if(DEBUGVEC)
        ierr = VecCheckCUDAStatus(cms[4],"on cudamalloc()");CHKERRQ(ierr);
        ierr = VecCheckCUDAStatus(ccs[4],"on copy chunks H2D in VecMPWD_SeqGPU");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();
    kernMAX<<<dimGrid,dimBlock,chunks*sizeof(double)>>>(devPartial,devChunks,devMax);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernel call to kernMAX");CHKERRQ(ierr);
    #endif
    ccs[5]=cudaMemcpy(max,devMax,sizeof(PetscScalar),cudaMemcpyDeviceToHost);/* copy back */
    cms[5]=cudaFree(devChunks);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(ccs[5],"on devMax copy D2H");CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(cms[5],"on cudaFree()");CHKERRQ(ierr);
    #endif

  }else{
    cudaDeviceSynchronize();
    ccs[6]=cudaMemcpy(max,devPartial,sizeof(PetscScalar),cudaMemcpyDeviceToHost);/* copy back */
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[6],"on devP copy D2H");CHKERRQ(ierr);
    #endif
  }

  for(i=0;i<chunks;i++){
    cudaStreamDestroy(pwdstream[i]);
  }
  cms[6] = cudaFree(devScratch);
  cms[7] = cudaFree(devPartial);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("max: %e\n",*max);
    #endif
    ierr = VecCheckCUDAStatus(cms[6],"on cudaFree()");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[7],"on cudaFree()");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}


extern __shared__ double mlist[];
#undef __FUNCT__
#define __FUNCT__ "kernMAX"
__global__ void  kernMAX(double* maxlist,int *chunks,double* max){
  int i,tid;
  tid = threadIdx.x;
  i = (blockDim.x+1)/2;
  mlist[threadIdx.x]=0.;
  if(tid<*chunks) mlist[tid]=maxlist[tid];
  __syncthreads();
  while(i>0){
    if(tid<i){
      mlist[tid] = (mlist[tid]>mlist[tid+i])?mlist[tid]:mlist[tid+i];
    }
    __syncthreads();
    i/=2;
  }
  if(tid==0)*max = mlist[0];
}



#undef __FUNCT__
#define __FUNCT__ "kernMAXPDIV"
__global__ void  kernMAXPDIV(double* devY,double* devX, int* segmentsize,
                             int* arrsize,int* offset,double* scratch,double* maxitem){
 /* w <- max(abs(x./y)) */
  __shared__ double chunkY[PDIVTCOUNT];
  __shared__ double chunkX[PDIVTCOUNT];
  __shared__ double chunkW[PDIVTCOUNT];
  __shared__ int n;    n   = *arrsize;
  __shared__ int seg;  seg = *segmentsize;
  __shared__ int off;  off = *offset;

  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  int j = (gridDim.x+1)/2;
  int item = seg*off+tid;

  if(item<n){
    chunkX[threadIdx.x]=devX[item];
    chunkY[threadIdx.x]=devY[item];
    if(chunkY[threadIdx.x]!=0){
      chunkW[threadIdx.x]=fabs(__ddiv_rn(chunkX[threadIdx.x],chunkY[threadIdx.x]));
      #if(DEBUGVEC && VERBOSE)
         printf("In kernMAXPDIV: chunkW[%d]: %e\n",threadIdx.x,chunkW[threadIdx.x]);
      #endif
    }else{
      chunkW[threadIdx.x]=fabs(chunkX[threadIdx.x]);
    }
  }else{
    chunkW[threadIdx.x]=0.0;
  }

  /* block Level reduction */
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i){
      chunkW[threadIdx.x]=(chunkW[threadIdx.x]>chunkW[threadIdx.x+i])?chunkW[threadIdx.x]:chunkW[threadIdx.x+i];
      #if(DEBUGVEC && VERBOSE)
         printf("In kernMAXPDIV: chunk2W[%d]: %e\n",threadIdx.x,chunkW[threadIdx.x]);
      #endif
    }
    __syncthreads();
    i/=2;
  }
  if(threadIdx.x==0)scratch[blockIdx.x]=chunkW[0];

  /* grid level reduction */
  while(j>0){
    if(threadIdx.x==0 && blockIdx.x<j){
      scratch[blockIdx.x]=(scratch[blockIdx.x]>scratch[blockIdx.x+j])?scratch[blockIdx.x]:scratch[blockIdx.x+j];
      #if(DEBUGVEC && VERBOSE)
         printf("In KernMAXPDIV: scratch[%d]: %e\n",blockIdx.x,scratch[blockIdx.x]);
      #endif
    }
    __syncthreads();
    j/=2;
  }
  if(tid==0)*maxitem=scratch[0];
}


#undef __FUNCT__
#define __FUNCT__ "VecPointwiseDivide_SeqGPU"
PetscErrorCode VecPointwiseDivide_SeqGPU(Vec w,Vec x,Vec y){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt bx,tx;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  Vec_SeqGPU *wd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid;
  dim3 dimBlock;
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
  bx=ceil((float)y->map->n/(float)PDIVTCOUNT);
  tx=PDIVTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;
  cudaDeviceSynchronize();
  kernPDIV<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  #if(DEBUGVEC) 
     ierr = VecCheckCUDAError("kernel call to kernPDIV");CHKERRQ(ierr); 
  #endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "kernPDIV"
__global__ void  kernPDIV(double* devY,double* devX, int* vlen, double* devW){
 /* w <- x./y */
  int tid;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunkY[PDIVTCOUNT];
  __shared__ double chunkX[PDIVTCOUNT];
  __shared__ double chunkW[PDIVTCOUNT];
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    if(chunkX[threadIdx.x]*chunkY[threadIdx.x]!=0){
      chunkW[threadIdx.x]=__ddiv_rn(chunkX[threadIdx.x],chunkY[threadIdx.x]);
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
  PetscErrorCode ierr;
  double *devScratch,*devPartial,zhost;
  PetscInt i,chunks=0,segment,partialsize,scratchsize,*devChunks;
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
  if(chunks>1){
    segment = (int) (CHUNKWIDTH);
    dimGrid.x=ceil((CHUNKWIDTH)/(float)THRNRMCNT);
  }else{
    segment = x->map->n;
    dimGrid.x=ceil(((float)segment)/(float)THRNRMCNT);
  }
  dimBlock.x  = THRNRMCNT;
  if(dimGrid.x&1)dimGrid.x++;
  scratchsize = chunks*dimGrid.x*sizeof(double);
  cms[0] = cudaMalloc((void**)&devScratch,scratchsize);
  ccs[0] = cudaMemsetAsync(devScratch,0,scratchsize,xd->stream);
  partialsize = chunks*sizeof(double);
  cms[1] = cudaMalloc((void**)&devPartial,partialsize);
  ccs[1] = cudaMemsetAsync(devPartial,0,partialsize,xd->stream);
  ccs[2]=cudaMemcpyAsync(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("NORM: chunks: %d, seg: %d, blks: %d, par: %d, scr: %d, chksize: %d, N: %d\n",chunks,segment,dimGrid.x,partialsize/8,scratchsize/8,(int)(CHUNKWIDTH),x->map->n);
    #endif
    ierr = VecCheckCUDAStatus(cms[0],"devScratch alloc in VecNorm_SeqGPU"); CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[0],"devScratch memset in VecNorm_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[1],"devPartial alloc in VecNorm_SeqGPU"); CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[1],"devPartial memset in VecNorm_SeqGPU");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(ccs[2],"on copy segment length H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
  #endif
  cudaDeviceSynchronize();
  for(i=0;i<chunks;i++){/* streaming async kernel calls */
    cms[2]=cudaStreamCreate(&(nrmstream[i]));
    ccs[3]=cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,nrmstream[i]);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(cms[2],"on cudaStreamCreate");CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(ccs[3],"on copy array length H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
    /* Overlapping execution */
    kernNorm2<<<dimGrid,dimBlock,0,nrmstream[i]>>>(xd->devptr,
                                                   xd->segment,
                                                   xd->length,
                                                   xd->offset,
                                                   devScratch+i*dimGrid.x,
                                                   (devPartial+i));
    #if(DEBUGVEC)
        ierr = VecCheckCUDAError("kernNorm2 launch in VecNorm_SeqGPU");CHKERRQ(ierr); 
    #endif
  }
  if(chunks>1){
    /* norm2 block reduction */
    dimGrid.x  = 1;
    dimBlock.x = chunks;
    if(dimBlock.x&1)dimBlock.x++;
    cms[3] = cudaMalloc((void**)&devChunks,sizeof(int));
    ccs[4] = cudaMemcpyAsync(devChunks,&chunks,sizeof(int),cudaMemcpyHostToDevice,xd->stream);
    #if(DEBUGVEC)
        ierr = VecCheckCUDAStatus(cms[3],"on cudamalloc()");CHKERRQ(ierr);
        ierr = VecCheckCUDAStatus(ccs[4],"on copy chunks H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    kernRedNorm<<<dimGrid,dimBlock,chunks*sizeof(double)>>>(devPartial,devChunks,xd->zval);
    #if(DEBUGVEC)
       ierr = VecCheckCUDAError("kernRedNorm_double launch in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif

    /* Copy back norm z */
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    ccs[5]=cudaMemcpy(&zhost,xd->zval,sizeof(double),cudaMemcpyDeviceToHost);/* copy back */
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    cms[4]=cudaFree(devChunks);
    #if(DEBUGVEC)
      ierr = VecCheckCUDAStatus(ccs[5],"on zval copy D2H");CHKERRQ(ierr);
      ierr = VecCheckCUDAStatus(cms[4],"on cudaFree()");CHKERRQ(ierr);
    #endif
  }else{
    cudaDeviceSynchronize();/* make sure everyone is caught up */
    ccs[6]=cudaMemcpy(&zhost,devPartial,sizeof(double),cudaMemcpyDeviceToHost);/* copy back z */
    #if(DEBUGVEC)
       ierr = VecCheckCUDAStatus(ccs[6],"on copy znorm devP D2H in VecNorm_SeqGPU");CHKERRQ(ierr);
    #endif
  }
  *z = PetscSqrtScalar(zhost);

  /* clean up resources */
  for(i=0;i<chunks;i++){
    cudaStreamDestroy(nrmstream[i]);
  }
  free(nrmstream);
  cms[3] = cudaFree(devPartial);
  cms[4] = cudaFree(devScratch);
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Znorm: %e, zhost: %e\n",*z,zhost);
    #endif
    ierr = VecCheckCUDAStatus(cms[3],"on cudaFree(devPartial)");CHKERRQ(ierr);
    ierr = VecCheckCUDAStatus(cms[4],"on cudaFree(devScratch)");CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}



extern __shared__ double zptr[];
#undef __FUNCT__
#define __FUNCT__ "kernRedNorm"
__global__ void kernRedNorm(double* arr, int* chunks,double* z){/* reduction kernel */
  int i = (blockDim.x+1)/2;
  zptr[threadIdx.x]=0.;
  if(threadIdx.x<*chunks)zptr[threadIdx.x]=arr[threadIdx.x];

  #if(DEBUGVEC && VERBOSE)
     printf("In kernRedNorm z[%d]: %e, arr[%d]: %e\n",threadIdx.x,zptr[threadIdx.x],threadIdx.x,arr[threadIdx.x]);
  #endif
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i){
      #if(DEBUGVEC && VERBOSE)
        printf("in kernRedNorm PRE z[%d]: %e, z[%d]: %e\n",threadIdx.x,zptr[threadIdx.x],threadIdx.x+i,zptr[threadIdx.x+i]);
      #endif
      zptr[threadIdx.x]+=zptr[threadIdx.x+i];
      #if(DEBUGVEC && VERBOSE)
        printf("in kernRedNorm POST z[%d]: %e, z[%d]: %e\n",threadIdx.x,zptr[threadIdx.x],threadIdx.x+i,zptr[threadIdx.x+i]);
      #endif
    }
    __syncthreads();
    i/=2;
  }/* end while */
  if(threadIdx.x==0){
    *z=zptr[0];
  }
  #if(DEBUGVEC && VERBOSE)
     if(threadIdx.x==0)printf("in kernRedNorm: z: %e, zptr[0]: %e\n",*z,zptr[0]); 
  #endif
}


#undef __FUNCT__
#define __FUNCT__ "kernNorm2"
__global__ void kernNorm2(double* devX,
                                 int* segmentsize,int* arrsize,
                                 int* offset,double *scratch, double *z){

  __shared__ double chunkX[THRNRMCNT];
  __shared__ int n;    n   = *arrsize;
  __shared__ int seg;  seg = *segmentsize;
  __shared__ int off;  off = *offset;

  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  int i = (blockDim.x+1)/2;
  int j = (gridDim.x+1)/2;
  int item = seg*off+tid;

  if(item<n){/* read in values to shared */
    chunkX[threadIdx.x]=devX[item]; /* offset values */
    #if(DEBUGVEC && VERBOSE)
       printf("in kernNorm2: chunkX[%d]: %e\n",threadIdx.x,chunkX[threadIdx.x]);
    #endif
  }else{
    chunkX[threadIdx.x]=0.0;
  }
  chunkX[threadIdx.x]*=chunkX[threadIdx.x];
  #if(DEBUGVEC && VERBOSE)
     printf("in kernNorm2: chunkX[%d]: %e\n",threadIdx.x,chunkX[threadIdx.x]);
  #endif
  __syncthreads();

  /* block level reduction */
  while(i>0){
     if(threadIdx.x<i){
       chunkX[threadIdx.x]+=chunkX[threadIdx.x+i];
     }
     __syncthreads();
     i/=2;
  }/* end while */
  __syncthreads();
  if(threadIdx.x==0){
    scratch[blockIdx.x]=chunkX[0];
    #if(DEBUGVEC && VERBOSE)
       printf("in kernNorm2: scratch[%d]: %e, chunk[0]: %e\n",blockIdx.x,scratch[blockIdx.x],chunkX[0]);
    #endif
  }
  __syncthreads();

  /* grid level reduction */
  while(j>0){
    if(threadIdx.x==0 && blockIdx.x<j){
      scratch[blockIdx.x]+=scratch[blockIdx.x+j];
      #if(DEBUGVEC && VERBOSE)
         printf("############## scratch[%d]: %e\n",blockIdx.x,scratch[blockIdx.x]);
      #endif
    }
    __syncthreads();
    j/=2;
  }/* end while */
  if(tid==0){
    *z=scratch[blockIdx.x];
    __threadfence();
    #if(DEBUGVEC && VERBOSE)
      printf(">>>>>>>>>>>>>> scratch[%d]: %e, z: %e\n",blockIdx.x,scratch[blockIdx.x],*z);
    #endif
    //correct semantics up to here
  }
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
  if(vd->syncState==VEC_GPU || flg1 || flg2){
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
  if(vd->syncState==VEC_CPU || vd->syncState==VEC_ALLOC){
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
  #if(DEBUGVEC)
    #if(VERBOSE)
       printf("Call to VecCopy_SeqGPU\n");
    #endif

    PetscBool same=PETSC_FALSE;
    cudaDeviceSynchronize();
    ierr = VecCompare_SeqGPU(s,d,&same,0,0);CHKERRQ(ierr);
    if(!same)SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Vector duplication failed.");
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
  /* using pinned memory */
  ierr = PinnedMalloc(&(seqgpu->cpuptr),V->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemzero(seqgpu->cpuptr,V->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
  seqgpu->syncState=VEC_ALLOC;
  //ierr = PetscNewLog(V,Vec_SeqGPU,&(V->data));CHKERRQ(ierr);

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
  #endif
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "VecDestroy_SeqGPU"
PetscErrorCode VecDestroy_SeqGPU(Vec v){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  if(vd && vd->syncState != VEC_UNALLOC){
      cms[0]=cudaFree(vd->devptr);  vd->devptr=PETSC_NULL;
      cms[1]=cudaFree(vd->length);  vd->length=PETSC_NULL;
      cms[2]=cudaFree(vd->segment); vd->segment=PETSC_NULL;
      cms[3]=cudaFree(vd->zval);    vd->zval=PETSC_NULL;
      cms[4] = cudaStreamDestroy(vd->stream);
      ierr = PinnedFree(vd->cpuptr); CHKERRQ(ierr);
      #if(DEBUGVEC)
        ierr=VecCheckCUDAStatus(cms[0],"destroying devptr in VecDestroy_SeqGPU"); CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[1],"destroying length in VecDestroy_SeqGPU"); CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[2],"destroying segment in VecDestroy_SeqGPU");CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[3],"destroying zval in VecDestroy_SeqGPU");   CHKERRQ(ierr);
        ierr=VecCheckCUDAStatus(cms[4],"destroying stream in VecDestroy_SeqGPU"); CHKERRQ(ierr);
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
  PetscErrorCode    ierr;
  PetscInt          i,n = xin->map->n;
  const char        *name;
  PetscViewerFormat format;
  PetscScalar *xv;

  PetscFunctionBegin;
  ierr = VecGetArray_SeqGPU(xin,&xv);CHKERRQ(ierr);
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_MATLAB) {
    ierr = PetscObjectGetName((PetscObject)xin,&name);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%s = [\n",name);CHKERRQ(ierr);
    for (i=0; i<n; i++) {
#if defined(PETSC_USE_COMPLEX)
      if (PetscImaginaryPart(xv[i]) > 0.0) {
        ierr = PetscViewerASCIIPrintf(viewer,"%18.16e + %18.16ei\n",PetscRealPart(xv[i]),PetscImaginaryPart(xv[i]));CHKERRQ(ierr);
      } else if (PetscImaginaryPart(xv[i]) < 0.0) {
        ierr = PetscViewerASCIIPrintf(viewer,"%18.16e - %18.16ei\n",PetscRealPart(xv[i]),-PetscImaginaryPart(xv[i]));CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",PetscRealPart(xv[i]));CHKERRQ(ierr);
      }
#else
      ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double) xv[i]);CHKERRQ(ierr);
#endif
    }
    ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);
  } else if (format == PETSC_VIEWER_ASCII_SYMMODU) {
    for (i=0; i<n; i++) {
#if defined(PETSC_USE_COMPLEX)
      ierr = PetscViewerASCIIPrintf(viewer,"%18.16e %18.16e\n",PetscRealPart(xv[i]),PetscImaginaryPart(xv[i]));CHKERRQ(ierr);
#else
      ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",xv[i]);CHKERRQ(ierr);
#endif
    }
  } else if (format == PETSC_VIEWER_ASCII_VTK || format == PETSC_VIEWER_ASCII_VTK_CELL) {
    /* 
       state 0: No header has been output
       state 1: Only POINT_DATA has been output
       state 2: Only CELL_DATA has been output
       state 3: Output both, POINT_DATA last
       state 4: Output both, CELL_DATA last 
    */
    static PetscInt stateId = -1;
    int outputState = 0;
    PetscBool  hasState;
    int doOutput = 0;
    PetscInt bs, b;

    if (stateId < 0) {
      ierr = PetscObjectComposedDataRegister(&stateId);CHKERRQ(ierr);
    }
    ierr = PetscObjectComposedDataGetInt((PetscObject) viewer, stateId, outputState, hasState);CHKERRQ(ierr);
    if (!hasState) {
      outputState = 0;
    }
    ierr = PetscObjectGetName((PetscObject) xin, &name);CHKERRQ(ierr);
    ierr = VecGetBlockSize(xin, &bs);CHKERRQ(ierr);
    if ((bs < 1) || (bs > 3)) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "VTK can only handle 3D objects, but vector dimension is %d", bs);
    }
    if (format == PETSC_VIEWER_ASCII_VTK) {
      if (outputState == 0) {
        outputState = 1;
        doOutput = 1;
      } else if (outputState == 1) {
        doOutput = 0;
      } else if (outputState == 2) {
        outputState = 3;
        doOutput = 1;
      } else if (outputState == 3) {
        doOutput = 0;
      } else if (outputState == 4) {
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "Tried to output POINT_DATA again after intervening CELL_DATA");
      }
      if (doOutput) {
        ierr = PetscViewerASCIIPrintf(viewer, "POINT_DATA %d\n", n/bs);CHKERRQ(ierr);
      }
    } else {
      if (outputState == 0) {
        outputState = 2;
        doOutput = 1;
      } else if (outputState == 1) {
        outputState = 4;
        doOutput = 1;
      } else if (outputState == 2) {
        doOutput = 0;
      } else if (outputState == 3) {
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "Tried to output CELL_DATA again after intervening POINT_DATA");
      } else if (outputState == 4) {
        doOutput = 0;
      }
      if (doOutput) {
        ierr = PetscViewerASCIIPrintf(viewer, "CELL_DATA %d\n", n);CHKERRQ(ierr);
      }
    }
    ierr = PetscObjectComposedDataSetInt((PetscObject) viewer, stateId, outputState);CHKERRQ(ierr);
    if (name) {
      if (bs == 3) {
        ierr = PetscViewerASCIIPrintf(viewer, "VECTORS %s double\n", name);CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer, "SCALARS %s double %d\n", name, bs);CHKERRQ(ierr);
      }
    } else {
      ierr = PetscViewerASCIIPrintf(viewer, "SCALARS scalars double %d\n", bs);CHKERRQ(ierr);
    }
    if (bs != 3) {
      ierr = PetscViewerASCIIPrintf(viewer, "LOOKUP_TABLE default\n");CHKERRQ(ierr);
    }
    for (i=0; i<n/bs; i++) {
      for (b=0; b<bs; b++) {
        if (b > 0) {
          ierr = PetscViewerASCIIPrintf(viewer," ");CHKERRQ(ierr);
        }
#if !defined(PETSC_USE_COMPLEX)
        ierr = PetscViewerASCIIPrintf(viewer,"%G",xv[i*bs+b]);CHKERRQ(ierr);
#endif
      }
      ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    }
  } else if (format == PETSC_VIEWER_ASCII_VTK_COORDS) {
    PetscInt bs, b;

    ierr = VecGetBlockSize(xin, &bs);CHKERRQ(ierr);
    if ((bs < 1) || (bs > 3)) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "VTK can only handle 3D objects, but vector dimension is %d", bs);
    }
    for (i=0; i<n/bs; i++) {
      for (b=0; b<bs; b++) {
        if (b > 0) {
          ierr = PetscViewerASCIIPrintf(viewer," ");CHKERRQ(ierr);
        }
#if !defined(PETSC_USE_COMPLEX)
        ierr = PetscViewerASCIIPrintf(viewer,"%G",xv[i*bs+b]);CHKERRQ(ierr);
#endif
      }
      for (b=bs; b<3; b++) {
        ierr = PetscViewerASCIIPrintf(viewer," 0.0");CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    }
  } else if (format == PETSC_VIEWER_ASCII_PCICE) {
    PetscInt bs, b;

    ierr = VecGetBlockSize(xin, &bs);CHKERRQ(ierr);
    if ((bs < 1) || (bs > 3)) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "PCICE can only handle up to 3D objects, but vector dimension is %d", bs);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"%D\n", xin->map->N/bs);CHKERRQ(ierr);
    for (i=0; i<n/bs; i++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%7D   ", i+1);CHKERRQ(ierr);
      for (b=0; b<bs; b++) {
        if (b > 0) {
          ierr = PetscViewerASCIIPrintf(viewer," ");CHKERRQ(ierr);
        }
#if !defined(PETSC_USE_COMPLEX)
        ierr = PetscViewerASCIIPrintf(viewer,"% 12.5E",xv[i*bs+b]);CHKERRQ(ierr);
#endif
      }
      ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    }
  } else {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)xin,viewer,"Vector Object");CHKERRQ(ierr);
    for (i=0; i<n; i++) {
      if (format == PETSC_VIEWER_ASCII_INDEX) {
        ierr = PetscViewerASCIIPrintf(viewer,"%D: ",i);CHKERRQ(ierr);
      }
#if defined(PETSC_USE_COMPLEX)
      if (PetscImaginaryPart(xv[i]) > 0.0) {
        ierr = PetscViewerASCIIPrintf(viewer,"%G + %G i\n",PetscRealPart(xv[i]),PetscImaginaryPart(xv[i]));CHKERRQ(ierr);
      } else if (PetscImaginaryPart(xv[i]) < 0.0) {
        ierr = PetscViewerASCIIPrintf(viewer,"%G - %G i\n",PetscRealPart(xv[i]),-PetscImaginaryPart(xv[i]));CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer,"%G\n",PetscRealPart(xv[i]));CHKERRQ(ierr);
      }
#else
      ierr = PetscViewerASCIIPrintf(viewer,"%G\n",(double) xv[i]);CHKERRQ(ierr);
#endif
    }
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = VecRestoreArray_SeqGPU(xin,&xv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}




















EXTERN_C_END
