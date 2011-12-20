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
__constant__ int     devN;//vector length
__constant__ double  dblScalarValue;//utility var
__constant__ double2 dblScalar2Value;//utility var
__constant__ float   fltScalarValue;//utility var
__constant__ float2  fltScalar2Value;//utility var


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
  cudaError_t cudastatus;
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

      cudastatus = cudaMalloc((void**)&devseeds,bx*sizeof(uint));
      ierr = VecCheckCUDAStatus(cudastatus,"error in cudaMalloc");CHKERRQ(ierr);

      cudastatus=cudaMemcpy(devseeds,seeds,bx*sizeof(uint),cudaMemcpyHostToDevice);
      ierr = VecCheckCUDAStatus(cudastatus,"on copy H2D in VecSetRandom_SeqGPU");CHKERRQ(ierr);
      xd->vstat.h2d_count++;
      xd->vstat.h2d_bytes+=bx*sizeof(uint);

      kernRandS<<<dimGrid,dimBlock>>>(devseeds);
      ierr = VecCheckCUDAError("kernRandS launch");CHKERRQ(ierr);
      ierr = PetscFree(seeds);CHKERRQ(ierr);
      cudaDeviceSynchronize();
      cudastatus = cudaFree(devseeds);
      ierr = VecCheckCUDAStatus(cudastatus,"in cudaFree()");CHKERRQ(ierr);
      seed_flag=PETSC_FALSE;
    }
    kernRand<<<dimGrid,dimBlock>>>(xd->devptr,xd->length);
    ierr = VecCheckCUDAError("kernRand launch");CHKERRQ(ierr);
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
  /* printf("Call to VecGetLocalSize_SeqGPU\n"); */
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
  /* printf("Call to VecGetSize_SeqGPU\n"); */
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
  //printf("Call to VecCopyOverDevice\n");
  cudaError_t cudastatus;
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
  cudastatus=cudaMemcpyAsync(dd->devptr,sd->devptr,
               s->map->n*sizeof(PetscScalar),cudaMemcpyDeviceToDevice,dd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2D in VecCopyOverDevice");CHKERRQ(ierr);
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
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  /* printf("Call to VecCopyBlockH2D\n"); */
  cudastatus=cudaMemcpyAsync(&(vd->devptr[offset]),y,
               blocksize*sizeof(PetscScalar),cudaMemcpyHostToDevice,vd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy H2D in VecCopyBlockH2D");CHKERRQ(ierr);
  vd->vstat.h2d_count++;
  vd->vstat.h2d_bytes+=blocksize*sizeof(PetscScalar);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCopyOverH2D"
PetscErrorCode VecCopyOverH2D(Vec v,PetscScalar *y){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  /* printf("Call to VecCopyOverH2D\n"); */
  cudastatus=cudaMemcpyAsync(vd->devptr,y,
               v->map->n*sizeof(PetscScalar),cudaMemcpyHostToDevice,vd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy H2D in VecCopyOverH2D");CHKERRQ(ierr);
  vd->vstat.h2d_count++;
  vd->vstat.h2d_bytes+=v->map->n*sizeof(PetscScalar);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCopyBlockD2H"
PetscErrorCode VecCopyBlockD2H(Vec v,PetscScalar *y,PetscInt offset, PetscInt blocksize){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  /* printf("Call to VecCopyBlockD2H\n"); */
  cudastatus=cudaMemcpyAsync(y,&(vd->devptr[offset]),
               blocksize*sizeof(PetscScalar),cudaMemcpyDeviceToHost,vd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H in VecCopyBlockD2H");CHKERRQ(ierr);
  vd->vstat.d2h_count++;
  vd->vstat.d2h_bytes+=blocksize*sizeof(PetscScalar);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCopyOverD2H"
PetscErrorCode VecCopyOverD2H(Vec v,PetscScalar *y){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  /* printf("Call to VecCopyOverD2H\n"); */
  cudastatus=cudaMemcpyAsync(y,vd->devptr,
               v->map->n*sizeof(PetscScalar),cudaMemcpyDeviceToHost,vd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H in VecCopyOverD2H");CHKERRQ(ierr);
  vd->vstat.d2h_count++;
  vd->vstat.d2h_bytes+=v->map->n*sizeof(PetscScalar);
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
  //printf("Call to VecSetValues_SeqGPU\n");
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
  cudaDeviceSynchronize();
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VecSet_SeqGPU"
PetscErrorCode VecSet_SeqGPU(Vec xin,PetscScalar alpha){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  dim3 dimgrid(ceil((float)xin->map->n/((float)TCOUNT)),1,1);
  dim3 dimblocks(TCOUNT,1,1);
  Vec_SeqGPU* xd = (Vec_SeqGPU*)xin->data;
  //printf("Call to VecSet_SeqGPU, alpha: %e\n",alpha);
  cudaDeviceSynchronize();
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,
            PETSC_ERR_MEM,"*** In VecSet_SeqGPU, Vec not allocated. ***\n");
  }else{
    cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
    kernSet<<<dimgrid,dimblocks>>>(xd->devptr,xd->length);
    ierr = VecCheckCUDAError("Call to kernSet. "); CHKERRQ(ierr);
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
    //printf("in kernSet: x[%d]: %e\n",tid,x[tid]);
  }
}




#undef __FUNCT__
#define __FUNCT__ "VecScale_SeqGPU"
PetscErrorCode VecScale_SeqGPU(Vec x, PetscScalar alpha){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  dim3 dimgrid(ceil((float)x->map->n/((float)TCOUNT)),1,1);
  dim3 dimblocks(TCOUNT,1,1);
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  //printf("VecScale_SeqGPU...alpha: %e\n",alpha);
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
    cudastatus = cudaMemsetAsync(xd->devptr,0,x->map->n*sizeof(double),xd->stream);
    ierr = VecCheckCUDAStatus(cudastatus,"error in cudaMemset");CHKERRQ(ierr);
  }else if (alpha != 1.0){
    cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
    kernScale<<<dimgrid,dimblocks,0,xd->stream>>>(xd->devptr,xd->length);
    ierr = VecCheckCUDAError("Call to kernScale. "); CHKERRQ(ierr);
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
    //if(x[tid]!=0)printf("Pre: kernScale: x[%d]: %e, alpha: %e\n",tid,x[tid],localdbl);
    arr[threadIdx.x] *= localdbl;
    x[tid] = arr[threadIdx.x];
    //if(x[tid]!=0)printf("kernScale: x[%d]: %e, alpha: %e\n",tid,x[tid],localdbl);
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
  cudaError_t cudastatus;
  double *devScratch,*devPartial,zhost;
  PetscInt i,chunks=0,segment,partialsize,scratchsize;
  cudaStream_t* dotstream;
  dim3 dimGrid, dimBlock;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;

  //printf("Call to VecDot_SeqGPU, chunkwidth: %f xlen: %d\n",(CHUNKWIDTH),x->map->n);

  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) /(float)(CHUNKWIDTH));
  //printf("Number of chunks in Dot: %d\n",chunks);
  dotstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));

  if(chunks>1){
    segment = (int) CHUNKWIDTH;
    dimGrid.x=ceil((CHUNKWIDTH)/(float)THRDOTCNT);
  }else{
    segment = x->map->n;
    dimGrid.x=ceil(((float)segment)/(float)THRDOTCNT);
  }
  dimBlock.x = THRDOTCNT;
  partialsize=chunks*sizeof(double);
  scratchsize=partialsize*dimGrid.x;
  /* set up on x stream */
  if(xd->syncState==VEC_CPU){
    printf("xd state VEC_CPU: copying to device.\n");
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  cudastatus = cudaMalloc((void**)&devScratch,scratchsize);/* scratch pad */
  ierr = VecCheckCUDAStatus(cudastatus,"devScratch alloc in VecDot_SeqGPU");CHKERRQ(ierr);
  cudastatus = cudaMemsetAsync(devScratch,0,scratchsize,xd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"devScratch memset in VecDot_SeqGPU");CHKERRQ(ierr);

  /* set up on y stream */
  if(yd->syncState==VEC_CPU){
    printf("yd state VEC_CPU: copying to device.\n");
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }
  cudastatus=cudaMalloc((void**)&devPartial,partialsize);/* partial results to be combined */
  ierr = VecCheckCUDAStatus(cudastatus,"devPartial alloc in VecDot_SeqGPU");CHKERRQ(ierr);
  cudastatus=cudaMemsetAsync(devPartial,0,partialsize,yd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"devPartial memset in VecDot_SeqGPU");CHKERRQ(ierr);

  cudaDeviceSynchronize();/* make sure everyone is ready */

  for(i=0;i<chunks;i++){  /* streaming async kernel calls */
    cudastatus=cudaStreamCreate(&(dotstream[i]));
    ierr = VecCheckCUDAStatus(cudastatus,"on cudaStreamCreate");CHKERRQ(ierr);
    cudastatus=cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,dotstream[i]);
    ierr = VecCheckCUDAStatus(cudastatus,"on copy array length H2D in VecDot_SeqGPU");CHKERRQ(ierr);
    cudastatus=cudaMemcpyAsync(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice,dotstream[i]);
    ierr = VecCheckCUDAStatus(cudastatus,"on copy segment size H2D in VecDot_SeqGPU");CHKERRQ(ierr);
    /* Overlapping execution */
    kernDot<<<dimGrid,dimBlock,0,dotstream[i]>>>(xd->devptr,yd->devptr,
                                                          xd->segment,
                                                          xd->length,
                                                          xd->offset,
                                                          (devScratch+i*dimGrid.x),
                                                 (devPartial+i));
    ierr = VecCheckCUDAError("kernDot launch in VecDot_SeqGPU");CHKERRQ(ierr);
    xd->vstat.h2d_count++;
    xd->vstat.h2d_bytes+=2*sizeof(int);
    yd->vstat.h2d_count++;
    yd->vstat.h2d_bytes+=2*sizeof(int);
  }

  /* dot product block reduction */
  dimGrid.x  = 1;
  dimBlock.x = chunks;
  cudaDeviceSynchronize();/* make sure everyone is caught up */
  kernRedDot<<<dimGrid,dimBlock,chunks*sizeof(double),xd->stream>>>(devPartial,xd->zval);
  ierr = VecCheckCUDAError("kernRedDot launch in VecDot_SeqGPU");CHKERRQ(ierr);
  cudaDeviceSynchronize();/* make sure everyone is caught up */

  /* Copy back dot z */
  cudastatus=cudaMemcpy(&zhost,xd->zval,sizeof(double),cudaMemcpyDeviceToHost);/* copy back z */
  ierr = VecCheckCUDAStatus(cudastatus,"on copy zdot D2H in VecDot_SeqGPU");CHKERRQ(ierr);
  *z=zhost;
  //printf("Zdot: %e\n",*z);

  /* clean up resources */
  cudastatus = cudaFree(devPartial);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  for(i=0;i<chunks;i++){
    cudastatus = cudaStreamDestroy(dotstream[i]);
    ierr = VecCheckCUDAStatus(cudastatus,"on cudaStreamCreate");CHKERRQ(ierr);
  }
  free(dotstream);
  cudastatus = cudaFree(devScratch);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


extern __shared__ double arrayDot[];
#undef __FUNCT__
#define __FUNCT__ "kernRedDot"
__global__ void kernRedDot(double* arr,double* z){/* reduction kernel */

  int i = (blockDim.x+1)/2;
  double* zptr=(double*)arrayDot;

  zptr[threadIdx.x]=arr[threadIdx.x];
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i){
      zptr[threadIdx.x]+=zptr[threadIdx.x+i];
    }
    __syncthreads();
    i/=2;
  }
  if(threadIdx.x==0){
    *z=zptr[0];
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
  // printf("VecMDot_SeqGPU\n");
  //printf("Number of vectors in MDot: %d\n",nv);
  for (i=0; i<nv; i++) {
    ierr = VecDot_SeqGPU(x,y[i],&val[i]);CHKERRQ(ierr);
    //cudaDeviceSynchronize();
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
  cudaError_t cudastatus;
  //printf("VecWAXPY_SeqGPU...");
  //printf("alpha: %e\n",alpha);
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
    ierr = VecCheckCUDAError("kernel call to kernWXPY");CHKERRQ(ierr);
  }else if(alpha==-1.0){
    kernWXMY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    ierr = VecCheckCUDAError("kernel call to kernWXMY");CHKERRQ(ierr);
  }else{
    cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
    kernWAXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
    ierr = VecCheckCUDAError("kernel call to kernWAXPY");CHKERRQ(ierr);
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
  /* printf("in kernWAXPY:alphaShared: %f, tid: %d, vlen: %d\n",alphaShared,tid,*vlen); */
  if(tid<*vlen){
    //if(devX[tid]!=0)printf("kernWAXPY: devX[%d]: %e\n",tid,devX[tid]);
    //if(devY[tid]!=0)printf("kernWAXPY: devY[%d]: %e\n",tid,devY[tid]);
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]+(chunkX[threadIdx.x]*alphaShared);
    devW[tid]=chunkW[threadIdx.x];
    //if(devW[tid]!=0)printf("kernWAXPY: devW[%d]: %e, alpha: %e\n",tid,devW[tid],alphaShared);
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

  /* printf("in kernWAXPY:alphaShared: %f, tid: %d, vlen: %d\n",alphaShared,tid,*vlen); */
  if(tid<*vlen){
    //if(devX[tid]!=0)printf("kernWXPY: devX[%d]: %e\n",tid,devX[tid]);
    //if(devY[tid]!=0)printf("kernWXPY: devY[%d]: %e\n",tid,devY[tid]);
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]+chunkX[threadIdx.x];
    devW[tid]=chunkW[threadIdx.x];
    //if(devW[tid]!=0)printf("kernWXPY: devW[%d]: %e\n",tid,devW[tid]);
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
  /* printf("in kernWAXPY:alphaShared: %f, tid: %d, vlen: %d\n",alphaShared,tid,*vlen); */
  if(tid<*vlen){
    //if(devX[tid]!=0)printf("kernWXMY: devX[%d]: %e\n",tid,devX[tid]);
    //if(devY[tid]!=0)printf("kernWXMY: devY[%d]: %e\n",tid,devY[tid]);
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]-chunkX[threadIdx.x];
    devW[tid]=chunkW[threadIdx.x];
    //if(devW[tid]!=0)printf("kernWXMY: devW[%d]: %e\n",tid,devW[tid]);
  }
}

#undef __FUNCT__
#define __FUNCT__ "VecMAXPY_SeqGPU"
PetscErrorCode VecMAXPY_SeqGPU(Vec x,PetscInt nv,const PetscScalar* alpha,Vec *y){
  /* y = y + sum(a[i]*x[i]) */
  PetscFunctionBegin;
  //printf("VecMAXPY_SeqGPU: alpha: %e\n",*alpha);
  PetscErrorCode ierr;
  PetscInt i;  PetscInt bx,tx;
  dim3 dimGrid;
  dim3 dimBlock;
  cudaError_t cudastatus;
  PetscScalar *devW;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=PETSC_NULL;

  cudastatus = cudaMalloc((void**)&devW,x->map->n*sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"error in device malloc");CHKERRQ(ierr);
  cudastatus = cudaMemset(devW,0,x->map->n*sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"error in device memset");CHKERRQ(ierr);

  /* assuming xwidth mem load isn't going to be an issue */
  bx=ceil((float)x->map->n/(float)AXPYTCOUNT);
  tx=AXPYTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;

  //printf("Number of vectors in MAXPY: %d\n",nv);
  //ierr = VecCheck_SeqGPU(x);CHKERRQ(ierr);
  for(i=0;i<nv;i++){
     if(y[i]->map->n!=x->map->n){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
    }
    yd=(Vec_SeqGPU*)y[i]->data;
    if(yd->syncState==VEC_CPU){/* synch x */
      ierr = VecCopyOverH2D(y[i],yd->cpuptr);CHKERRQ(ierr);
      yd->syncState=VEC_SYNCHED;
    }
    cudaDeviceSynchronize();
    cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha[i],sizeof(double),0,cudaMemcpyHostToDevice);
    ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
    //printf("Alpha[%d]: %e\n", i, alpha[i]);
    if(alpha[i]==0){
      continue;
    }else if(alpha[i]==1.){
      /* assuming width mem load isn't going to be an issue */
      kernXPY<<<dimGrid,dimBlock>>>(devW,yd->devptr,yd->length);
      ierr = VecCheckCUDAError("kernel call to kernXPY");CHKERRQ(ierr);
    }else{
      /* assuming width mem load isn't going to be an issue */
      kernAXPY<<<dimGrid,dimBlock>>>(devW,yd->devptr,yd->length);
      ierr = VecCheckCUDAError("kernel call to kernAXPY");CHKERRQ(ierr);
    }
  }
  if(xd->syncState==VEC_CPU){/* synch x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  cudaDeviceSynchronize();
  kernXPY<<<dimGrid,dimBlock>>>(xd->devptr,devW,xd->length);
  ierr = VecCheckCUDAError("kernel call to kernXPY");CHKERRQ(ierr);

  cudastatus = cudaFree(devW);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree");CHKERRQ(ierr);
  xd->syncState=VEC_GPU;
  //ierr = VecCheck_SeqGPU(x);CHKERRQ(ierr);
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
  cudaError_t cudastatus;
  //printf("VecAXPY_SeqGPU\n");

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
  cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
  ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
  cudaDeviceSynchronize();
  if(alpha==1.){
    /* assuming width mem load isn't going to be an issue */
    bx=ceil((float)x->map->n/(float)AXPYTCOUNT);
    tx=AXPYTCOUNT;
    dimGrid.x=bx; dimGrid.y=1;
    dimBlock.x=tx; dimBlock.y=1;
    kernXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,yd->length);
    ierr = VecCheckCUDAError("kernel call to kernXPY");CHKERRQ(ierr);
  }else if(alpha!=0){
    /* assuming width mem load isn't going to be an issue */
    bx=ceil((float)x->map->n/(float)AXPYTCOUNT);
    tx=AXPYTCOUNT;
    dimGrid.x=bx; dimGrid.y=1;
    dimBlock.x=tx; dimBlock.y=1;
    kernAXPY<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,yd->length);
    ierr = VecCheckCUDAError("kernel call to kernAXPY");CHKERRQ(ierr);
  }
  //cudaDeviceSynchronize();
  yd->syncState=VEC_GPU;
  //ierr = VecCheck_SeqGPU(y);CHKERRQ(ierr);
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
  PetscErrorCode ierr;
  int blocks=ceil((float)x->map->n/(float)AXPBYPCZTCOUNT);/* assuming shared memory size is not an issue */
  int threads=AXPBYPCZTCOUNT;
  cudaError_t cudastatus;
  Vec_SeqGPU* devX = (Vec_SeqGPU*)x->data;
  Vec_SeqGPU* devY = (Vec_SeqGPU*)y->data;
  Vec_SeqGPU* devZ = (Vec_SeqGPU*)z->data;

  double2 alphabeta;
  alphabeta.x = alpha;
  alphabeta.y = beta;
  cudaDeviceSynchronize();
  cudastatus=cudaMemcpyToSymbol("dblScalar2Value",(void*)&alphabeta,sizeof(double2),0,cudaMemcpyHostToDevice);
  ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
  cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&gamma,sizeof(double),0,cudaMemcpyHostToDevice);
  ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);

  dim3 dimGrid; dimGrid.x=blocks; dimGrid.y=1;
  dim3 dimBlock; dimBlock.x=threads; dimBlock.y=1;
  kernAXPBYPCZ<<<dimGrid,dimBlock>>>(devX->devptr,devY->devptr,devZ->devptr,devX->length);
  ierr = VecCheckCUDAError("launch kernAXPBYPCZ");CHKERRQ(ierr);
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
  printf("VecPointwiseMult_SeqGPU\n");
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
  cudaDeviceSynchronize();
  /* assuming width mem load isn't going to be an issue */
  bx=ceil((float)y->map->n/(float)PMULTCOUNT);
  tx=PMULTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;
  kernPMULT<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  ierr = VecCheckCUDAError("kernel call to kernPMULT");CHKERRQ(ierr);
  //cudaDeviceSynchronize();
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
  printf("VecMaxPointwiseDivide_SeqGPU...");
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  PetscInt i,bx,tx;
  PetscScalar *maxlist=PETSC_NULL;
  PetscScalar *devmaxlist=PETSC_NULL;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid;
  dim3 dimBlock;
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
  cudaDeviceSynchronize();
  /* assuming width mem load isn't going to be an issue */
  bx=ceil((float)y->map->n/(float)PDIVTCOUNT);
  tx=PDIVTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;

  ierr = PetscMalloc(bx*sizeof(PetscScalar),&maxlist);CHKERRQ(ierr);

  cudastatus=cudaMalloc((void**)&devmaxlist,bx*sizeof(PetscScalar));
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H");CHKERRQ(ierr);

  kernMAXPDIV<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,devmaxlist);
  ierr = VecCheckCUDAError("kernel call to kernPDIV");CHKERRQ(ierr);

  cudastatus=cudaMemcpy(maxlist,devmaxlist,bx*sizeof(PetscScalar),cudaMemcpyDeviceToHost);/* copy back */
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H");CHKERRQ(ierr);

  *max = maxlist[0];
  if(bx>1){/* final collapse */
    for(i=1;i<bx;i++){
      if(maxlist[i]>*max){
        *max=maxlist[i];
      }
    }
  }
  printf("max: %f\n",*max);
  ierr = PetscFree(maxlist);CHKERRQ(ierr);
  cudastatus = cudaFree(devmaxlist);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

extern __shared__ double maxshared[]
#undef __FUNCT__
#define __FUNCT__ "kernMAX"
__global__ void  kernMAX(int* vlen, double* maxlist,double* max){
  int i,tid;
  tid = threadIdx.x;
  i = (blockDim.x+1)/2;
  __shared__  double* slist;
  slist = maxshared;

  slist[tid]=maxlist[tid];
  while(i<0){
    if(tid>i){
      slist[tid] = (slist[tid]>slist[tid+i])?slist[tid]:slist[tid+i];
    }
    __synchthreads();
    i/=2;
  }
  if(tid==0)*max = slist[0];
}




#undef __FUNCT__
#define __FUNCT__ "kernMAXPDIV"
__global__ void  kernMAXPDIV(double* devY,double* devX, int* vlen, double* maxlist){

 /* w <- max(abs(x./y)) */
  int i,tid;
  i = (PDIVTCOUNT+1)/2;
  tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunkY[PDIVTCOUNT];
  __shared__ double chunkX[PDIVTCOUNT];
  __shared__ double chunkW[PDIVTCOUNT];
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    if(chunkY[threadIdx.x]!=0){
      chunkW[threadIdx.x]=fabs(__ddiv_rn(chunkX[threadIdx.x],chunkY[threadIdx.x]));
    }else{
      chunkW[threadIdx.x]=fabs(chunkX[threadIdx.x]);
    }
  }else{
    chunkW[threadIdx.x]=0.0;
  }
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i && chunkW[threadIdx.x]<chunkW[threadIdx.x+i]){
      chunkW[threadIdx.x]=chunkW[threadIdx.x+i];
    }
    i/=2;
    __syncthreads();
  }
  if(threadIdx.x==0)maxlist[blockIdx.x]=chunkW[0];
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
  printf("Call to VecPointwiseDivide_SeqGPU\n");
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
  bx=ceil((float)y->map->n/(float)PDIVTCOUNT);
  tx=PDIVTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;

  kernPDIV<<<dimGrid,dimBlock>>>(yd->devptr,xd->devptr,xd->length,wd->devptr);
  ierr = VecCheckCUDAError("kernel call to kernPDIV");CHKERRQ(ierr);
  //cudaDeviceSynchronize();
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
  //printf("VecDotNorm2_SeqGPU\n");
  ierr = VecDot(s,t,dp); CHKERRQ(ierr);
  ierr = VecNorm(t,NORM_2,nm); CHKERRQ(ierr);
  //printf("dp: %e, nm: %e\n",*dp,*nm);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecNorm_SeqGPU"
PetscErrorCode VecNorm_SeqGPU(Vec x,NormType type,PetscReal* z){
  /* NormType: NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4 */
  /* dealing with NORM_2 for now... */
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  double *devScratch,*devPartial,zhost;
  PetscInt i,chunks=0,segment,partialsize,scratchsize;
  cudaStream_t* nrmstream;
  dim3 dimGrid, dimBlock;
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;

  /* figure out how many chunks will be needed */
  chunks = ceil( ((float)x->map->n) /(float)(CHUNKWIDTH));
  nrmstream = (cudaStream_t*)malloc(chunks*sizeof(cudaStream_t));

  if(chunks>1){
    segment = (int) CHUNKWIDTH;
    dimGrid.x=ceil((CHUNKWIDTH)/(float)THRNRMCNT);
  }else{
    segment = x->map->n;
    dimGrid.x=ceil(((float)segment)/(float)THRNRMCNT);
  }
  dimBlock.x  = THRNRMCNT;
  partialsize = chunks*sizeof(double);
  scratchsize = partialsize*dimGrid.x;

  printf("chunks: %d, segmentsize: %d, dimGrid.x: %d, partialsize: %d, scratchsize: %d, CHUNKWIDTH: %d\n",
         chunks,xd->segment,dimGrid.x,partialsize,scratchsize,(int)(CHUNKWIDTH));
  cudastatus = cudaMalloc((void**)&devScratch,scratchsize);
  ierr = VecCheckCUDAStatus(cudastatus,"devScratch alloc in VecNorm_SeqGPU");CHKERRQ(ierr);
  cudastatus = cudaMemsetAsync(devScratch,0,scratchsize,xd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"devScratch memset in VecNorm_SeqGPU");CHKERRQ(ierr);

  cudastatus=cudaMalloc((void**)&devPartial,partialsize);
  ierr = VecCheckCUDAStatus(cudastatus,"devPartial alloc in VecNorm_SeqGPU");CHKERRQ(ierr);
  cudastatus = cudaMemsetAsync(devPartial,0,partialsize,xd->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"devPartial memset in VecNorm_SeqGPU");CHKERRQ(ierr);

  if(xd->syncState==VEC_CPU){
    printf("xd state VEC_CPU: copying to device.\n");
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  cudaDeviceSynchronize();/* make sure everyone is ready to go */

  for(i=0;i<chunks;i++){/* streaming async kernel calls */
    cudastatus=cudaStreamCreate(&(nrmstream[i]));
    ierr = VecCheckCUDAStatus(cudastatus,"on cudaStreamCreate");CHKERRQ(ierr);
    cudastatus=cudaMemcpyAsync(xd->offset,&i,sizeof(int),cudaMemcpyHostToDevice,nrmstream[i]);
    ierr = VecCheckCUDAStatus(cudastatus,"on copy array length H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
    cudastatus=cudaMemcpyAsync(xd->segment,&segment,sizeof(int),cudaMemcpyHostToDevice,nrmstream[i]);
    ierr = VecCheckCUDAStatus(cudastatus,"on copy segment size H2D in VecNorm_SeqGPU");CHKERRQ(ierr);
    /* Overlapping execution */
    kernNorm2_double<<<dimGrid,dimBlock,0,nrmstream[i]>>>(xd->devptr,
                                                          xd->segment,
                                                          xd->length,
                                                          xd->offset,
                                                          (devScratch+i*dimGrid.x),
                                                          (devPartial+i));
    ierr = VecCheckCUDAError("kernNorm2 launch in VecNorm_SeqGPU");CHKERRQ(ierr);

    xd->vstat.h2d_count++;
    xd->vstat.h2d_bytes+=2*sizeof(int);
  }

  /* norm2 block reduction */
  dimGrid.x  = 1;
  dimBlock.x = chunks;
  cudaDeviceSynchronize();/* make sure everyone is caught up */

  kernRedNorm_double<<<dimGrid,dimBlock,chunks*sizeof(double),xd->stream>>>(devPartial,xd->zval);
  ierr = VecCheckCUDAError("kernRedNorm_double launch in VecNorm_SeqGPU");CHKERRQ(ierr);
  cudaDeviceSynchronize();/* make sure everyone is caught up */

  /* Copy back norm z */
  cudastatus=cudaMemcpy(&zhost,xd->zval,sizeof(double),cudaMemcpyDeviceToHost);/* copy back z */
  ierr = VecCheckCUDAStatus(cudastatus,"on copy znorm D2H in VecNorm_SeqGPU");CHKERRQ(ierr);
  *z = PetscSqrtScalar(zhost);
  xd->vstat.h2d_count++;
  xd->vstat.h2d_bytes+=sizeof(double);
  //printf("Znorm: %e, zhost: %e\n",*z,zhost);

  /* clean up resources */
  cudastatus = cudaFree(devPartial);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  for(i=0;i<chunks;i++){
    cudastatus = cudaStreamDestroy(nrmstream[i]);
    ierr = VecCheckCUDAStatus(cudastatus,"on cudaStreamCreate");CHKERRQ(ierr);
  }
  free(nrmstream);
  cudastatus = cudaFree(devScratch);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


extern __shared__ double arrayNorm_double[];
#undef __FUNCT__
#define __FUNCT__ "kernRedNorm_double"
__global__ void kernRedNorm_double(double* arr,double* z){/* reduction kernel */

  int i = (blockDim.x+1)/2;
  double* zptr=(double*)arrayNorm_double;

  zptr[threadIdx.x]=arr[threadIdx.x];
  __syncthreads();
  while(i>0){
    if(threadIdx.x<i){
      zptr[threadIdx.x]+=zptr[threadIdx.x+i];
    }
    __syncthreads();
    i/=2;
  }/* end while */
  if(threadIdx.x==0){
    *z=zptr[0];
  }
}

#undef __FUNCT__
#define __FUNCT__ "kernNorm2_double"
__global__ void kernNorm2_double(double* devX,
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

  /* block level reduction */
  if(item<n){/* read in values to shared */
    chunkX[threadIdx.x]=devX[item]; /* offset values */
  }else{
    chunkX[threadIdx.x]=0.;
  }

  chunkX[threadIdx.x]*=chunkX[threadIdx.x];
  __syncthreads();

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
  }/* end while */
  if(tid==0){
     *z=scratch[blockIdx.x];
  }
  return;
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
  //printf("Call to VecGetArray_SeqGPU\n");
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
  /* printf("Call to VecRestoreArray_SeqGPU\n"); */
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
  //printf("VecCreateSeqGPU\n");
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
  cudaError_t cudastatus;

  if(d->map->n!=s->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
   }

  //printf("Call to VecCopy_SeqGPU\n");
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

  //cudastatus = cudaMemcpy(dd->length,sd->length,sizeof(int),cudaMemcpyDeviceToDevice);
  cudastatus = cudaMemcpyAsync(dd->zval,sd->zval,sizeof(double),cudaMemcpyDeviceToDevice,dd->stream);
  ierr = VecCopyOverDevice(d,s); CHKERRQ(ierr);
  dd->syncState=sd->syncState;/* synch signal copy */
  cudaDeviceSynchronize();
  //PetscBool same=PETSC_FALSE;
  //ierr = VecCompare_SeqGPU(s,d,&same,0,0);CHKERRQ(ierr);
  //if(!same)SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Vector duplication failed.");
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
  //printf("Call to VecDuplicate_SeqGPU\n");
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
    cudaDeviceSynchronize();
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
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PinnedMalloc"
static PetscErrorCode  PinnedMalloc(PetscScalar** x,PetscInt n){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  cudastatus=cudaHostAlloc((void**)x,n,0);
  ierr=VecCheckCUDAStatus(cudastatus,"in PinnedMalloc");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PinnedFree"
static PetscErrorCode  PinnedFree(PetscScalar* x){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  cudastatus=cudaFreeHost(x);
  ierr=VecCheckCUDAStatus(cudastatus,"in PinnedFree");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecCreate_SeqGPU"
PetscErrorCode  VecCreate_SeqGPU(Vec V){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  PetscMPIInt    size;
  Vec_SeqGPU* seqgpu=PETSC_NULL;
  //printf("Call to VecCreate_SeqGPU\n");
  /*  ierr = PetscNewLog(V,Vec_SeqGPU,&(V->data));CHKERRQ(ierr); */
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
  seqgpu->lifetime       = VEC_PERSIST;

  seqgpu->vstat.h2d_count=0;
  seqgpu->vstat.d2h_count=0;
  seqgpu->vstat.h2d_bytes=0;
  seqgpu->vstat.d2h_bytes=0;

  seqgpu->unplacedarray=PETSC_NULL;
  seqgpu->array_allocated=PETSC_NULL;
  seqgpu->array=PETSC_NULL;
  cudaDeviceSynchronize();

  /* create an associated stream */
  cudastatus = cudaStreamCreate(&(seqgpu->stream));
  ierr = VecCheckCUDAStatus(cudastatus,"on device cudaStreamCreate");CHKERRQ(ierr);

  /* allocate the variable for vector size */
  cudastatus=cudaMalloc((void**)&(seqgpu->length),sizeof(int));
  ierr = VecCheckCUDAStatus(cudastatus,"**** Alloc devlength in VecCreate_SeqGPU");CHKERRQ(ierr);
  /* send vec length size to device */
  cudastatus=cudaMemcpyAsync((void*)seqgpu->length,
               (void*)&(V->map->n),sizeof(int),cudaMemcpyHostToDevice,seqgpu->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"**** Copy H2D devlength in VecCreate_SeqGPU");CHKERRQ(ierr);
  seqgpu->vstat.h2d_count++;
  seqgpu->vstat.h2d_bytes+=sizeof(int);

  /* allocate the vector on device */
  cudastatus=cudaMalloc((void**)&(seqgpu->devptr),V->map->n*sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"**** Alloc of devptr in VecCreate_SeqGPU");CHKERRQ(ierr);

  cudastatus=cudaMemsetAsync((void*)seqgpu->devptr,0,V->map->n*sizeof(double),seqgpu->stream);
  ierr = VecCheckCUDAStatus(cudastatus,"on device cudaMemSet");CHKERRQ(ierr);

  /* allocate the variable for vector offsets */
  cudastatus=cudaMalloc((void**)&(seqgpu->offset),sizeof(int));
  ierr = VecCheckCUDAStatus(cudastatus,"**** Alloc devoffset in VecCreate_SeqGPU");CHKERRQ(ierr);

  /* allocate the variable for vector segment length */
  cudastatus=cudaMalloc((void**)&(seqgpu->segment),sizeof(int));
  ierr = VecCheckCUDAStatus(cudastatus,"**** Alloc dev segment in VecCreate_SeqGPU");CHKERRQ(ierr);

  /* allocate the variable for vector single value result */
  cudastatus=cudaMalloc((void**)&(seqgpu->zval),sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"**** Alloc dev zval in VecCreate_SeqGPU");CHKERRQ(ierr);

  /* using pinned memory */
  ierr = PinnedMalloc(&(seqgpu->cpuptr),V->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
  //ierr = PetscMalloc(V->map->n*sizeof(PetscScalar),&(seqgpu->cpuptr));CHKERRQ(ierr);
  ierr = PetscMemzero(seqgpu->cpuptr,V->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
  seqgpu->syncState=VEC_ALLOC;
  
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "VecDestroy_SeqGPU"
PetscErrorCode VecDestroy_SeqGPU(Vec v){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus;
  Vec_SeqGPU* vd = (Vec_SeqGPU*)v->data;
  PetscValidHeaderSpecific(v,VEC_CLASSID,1);
  //printf("VecDestroy_SeqGPU vstats: \n");
  //printf("...................................\n");
  //printf("H2D transfers: %d, byte count: %d\n",vd->vstat.h2d_count,vd->vstat.h2d_bytes);
  //printf("D2H transfers: %d, byte count: %d\n",vd->vstat.d2h_count,vd->vstat.d2h_bytes);
  //printf("...................................\n");
  /* static int counter = 1; */
  if(vd && vd->syncState != VEC_UNALLOC){
    if(vd->devptr){
      cudastatus=cudaFree(vd->devptr);
      ierr=VecCheckCUDAStatus(cudastatus,"destroying vd->devptr in VecDestroy_SeqGPU");CHKERRQ(ierr);
      vd->devptr=PETSC_NULL;
    }
    if(vd->length){
      cudastatus=cudaFree(vd->length);
      ierr=VecCheckCUDAStatus(cudastatus,"destroying vd->length in VecDestroy_SeqGPU");CHKERRQ(ierr);
      vd->length=PETSC_NULL;
    }
    if(vd->segment){
      cudastatus=cudaFree(vd->segment);
      ierr=VecCheckCUDAStatus(cudastatus,"destroying vd->segment in VecDestroy_SeqGPU");CHKERRQ(ierr);
      vd->segment=PETSC_NULL;
    }
    if(vd->zval){
      cudastatus=cudaFree(vd->zval);
      ierr=VecCheckCUDAStatus(cudastatus,"destroying vd->zval in VecDestroy_SeqGPU");CHKERRQ(ierr);
      vd->zval=PETSC_NULL;
    }
    if(vd->cpuptr){
      ierr = PinnedFree(vd->cpuptr); CHKERRQ(ierr);
    }
    vd->syncState = VEC_UNALLOC;
  }

  cudastatus = cudaStreamDestroy(vd->stream);
  ierr = VecCheckCUDAError("call to cudaStreamDestroy");CHKERRQ(ierr);

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
