#include <petscconf.h>
#include <petscsys.h>
PETSC_CUDA_EXTERN_C_BEGIN
#include <string.h>
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
    fprintf(stderr, "Cuda error: %s: %s.\n", msg,cudaGetErrorString(err));
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
      fprintf(stderr, "Cuda error!: %s: %s.\n",msg,cudaGetErrorString(cs));
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

  printf("Call to VecSetRandom_SeqGPU\n");
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


/*

#undef __FUNCT__
#define __FUNCT__ "VecSetRandom_SeqGPU"
PetscErrorCode VecSetRandom_SeqGPU(Vec x,PetscRandom r){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSetRandom_Seq(x,r);CHKERRQ(ierr);
  printf("Call to VecSetRandom_SeqGPU (***EMPTY***)\n");
  PetscFunctionReturn(0);
}
*/




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





/*

#undef __FUNCT__
#define __FUNCT__ "VecResetArray_Seq"
static PetscErrorCode VecResetArray_Seq(Vec vin){
  //Vec_Seq *v = (Vec_Seq *)vin->data;

  PetscFunctionBegin;
  printf("Call to VecResetArray_Seq\n");
  PetscFunctionReturn(0);
}

*/





#undef __FUNCT__
#define __FUNCT__ "VecSetValues_SeqGPU"
/*@
   VecSetValues - Inserts or adds values into certain locations of a vector.
@*/
PetscErrorCode VecSetValues_SeqGPU(Vec x,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode iora){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt i;
  PetscScalar yval=0;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  //printf("Call to VecSetValues_SeqGPU\n");
  if(xd->syncState==VEC_CPU || xd->syncState==VEC_SYNCHED){
    if(iora==INSERT_VALUES){
      for(i=0;i<ni;i++){
         yval = y[i];
         xd->cpuptr[i]=yval;
      }
      ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
      xd->syncState=VEC_SYNCHED;
    }else{
      /* ADD_VALUES not supported now */
      printf("Call to VecSetValues_SeqGPU: ADD_VALUES (*** EMPTY ***)\n");
    }
  }else{
      if(iora==INSERT_VALUES){/* not efficient */
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
  PetscInt bx, tx;
  Vec_SeqGPU* dd = (Vec_SeqGPU*)d->data;
  Vec_SeqGPU* sd = (Vec_SeqGPU*)s->data;
  //printf("Call to VecCopyOverDevice\n");
  dim3 dimGrid;
  dim3 dimBlock;

  if(s->map->n!=d->map->n){
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
  }

  if(sd->syncState==VEC_CPU){/* synch y */
    ierr = VecCopyOverH2D(s,sd->cpuptr);CHKERRQ(ierr);
    sd->syncState=VEC_SYNCHED;
  }

  /* assuming width mem load isn't going to be an issue */
  bx=ceil((float)d->map->n/(float)CPYTCOUNT);
  tx=CPYTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;
  kernCODevice<<<dimGrid,dimBlock>>>(dd->devptr,sd->devptr,sd->length);
  ierr = VecCheckCUDAError("kernel call to kernCODevice");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "kernCODevice"
__global__ void kernCODevice(double* devY,double* devX, int *n){

  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  __shared__ double chunkX[CPYTCOUNT];
  if(tid<*n){
    //if(devX[tid]!=0)printf("devX[%d]: %e, len: %d\n",tid,devX[tid],*n);
    chunkX[threadIdx.x]=devX[tid];
    devY[tid]=chunkX[threadIdx.x];
    //if(devY[tid]!=0)printf("devY[%d]: %e, len: %d\n",tid,devY[tid],*n);
  }
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
  cudastatus=cudaMemcpy(&(vd->devptr[offset]),y,blocksize*sizeof(PetscScalar),cudaMemcpyHostToDevice);
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
  cudastatus=cudaMemcpy(vd->devptr,y,v->map->n*sizeof(PetscScalar),cudaMemcpyHostToDevice);
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
  cudastatus=cudaMemcpy(y,&(vd->devptr[offset]),blocksize*sizeof(PetscScalar),cudaMemcpyDeviceToHost);
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
  cudastatus=cudaMemcpy(y,vd->devptr,v->map->n*sizeof(PetscScalar),cudaMemcpyDeviceToHost);
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H in VecCopyOverD2H");CHKERRQ(ierr);
  vd->vstat.d2h_count++;
  vd->vstat.d2h_bytes+=v->map->n*sizeof(PetscScalar);
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VecSet_SeqGPU"
PetscErrorCode VecSet_SeqGPU(Vec xin,PetscScalar alpha){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //PetscInt i=0;
  cudaError_t cudastatus;
  dim3 dimgrid(ceil((float)xin->map->n/((float)TCOUNT)),1,1);
  dim3 dimblocks(TCOUNT,1,1);
  Vec_SeqGPU* xd = (Vec_SeqGPU*)xin->data;
  //printf("Call to VecSet_SeqGPU\n");
  if(xd->syncState==VEC_UNALLOC){
    SETERRQ(PETSC_COMM_SELF,
            PETSC_ERR_MEM,"*** In VecSet_SeqGPU, Vec not allocated. ***\n");
  }else{
    if(alpha==0){
      cudastatus=cudaMemset((void*)xd->devptr,0,xin->map->n*sizeof(double));
      ierr = VecCheckCUDAStatus(cudastatus,"on device cudaMemSet VecSet_SeqGPU");CHKERRQ(ierr);
    }else{
      cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
      ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
      kernSet<<<dimgrid,dimblocks>>>(xd->devptr,xd->length);
      ierr = VecCheckCUDAError("Call to kernSet. "); CHKERRQ(ierr);
    }
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

  if(alpha==0){
    cudastatus=cudaMemset((void*)xd->devptr,0,x->map->n*sizeof(double));
    ierr = VecCheckCUDAStatus(cudastatus,"on device cudaMemSet VecSet_SeqGPU");CHKERRQ(ierr);
  }else if (alpha != 1.0){
    cudastatus=cudaMemcpyToSymbol("dblScalarValue",(void*)&alpha,sizeof(double),0,cudaMemcpyHostToDevice);
    ierr = VecCheckCUDAStatus(cudastatus,"error in symbol copy to device");CHKERRQ(ierr);
    kernScale<<<dimgrid,dimblocks>>>(xd->devptr,xd->length);
    ierr = VecCheckCUDAError("Call to kernScale. "); CHKERRQ(ierr);
  }
  cudaDeviceSynchronize();
  fflush(NULL);
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






#undef __FUNCT__
#define __FUNCT__ "VecCheck_SeqGPU"
PetscErrorCode VecCheck_SeqGPU(Vec x){
  PetscFunctionBegin;/*
  PetscErrorCode ierr;
  dim3 dimgrid(ceil((float)x->map->n/((float)TCOUNT)),1,1);
  dim3 dimblocks(TCOUNT,1,1);
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  printf("******************************************\n");
  kernCheck<<<dimgrid,dimblocks>>>(xd->devptr,xd->length);
  ierr = VecCheckCUDAError("Call to kernCheck. "); CHKERRQ(ierr);
  cudaDeviceSynchronize();
  printf("******************************************\n");
  fflush(NULL);*/
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "kernCheck"
__global__ void kernCheck(double* x, int* n){
  int tid = threadIdx.x + blockDim.x*blockIdx.x;
  if(tid<*n){
    if(x[tid]!=0)printf("kernCheck: x[%d]: %e, length: %d\n",tid,x[tid],*n);
  }
}








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
  double *devS;
  PetscScalar *s=PETSC_NULL;
  PetscInt i,blks,thds;
  //printf("Call to VecDot_SeqGPU\n");
  blks = ceil((float)x->map->n/(float)DOTTCOUNT);
  thds = DOTTCOUNT;
  //  printf("Blocks: %d, Threads: %d\n", blks,thds);
  ierr = PetscMalloc(blks*sizeof(PetscScalar),(void**)&s);CHKERRQ(ierr);
  cudastatus=cudaMalloc((void**)&devS,blks*sizeof(double));/* could probably make this static */
  ierr = VecCheckCUDAStatus(cudastatus,"s alloc in VecDot_SeqGPU");CHKERRQ(ierr);
  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)y->data;
  dim3 dimGrid(blks,1);
  dim3 dimBlock(thds,1);

  if(xd->syncState==VEC_CPU){
    printf("xd state VEC_CPU: copying to device.\n");
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }

  if(yd->syncState==VEC_CPU){
    printf("yd state VEC_CPU: copying to device.\n");
    ierr = VecCopyOverH2D(y,yd->cpuptr);CHKERRQ(ierr);
    yd->syncState=VEC_SYNCHED;
  }

  kernDot<<<dimGrid,dimBlock>>>((double*)xd->devptr,(double*)yd->devptr,(int*)xd->length,(double*)devS);
  ierr = VecCheckCUDAError("kern launch in VecDot_SeqGPU");CHKERRQ(ierr);


  /* implicit barrier here */
  cudastatus=cudaMemcpy(s,devS,blks*sizeof(PetscScalar),cudaMemcpyDeviceToHost);/* copy back s */
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H in VecDot_SeqGPU");CHKERRQ(ierr);

  *z=0;
  for(i=0;i<blks;i++){/* last reduction done on CPU */
    *z+=s[i];
  }
  //printf("dot product Z: %e\n",*z,s[0],x->map->n);
  ierr = PetscFree(s);CHKERRQ(ierr);
  cudastatus = cudaFree(devS);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}






__global__ void kernDot(double* devX, double* devY, int* n, double* s){
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  //int localn=*n;
  int i = (DOTTCOUNT+1)/2;

  __shared__ double chunkX[DOTTCOUNT];
  __shared__ double chunkY[DOTTCOUNT];

  if(tid<*n){
    /* read in values to shared */
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkX[threadIdx.x]*=chunkY[threadIdx.x];
  }else{
    chunkX[threadIdx.x]=0;
  }
  __syncthreads();
  while(i>0){
     if(threadIdx.x<i){
       chunkX[threadIdx.x]+=chunkX[threadIdx.x+i];
     }
     __syncthreads();
     i/=2;
  }/* end while */
  if(threadIdx.x==0){
    s[blockIdx.x]=chunkX[0];
  }
}





#undef __FUNCT__
#define __FUNCT__ "VecMDot_SeqGPU"
/*@
   VecMDot - Computes vector multiple dot products. 

   Collective on Vec

   Input Parameters:
+  x - one vector
.  nv - number of vectors
-  y - array of vectors. 

   Output Parameter:
.  val - array of the dot products (does not allocate the array)

@*/
PetscErrorCode  VecMDot_SeqGPU(Vec x,PetscInt nv,const Vec y[],PetscScalar val[]){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt       i;
  // printf("VecMDot_SeqGPU\n");
  for (i=0; i<nv; i++) {
    ierr = VecDot_SeqGPU(x,y[i],&val[i]);CHKERRQ(ierr);
    if(PetscIsInfOrNanScalar(val[i])){
      SETERRQ1(((PetscObject)x)->comm,PETSC_ERR_FP,"Infinite or not-a-number generated in mdot, entry %D",i);
    }
  }
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

  kernCopyLen<<<1,1>>>(dd->length,sd->length);
  ierr = VecCheckCUDAError("call to kernCopyLen");CHKERRQ(ierr);

  ierr = VecCopyOverDevice(d,s); CHKERRQ(ierr);
  dd->syncState=sd->syncState;/* synch signal copy */

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
#define __FUNCT__ "VecAXPBY_SeqGPU"
PetscErrorCode VecAXPBY_SeqGPU(Vec yin,PetscScalar beta,PetscScalar alpha,Vec xin){
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
  cudaDeviceSynchronize();
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
    __syncthreads();
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
    __syncthreads();
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
;
  /* printf("in kernWAXPY:alphaShared: %f, tid: %d, vlen: %d\n",alphaShared,tid,*vlen); */
  if(tid<*vlen){
    //if(devX[tid]!=0)printf("kernWXMY: devX[%d]: %e\n",tid,devX[tid]);
    //if(devY[tid]!=0)printf("kernWXMY: devY[%d]: %e\n",tid,devY[tid]);
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkW[threadIdx.x]=chunkY[threadIdx.x]-chunkX[threadIdx.x];
    __syncthreads();
    devW[tid]=chunkW[threadIdx.x];
    //if(devW[tid]!=0)printf("kernWXMY: devW[%d]: %e\n",tid,devW[tid]);
  }
}

















#undef __FUNCT__
#define __FUNCT__ "VecMAXPY_SeqGPU"
PetscErrorCode VecMAXPY_SeqGPU(Vec x,PetscInt nv,const PetscScalar* alpha,Vec *y){
  /* y = y + sum(a[i]*x[i]) */
  PetscFunctionBegin;
  //printf("VecMAXPY_SeqGPU\n");
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
  cudastatus=cudaMemset((void*)devW,0,x->map->n*sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaMemset");CHKERRQ(ierr);

  /* assuming xwidth mem load isn't going to be an issue */

  bx=ceil((float)x->map->n/(float)AXPYTCOUNT);
  tx=AXPYTCOUNT;
  dimGrid.x=bx; dimGrid.y=1;
  dimBlock.x=tx; dimBlock.y=1;

  for(i=0;i<nv;i++){
    if(y[i]->map->n!=x->map->n){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"Vector size mismatch.");
    }
    yd=(Vec_SeqGPU*)y[i]->data;
    if(yd->syncState==VEC_CPU){/* synch x */
      ierr = VecCopyOverH2D(y[i],yd->cpuptr);CHKERRQ(ierr);
      yd->syncState=VEC_SYNCHED;
    }
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
    cudaDeviceSynchronize();
  }

  if(xd->syncState==VEC_CPU){/* synch x */
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }
  kernXPY<<<dimGrid,dimBlock>>>(xd->devptr,devW,xd->length);
  ierr = VecCheckCUDAError("kernel call to kernXPY");CHKERRQ(ierr);
  cudaDeviceSynchronize();

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
    //if(devX[tid]!=0)printf("kernXPY: devX[%d]: %e\n",tid,devX[tid]);
    //if(devX[tid]!=0)printf("kernXPY: PRE: devY[%d]: %e\n",tid,devY[tid]);
    chunkY[threadIdx.x]=devY[tid];
    chunkY[threadIdx.x]+=chunkX[threadIdx.x];
    devY[tid]=chunkY[threadIdx.x];
    //if(devY[tid]!=0)printf("kernXPY: POST: devY[%d]: %e\n",tid,devY[tid]);
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
  //if(tid==0)printf("in kernAXPY:alphaShared: %f\n",alphaShared);
  if(tid<*vlen){
    chunkX[threadIdx.x]=devX[tid];
    //if(devX[tid]!=0)printf("kernAXPY: devIN[%d]: %e, len: %d\n",tid,devX[tid],*vlen);
    //if(devY[tid]!=0)printf("kernAXPY: PRE: devOUT[%d]: %e\n",tid,devY[tid]);
    chunkY[threadIdx.x]=devY[tid];
    chunkY[threadIdx.x]+=chunkX[threadIdx.x]*alphaShared;
    devY[tid]=chunkY[threadIdx.x];
    //if(devY[tid]!=0)printf("kernAXPY: POST: devOUT[%d]: %e\n",tid,devY[tid]);
  }
}



#undef __FUNCT__
#define __FUNCT__ "VecPointwiseMult_SeqGPU"
PetscErrorCode VecPointwiseMult_SeqGPU(Vec w,Vec x,Vec y){
  PetscFunctionBegin;
  //printf("VecPointwiseMult_SeqGPU\n");
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
  //printf("VecMaxPointwiseDivide_SeqGPU...");
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
  //printf("max: %f\n",*max);
  ierr = PetscFree(maxlist);CHKERRQ(ierr);
  cudastatus = cudaFree(devmaxlist);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
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
    //printf("kernMAXDIV: devX[%d]: %e\n",tid,devX[tid]);
    //printf("kernMAXDIV: devY[%d]: %e\n",tid,devY[tid]);
    if(chunkY[threadIdx.x]!=0){
      chunkW[threadIdx.x]=fabs(__ddiv_rn(chunkX[threadIdx.x],chunkY[threadIdx.x]));
    }else{
      chunkW[threadIdx.x]=fabs(chunkX[threadIdx.x]);
    }
    //printf("kernMAXDIV: d[%d]: %e\n",threadIdx.x,chunkW[threadIdx.x]);

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
  //printf("Call to VecPointwiseDivide_SeqGPU\n");
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




#undef __FUNCT__
#define __FUNCT__ "VecNorm_SeqGPU"
PetscErrorCode VecNorm_SeqGPU(Vec x,NormType type,PetscReal* z){
  /* NormType: NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4 */
  /* dealing with NORM_2 for now... */
  PetscFunctionBegin;
  PetscErrorCode ierr;
  cudaError_t cudastatus=cudaSuccess;
  double *deviceS=PETSC_NULL;
  PetscScalar *s=PETSC_NULL;
  PetscInt i,blks,thds;
  //printf("VecNorm_SeqGPU\n");
  blks = ceil((float)x->map->n/(float)NORMTCOUNT);
  thds = NORMTCOUNT;
  //printf("Blocks: %d, Threads: %d\n", blks,thds);
  ierr = PetscMalloc(blks*sizeof(PetscScalar),(void**)&s);CHKERRQ(ierr);
  cudaDeviceSynchronize();

  cudastatus=cudaMalloc((void**)&deviceS,blks*sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"deviceS alloc in VecNorm_SeqGPU");CHKERRQ(ierr);

  Vec_SeqGPU *xd=(Vec_SeqGPU*)x->data;
  dim3 dimGrid(blks,1);
  dim3 dimBlock(thds,1);

  if(xd->syncState==VEC_CPU){
    ierr = VecCopyOverH2D(x,xd->cpuptr);CHKERRQ(ierr);
    xd->syncState=VEC_SYNCHED;
  }

  kernNorm2<<<dimGrid,dimBlock>>>((double*)xd->devptr,(int*)xd->length,(double*)deviceS);
  ierr = VecCheckCUDAError("kern launch in VecDot_SeqGPU");CHKERRQ(ierr);

  /* implicit device barrier below */
  cudastatus=cudaMemcpy(s,deviceS,blks*sizeof(PetscScalar),cudaMemcpyDeviceToHost);/* copy back s */
  ierr = VecCheckCUDAStatus(cudastatus,"on copy D2H in VecDot_SeqGPU");CHKERRQ(ierr);

  cudaDeviceSynchronize();
  cudastatus = cudaFree(deviceS);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);

  *z=0;
  for(i=0;i<blks;i++){/* last reduction done on CPU */
    *z+=s[i];
  }
  ierr = PetscFree(s);CHKERRQ(ierr);
  z[0]=sqrt(z[0]);/* norm2 in 0 norm1_2 in 0 and 1 respectively */
  //printf("NORM2: %e\n",z[0]);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "kernNorm2"
__global__ void kernNorm2(double* devX,int* n, double* s){
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  //int localn=*n;
  int i = (NORMTCOUNT+1)/2;

  __shared__ double chunkX[NORMTCOUNT];
  if(tid<*n){
    /* read in values to shared */
    //if(devX[tid]!=0)printf("kernNorm2: devX[%d]: %e\n",tid,devX[tid]);
    chunkX[threadIdx.x]=devX[tid];
    chunkX[threadIdx.x]*=chunkX[threadIdx.x];
  }else{
    chunkX[threadIdx.x]=0;
  }

  __syncthreads();

  //if(chunkX[threadIdx.x]>0)printf("chunkX[%d]: %e\n",threadIdx.x,chunkX[threadIdx.x]);
  while(i>0){
    if(threadIdx.x<i){
       chunkX[threadIdx.x]+=chunkX[threadIdx.x+i];
    }
    __syncthreads();
    i/=2;
  }/* end while */
  if(threadIdx.x==0){
    s[blockIdx.x]=chunkX[0];
    //printf("devS[%d]: %e\n",blockIdx.x,s[blockIdx.x]);
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
  z[0]=sqrt(z[0]);
  printf("ZNORM: %f\n\n",*z);
  PetscFunctionReturn(0);
}*/


/*
#undef __FUNCT__
#define __FUNCT__ "kernReduceAbsSum"
PetscErrorCode kernReduceAbsSum(double * x, PetscReal* z){

}
*/



#undef __FUNCT__
#define __FUNCT__ "VecResetArray_SeqGPU"
PetscErrorCode VecResetArray_SeqGPU(Vec vin){
  PetscFunctionBegin;
  printf("VecResetArray_SeqGPU (***EMPTY***)\n");
  PetscFunctionReturn(0);
}




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

  int i,flg1=0,flg2=0;
  PetscStackCheckByName(4,"DMDAVecGetArray",flg1);
  PetscStackCheckByName(6,"DMGlobalToLocalBegin",flg2);
  if(vd->syncState==VEC_GPU || flg1 || flg2){/* may not need VEC_GPU element */
    ierr = VecCopyOverD2H(v,vd->cpuptr); CHKERRQ(ierr);
    vd->syncState = VEC_CPU;
    //for(i=0;i<v->map->n;i++){
    //   if(vd->cpuptr[i]!=0)printf("get[%d]: %e\n",i,vd->cpuptr[i]);
    //}
  }
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
    int i=0;
    if(a){
      ierr = VecCopyOverH2D(v,*a);CHKERRQ(ierr);
      //for(i=0;i<v->map->n;i++){
      //  if(vd->cpuptr[i]!=0)printf("put *a[%d]: %e\n",i,(*a)[i]);
      //}
      vd->syncState=VEC_GPU;
    }else{
      ierr = VecCopyOverH2D(v,vd->cpuptr);CHKERRQ(ierr);
      // for(i=0;i<v->map->n;i++){
      //   if(vd->cpuptr[i]!=0)printf("put cpuptr[%d]: %e\n",i,vd->cpuptr[i]);
      //}
      vd->syncState=VEC_SYNCHED;
    }
  }
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "VecPlaceArray_SeqGPU"
PetscErrorCode VecPlaceArray_SeqGPU(Vec vin,const PetscScalar *a){
  //PetscErrorCode ierr;
  PetscFunctionBegin;
  printf("VecPlaceArray_SeqGPU (***EMPTY***)\n");
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "VecReplaceArray_SeqGPU"
PetscErrorCode VecReplaceArray_SeqGPU(Vec vin,const PetscScalar *a){
  //PetscErrorCode ierr;
  PetscFunctionBegin;
  printf("VecReplaceArray_SeqGPU (***EMPTY***)\n");
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
#define __FUNCT__ "VecDotNorm2_SeqGPU"
PetscErrorCode VecDotNorm2_SeqGPU(Vec s, Vec t, PetscScalar *dp, PetscScalar *nm){
  /* PetscErrorCode ierr; */
  /* PetscScalar zero = 0.0,n=s->map->n; */
  PetscFunctionBegin;
  printf("VecDotNorm2_SeqGPU (***EMPTY***)\n");
  /* ierr = PetscLogFlops(4.0*n);CHKERRQ(ierr); */
  PetscFunctionReturn(0);
}






/*

#undef __FUNCT__
#define __FUNCT__ "VecAXPBYPCZ_SeqGPU"
PetscErrorCode VecAXPBYPCZ_SeqGPU(Vec zin,PetscScalar alpha,PetscScalar beta,PetscScalar gamma,Vec xin,Vec yin)
{
  PetscErrorCode     ierr;
  PetscInt           n = zin->map->n;

  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
*/

#undef __FUNCT__
#define __FUNCT__ "VecAXPBYPCZ_SeqGPU"
PetscErrorCode VecAXPBYPCZ_SeqGPU(Vec x, PetscScalar alpha, PetscScalar beta,\
                           PetscScalar gamma, Vec y, Vec z){

  PetscFunctionBegin;
  PetscErrorCode ierr;
  int blocks=ceil(x->map->n/32);/* assuming shared memory size is not an issue */
  int threads=32;
  cudaError_t cudastatus;
  Vec_SeqGPU* devX = (Vec_SeqGPU*)x->data;
  Vec_SeqGPU* devY = (Vec_SeqGPU*)y->data;
  Vec_SeqGPU* devZ = (Vec_SeqGPU*)z->data;

  double2 params[2];
  double2 *devparams=PETSC_NULL;

  cudastatus = cudaMalloc((void**)&devparams,2*sizeof(double2));
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaMalloc()");CHKERRQ(ierr);
  params[0].x=alpha;
  params[0].y=beta;
  params[1].x=gamma;
  cudastatus=cudaMemcpy(devparams,params,2,cudaMemcpyHostToDevice);
  dim3 dimGrid; dimGrid.x=blocks; dimGrid.y=1;
  dim3 dimBlock; dimBlock.x=threads; dimBlock.y=1;
  kernAXPBYPCZ<<<dimGrid,dimBlock>>>(devparams,devX->devptr,devY->devptr,devZ->devptr);
  cudastatus = cudaFree(devparams);
  ierr = VecCheckCUDAStatus(cudastatus,"on cudaFree()");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}






__global__ void kernAXPBYPCZ(double2* devparams, double* devX, double* devY, double* devZ){
  /* x <- alpha*x + beta*y + gamma*z */
  int tid = blockIdx.x*blockDim.x+threadIdx.x;
  double2 alphabeta = devparams[0];
  double2 gamma = devparams[1];
  double work;
  int localn = devN;

  __shared__ double chunkX[32];
  __shared__ double chunkY[32];
  __shared__ double chunkZ[32];

  if(tid<localn){
    /* read in values to shared */
    chunkX[threadIdx.x]=devX[tid];
    chunkY[threadIdx.x]=devY[tid];
    chunkZ[threadIdx.x]=devZ[tid];

    /* do flops */
    if(gamma.x){
      work=gamma.x*chunkZ[threadIdx.x];
    }
    if(alphabeta.y){
      work+=alphabeta.y*chunkY[threadIdx.x];
    }
    if(alphabeta.x){
      work+=alphabeta.x*chunkX[threadIdx.x];
    }
    chunkX[threadIdx.x]+=work;

    /* write back */
    devX[tid]=chunkX[threadIdx.x];
  }
  return;
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
  }
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
  printf("VecDestroy_SeqGPU vstats: \n");
  printf("...................................\n");
  printf("H2D transfers: %d, byte count: %d\n",vd->vstat.h2d_count,vd->vstat.h2d_bytes);
  printf("D2H transfers: %d, byte count: %d\n",vd->vstat.d2h_count,vd->vstat.d2h_bytes);
  printf("...................................\n");
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
    if(vd->cpuptr){
      ierr = PetscFree(vd->cpuptr);CHKERRQ(ierr);
    }
    vd->syncState = VEC_UNALLOC;
  }

  ierr = PetscObjectDepublish(v);CHKERRQ(ierr);
#if defined(PETSC_USE_LOG)
  PetscLogObjectState((PetscObject)v,"Length=%D",v->map->n);
#endif
  ierr = PetscFree(v->data);CHKERRQ(ierr);
  //  ierr = VecDestroy_Seq(v);CHKERRQ(ierr);
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
#define __FUNCT__ "VecSetDimensions_SeqGPU"
PetscErrorCode  VecSetDimensions_SeqGPU(Vec x,PetscInt ndims,dim3 dimsize){
  PetscFunctionBegin;
  Vec_SeqGPU* xd = (Vec_SeqGPU*)x->data;
  if(ndims==3){
    xd->ndims=ndims;
    xd->dimsize.x=dimsize.x;
    xd->dimsize.y=dimsize.y;
    xd->dimsize.z=dimsize.z;
  }else if(ndims==2){
    xd->ndims=ndims;
    xd->dimsize.x=dimsize.x;
    xd->dimsize.y=dimsize.y;
    xd->dimsize.z=0;
  }else if(ndims==1){
    xd->ndims=ndims;
    if(x->map->n!=dimsize.x){
      SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"VecSetDimensions_SeqGPU does not support memory resizing.");
    }else{
      xd->dimsize.x=dimsize.x;
      xd->dimsize.y=0;
      xd->dimsize.z=0;
    }
  }else{
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_MEM,"VECSEQGPU doesn't suport given number of dimensions.");
  }
  xd->dimsetflag=PETSC_TRUE;
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




  /* allocate the variable for vector size */
  cudastatus=cudaMalloc((void**)&(seqgpu->length),sizeof(int));
  ierr = VecCheckCUDAStatus(cudastatus,"**** Alloc devlength in VecCreate_SeqGPU");CHKERRQ(ierr);

  /* send vec length size to device */
  cudastatus=cudaMemcpy((void*)seqgpu->length,(void*)&(V->map->n),sizeof(int),cudaMemcpyHostToDevice);
  ierr = VecCheckCUDAStatus(cudastatus,"**** Copy H2D devlength in VecCreate_SeqGPU");CHKERRQ(ierr);
  seqgpu->vstat.h2d_count++;
  seqgpu->vstat.h2d_bytes+=sizeof(int);

  /* allocate the vector on device */
  cudastatus=cudaMalloc((void**)&(seqgpu->devptr),V->map->n*sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"**** Alloc of devptr in VecCreate_SeqGPU");CHKERRQ(ierr);
  cudastatus=cudaMemset((void*)seqgpu->devptr,0,V->map->n*sizeof(double));
  ierr = VecCheckCUDAStatus(cudastatus,"on device cudaMemSet");CHKERRQ(ierr);

  seqgpu->ndims=1;/* default number of dimensions */
  seqgpu->dimsize.x=V->map->n;
  seqgpu->dimsize.y=0;
  seqgpu->dimsize.z=0;
  seqgpu->dimsetflag=PETSC_FALSE;

  ierr = PetscMalloc(V->map->n*sizeof(PetscScalar),&(seqgpu->cpuptr));CHKERRQ(ierr);
  ierr = PetscMemzero(seqgpu->cpuptr,V->map->n*sizeof(PetscScalar));CHKERRQ(ierr);
  seqgpu->syncState=VEC_ALLOC;

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
