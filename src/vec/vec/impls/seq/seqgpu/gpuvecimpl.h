#ifndef __GPUVECIMPL
#define __GPUVECIMPL

#include <petsc-private/vecimpl.h>
#include <petscconf.h>
#include <petscsys.h>
#include <petscvec.h>

#include <../src/vec/vec/impls/dvecimpl.h>

EXTERN_C_BEGIN
#include <cuda.h>

#undef  DEBUGVEC
#define DEBUGVEC        0
#undef  VERBOSE
#define VERBOSE         0
#define MPLIER          4.0
#define CHUNKWIDTH      MPLIER*65536.0
#define TCOUNT          128
#define MAXTHREADS      1024
#define MAXBLOCKS       65536
#define THRDOTCNT       128
#define THRDOTCNT2      32
#define AXPYTCOUNT      128
#define AXPBYPCZTCOUNT  128
#define THRNRMCNT       128
#define THRNRMCNT2      32
#define PDIVTCOUNT      128
#define PDIVTCOUNT2     32
#define PMULTCOUNT      128
#define CPYTCOUNT       128
#define NNN		624
#define MMM		397
#define INIT_MULT	1812433253	/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
#define	ARRAY_SEED	19650218	/* Seed for initial setup before incorp array seed */
#define MATRIX_A	0x9908b0df	/* Constant vector a */
#define UPPER_MASK	0x80000000	/* Most significant w-r bits */
#define LOWER_MASK	0x7fffffff	/* Least significant r bits */
#define	TEMPER1		0x9d2c5680
#define	TEMPER2		0xefc60000




__shared__ int	mtNexts;			/* Start of next block of seeds */
__shared__ uint	mtNexti;	/* Indirect on above to save one global read time/call */
__device__ uint	s_seeds[NNN];




/*
struct copyCounters{
  PetscInt h2d_count;
  PetscInt d2h_count;
  PetscInt h2d_bytes;
  PetscInt d2h_bytes;
};
 */

static PetscErrorCode VecView_Seq_ASCII(Vec ,PetscViewer);
static PetscErrorCode PinnedMalloc(PetscScalar** x,PetscInt n);
static PetscErrorCode PinnedFree(PetscScalar* x);
extern PetscErrorCode VecDotNorm2_SeqGPU(Vec,Vec,PetscScalar *, PetscScalar *);
extern PetscErrorCode VecPointwiseDivide_SeqGPU(Vec,Vec,Vec);
extern PetscErrorCode VecMaxPointwiseDivide_SeqGPU(Vec,Vec,PetscReal*);
extern PetscErrorCode VecWAXPY_SeqGPU(Vec,PetscScalar,Vec,Vec);
extern PetscErrorCode VecMDot_SeqGPU(Vec,PetscInt,const Vec[],PetscScalar *);
extern PetscErrorCode VecSet_SeqGPU(Vec,PetscScalar);
extern PetscErrorCode VecMAXPY_SeqGPU(Vec,PetscInt,const PetscScalar *,Vec *);
extern PetscErrorCode VecAXPBYPCZ_SeqGPU(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
extern PetscErrorCode VecPointwiseMult_SeqGPU(Vec,Vec,Vec);
extern PetscErrorCode VecPlaceArray_SeqGPU(Vec,const PetscScalar *);
extern PetscErrorCode VecResetArray_SeqGPU(Vec);
extern PetscErrorCode VecReplaceArray_SeqGPU(Vec,const PetscScalar *);
extern PetscErrorCode VecDot_SeqGPU(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecTDot_SeqGPU(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecScale_SeqGPU(Vec,PetscScalar);
extern PetscErrorCode VecCopy_SeqGPU(Vec,Vec);
extern PetscErrorCode VecSwap_SeqGPU(Vec,Vec);
extern PetscErrorCode VecAXPY_SeqGPU(Vec,PetscScalar,Vec);
extern PetscErrorCode VecAXPBY_SeqGPU(Vec,PetscScalar,PetscScalar,Vec);
extern PetscErrorCode VecDuplicate_SeqGPU(Vec,Vec *);
extern PetscErrorCode VecNorm_SeqGPU(Vec,NormType,PetscReal*);
extern PetscErrorCode VecCreate_SeqGPU(Vec);
extern PetscErrorCode VecView_SeqGPU(Vec,PetscViewer);
extern PetscErrorCode VecDestroy_SeqGPU(Vec);
extern PetscErrorCode VecDestroyVecs_SeqGPU(PetscInt,Vec*);
//extern PetscErrorCode VecAYPX_SeqGPU(Vec,PetscScalar,Vec);
extern PetscErrorCode VecSetRandom_SeqGPU(Vec,PetscRandom);
extern PetscErrorCode VecSetValues_SeqGPU(Vec,PetscInt,const PetscInt*,const PetscScalar *y, InsertMode);
extern PetscErrorCode VecCopyOverH2D(Vec,PetscScalar*);
extern PetscErrorCode VecCopyOverD2H(Vec,PetscScalar*);
extern PetscErrorCode VecCopyBlockD2H(Vec,PetscScalar*,PetscInt,PetscInt);
extern PetscErrorCode VecCopyBlockH2D(Vec,PetscScalar*,PetscInt,PetscInt);
extern PetscErrorCode VecCopyOverDevice(Vec,Vec);
extern PetscErrorCode VecCopyBlockDevice(Vec, Vec, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode VecView_SeqGPU(Vec,PetscViewer);
extern PetscErrorCode VecGPUAllocateCheck_Public(Vec);
extern PetscErrorCode VecCopyToGPUSome_Public(Vec,PetscInt);
extern PetscErrorCode VecGetArray_SeqGPU(Vec,PetscScalar**);
extern PetscErrorCode VecRestoreArray_SeqGPU(Vec,PetscScalar**);
extern PetscErrorCode VecCheckCUDAError(const char *);
extern PetscErrorCode VecCheckCUDAStatus(cudaError_t ,const char *);
extern PetscErrorCode VecGetLocalSize_SeqGPU(Vec , PetscInt *);
extern PetscErrorCode VecGetSize_SeqGPU(Vec , PetscInt *);
extern PetscErrorCode VecCompare_SeqGPU(Vec,Vec, PetscBool*, PetscInt, PetscInt);
extern PetscErrorCode VecCheck_SeqGPU(Vec);
extern __global__ void kernAXPBYPCZ(double*, double*, double*,int*);
extern __global__ void kernCheck(double*, int*);
extern __global__ void kernCopyLen(int*, int*);
extern __global__ void kernCODevice(double*,double*, int*);
extern __global__ void kernCompare(double*, double*,int*, int*, int*);
extern __global__ void kernSet(double*, int*);
extern __global__ void kernScale(double*, int*);
extern __global__ void kernCopy(double*, double*, int*, int*);
extern __global__ void kernDot(double*,double*,int*,int*,int*,double*);
extern __global__ void kernRedDot(int*,double*,double*);
extern __global__ void kernAXPY(double*, double*, int*,double*);
extern __global__ void kernWAXPY(double*, double*, int*, double*);
extern __global__ void kernWXPY(double*, double*, int*, double*);
extern __global__ void kernWXMY(double*, double*, int*, double*);
extern __global__ void kernXPY(double*, double*, int*);
extern __global__ void kernNorm2(double*,int*,int*,int*,double*);
extern __global__ void kernRedNorm(int*,double*,double*);
extern __global__ void kernPDIV(double*,double*,int*,double* );
extern __global__ void kernPMULT(double*,double*,int*,double* );
extern __global__ void kernMAXPDIV(double*,double*,int*,int*,int*,double*);
extern __global__ void kernMAX(int*,double*,double*);
extern __global__ void kernRand(double*, int*);
extern __global__ void kernRandS(uint *);

__device__ void mt19937si(uint);
__device__ void mt19937sai(uint*,uint);
__device__ uint mt19937s(void);
__device__ uint mt19937sl(void);


extern PetscBool  synchronizeGPU;



#define CHKERRGPU(err) if (((int)err) != (int)CUBLAS_STATUS_SUCCESS) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"CUDA error %s",cudaGetErrorString(err))

#define WaitForGPU() synchronizeGPU ? cudaThreadSynchronize() : 0


typedef enum {VEC_SINGLE,VEC_PERSIST,VEC_DEALLOC,VEC_COLLECT} VecUsageGPUFlag;
typedef enum {VEC_UNALLOC,VEC_ALLOC,VEC_GPU, VEC_CPU,VEC_SYNCHED} VecGPUFlag;
typedef struct{
  VECHEADER
    //VecUsageGPUFlag     lifetime;
  VecGPUFlag          syncState;
  //PetscBool           dimsetflag;
  PetscInt*           length;
  PetscScalar*        scalar;
  PetscScalar*        cpuptr;
  PetscScalar*        devptr;
  PetscScalar*        zval;
  cudaStream_t        stream;
  PetscInt*           offset;
  PetscInt*           segment;
  //struct copyCounters vstat;
 }Vec_SeqGPU;




#undef __FUNCT__
#define __FUNCT__ "VecGPUAllocateCheck"
PETSC_STATIC_INLINE PetscErrorCode VecGPUAllocateCheck(Vec v)
{
  //PetscErrorCode ierr;
  //Vec_SeqGPU        *s = (Vec_SeqGPU*)v->data;;

  PetscFunctionBegin;
  printf("Call to VecGPUAllocateCheck_SeqGPU\n");
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "VecCopyToGPU"
/* Copies a vector from the CPU to the GPU unless we already have an up-to-date copy on the GPU */
PETSC_STATIC_INLINE PetscErrorCode VecCopyToGPU(Vec v)
{
  //PetscBLASInt   cn = v->map->n;
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  printf("Call to VecCopyToGPU\n");
  // if (v->valid_GPU_array == PETSC_CUSP_CPU){
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
  //    v->valid_GPU_array = PETSC_CUSP_SYNCHED;
  // }
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VecCopyToGPUSome"
PETSC_STATIC_INLINE PetscErrorCode VecCopyToGPUSome(Vec v,PetscInt x)
{
  //Vec_Seq        *s = (Vec_Seq *)v->data;
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  printf("Call to VecCopyToGPUSome\n");
  //  if (v->valid_GPU_array == PETSC_CUSP_CPU) {
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
  //}
//v->valid_GPU_array = PETSC_CUSP_GPU;
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VecGetArrayReadWrite"
PETSC_STATIC_INLINE PetscErrorCode VecGPUGetArrayReadWrite(Vec v, PetscScalar *va)
{
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  printf("Call to VecGetArrayReadWrite\n");
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VecGPURestoreArrayReadWrite"
PETSC_STATIC_INLINE PetscErrorCode VecGPURestoreArrayReadWrite(Vec v,  PetscScalar *va)
{
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  printf("Call to VecGetRestoreReadWrite\n");
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VecGPUGetArrayRead"
PETSC_STATIC_INLINE PetscErrorCode VecGPUGetArrayRead(Vec v, PetscScalar *va)
{
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  printf("Call to VecGetArrayRead\n");
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VecGPURestoreArrayRead"
PETSC_STATIC_INLINE PetscErrorCode VecGPURestoreArrayRead(Vec v,  PetscScalar *va)
{
  PetscFunctionBegin;
  printf("Call to VecGetRestoreRead\n");
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VecGPUGetArrayWrite"
PETSC_STATIC_INLINE PetscErrorCode VecGPUGetArrayWrite(Vec v,PetscScalar *va)
{
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  printf("Call to VecGetArrayWrite\n");
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "VecGPURestoreArrayWrite"
PETSC_STATIC_INLINE PetscErrorCode VecGPURestoreArrayWrite(Vec v,PetscScalar *va)
{
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  printf("Call to VecRestoreArrayWrite\n");
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  PetscFunctionReturn(0);
}

EXTERN_C_END
#endif



