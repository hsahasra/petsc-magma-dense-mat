
#ifndef __GPUVECIMPL
#define __GPUVECIMPL


/* #include <petscsys.h>
#include <petscvec.h>
 #include <petscconf.h> */

#include <petsc-private/vecimpl.h>
#include <../src/vec/vec/impls/dvecimpl.h>

#define MyPetscStackCheckByName(str,flg)                                 \
  do{ \
    if(PetscStackActive){                                                    \
       flg=0;                                                          \
       int ss = petscstack->currentsize;                               \
       int i=0;                                                        \
       while(i<ss && ss>0){                                           \
         if(strcmp(petscstack->function[i],str) == 0){              \
           flg = 1;                            \
           break;                              \
         }                                     \
         i++;                                  \
       }                                       \
 }else{                                        \
  printf(">>>>>> No petscstack present!!!\n"); \
 }                                             \
 }while(0)

EXTERN_C_BEGIN
#include <sys/time.h>
#include <math.h>
#include <cuda.h>

/*
 These macros define debugging ane output needed for debugging
 VTIMER only has implementation in WAXPY, VecMaxPWDivide,
 NORM2 and DOT
*/
#undef  DEBUGVEC
#define DEBUGVEC        0
#undef  VVERBOSE
#define VVERBOSE        0
#undef  VTIMER
#define VTIMER          0
#define ITSHOW          4000 /* reduces the number of printfs */

/*
  The macros below toggles the manual or Orio code for the
  kernels DOT, NORM2, WAXPY, and AXPY
  0 toggles Orio kernels, 1 toggles manual tuned kernels
*/
#define VMANDOT         1
#define VMANNRM         1
#define VMANWXP         1
#define VMANXPY         1


/*
  These macros below control the tuning of the manual kernels
*/
#define CHUNKWIDTH      65536.0 /* active threads */
#define TCOUNT          128

/* Hardware device maximums */
#define MAXTHREADS      1024
#define MAXBLOCKS       65536

/* Dot product */
#define DOTMPLIER       1.0 /* element-thread work multiplier */
#define THRDOTCNT       256 /* first reduction */
#define THRDOTCNT2      256 /* second reduction */

/* AXPY, WAXPY, MAXPY, AXPBYPCZ */
#define AXPYTCOUNT      256
#define AXPBYPCZTCOUNT  128

/* Norm2 and Infinity Norm */
#define NRMMPLIER       8.0 /* element-thread work multiplier */
#define THRNRMCNT       256 /* first reduction */
#define THRNRMCNT2      256 /* second reduction */

/* Max pointwise divide */
#define MAXMPLIER       1.0 /* element-thread work multiplier */
#define PDIVTCOUNT      256 /* first reduction */
#define PDIVTCOUNT2     256 /* second reduction */

/* Pointwise multiply */
#define PMULTCOUNT      128

/* Copy over (unused) */
#define CPYTCOUNT       128
/* --------------------- */


/* Random number generator macros. Newer SDKs (4.1+)have cuRAND
   libraries, which implement these routines, so this may be obsolete
*/
#define NNN		624
#define MMM		397
#define INIT_MULT	1812433253 /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
#define	ARRAY_SEED	19650218   /* Seed for initial setup before incorp array seed */
#define MATRIX_A	0x9908b0df /* Constant vector a */
#define UPPER_MASK	0x80000000 /* Most significant w-r bits */
#define LOWER_MASK	0x7fffffff /* Least significant r bits */
#define	TEMPER1		0x9d2c5680
#define	TEMPER2		0xefc60000

__shared__ int	mtNexts;	   /* Start of next block of seeds */
__shared__ uint	mtNexti;	   /* Indirect on above to save one global read time/call */
__device__ uint	s_seeds[NNN];
/* --------------------------------------------------------------------- */



PetscErrorCode VecView_SeqGPU_ASCII(Vec ,PetscViewer);

/* NVidia device mallocs specific for pinned memory, required for streaming */
PetscErrorCode PinnedMalloc(PetscScalar** x,PetscInt n);
PetscErrorCode PinnedFree(PetscScalar* x);

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
extern PetscErrorCode VecDuplicateVecs_SeqGPU(Vec,PetscInt,Vec **);
extern PetscErrorCode VecNorm_SeqGPU(Vec,NormType,PetscReal*);
extern PetscErrorCode VecSwap_SeqGPU(Vec,Vec);
extern PetscErrorCode VecCreate_SeqGPU(Vec);

extern PetscErrorCode VecView_SeqGPU(Vec,PetscViewer);
extern PetscErrorCode VecDestroy_SeqGPU(Vec);
extern PetscErrorCode VecDestroyVecs_SeqGPU(PetscInt,Vec*);
extern PetscErrorCode VecAYPX_SeqGPU(Vec,PetscScalar,Vec);
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
extern PetscErrorCode VecAXPY(Vec,PetscReal,Vec);
extern PetscErrorCode VecMax_SeqGPU(Vec,PetscInt*,PetscScalar*);
extern PetscErrorCode VecMin_SeqGPU(Vec,PetscInt*,PetscScalar*);

extern __global__ void kernAXPBYPCZ(double*, double*, double*,int*);
extern __global__ void kernCheck(double*, int*);
extern __global__ void kernCopyLen(int*, int*);
extern __global__ void kernCODevice(double*,double*, int*);
extern __global__ void kernCompare(double*, double*,int*, int*, int*);
extern __global__ void kernSet(double*,double, int);
extern __global__ void kernAddValues(double* x, int n, int* xi, int ni,double *y);
extern __global__ void kernScale(double*,double, int);
extern __global__ void kernCopy(double*, double*, int*, int*);
extern __global__ void kernSwap(int,double*,double*,double*);
extern __global__ void kernAXPY(double*, double*,double, int);
extern __global__ void kernWYMX(double*, double*, int, double*);
extern __global__ void kernXPY(double*, double*, int);



extern __global__ void kernPDIV(double*,double*,int*,double* );
extern __global__ void kernPMULT(double*,double*,int*,double* );
extern __global__ void kernMAXPDIV(double*,double*,int,int,int,double*);
extern __global__ void kernMAX(int,double*);
extern __global__ void kernRand(double*, int*);
extern __global__ void kernRandS(uint *);


__device__ void orcu_warpReduce32(int, volatile double*);
__device__ void orcu_warpReduce64(int, volatile double*);


#if(VMANDOT)
extern __global__ void kernDot(double*,double*,int,int,int,double*);
extern __global__ void kernRedDot(int,double*);
#else
extern __global__ void orcu_dotkernel_1e5(int, double*, double*, double*);
extern __global__ void orcu_dotblksum_1e5(int, double*);
extern __global__ void orcu_dotkernel_1e6(int, double*, double*, double*);
extern __global__ void orcu_dotblksum_1e6(int, double*);
extern __global__ void orcu_dotkernel_1e7(int, double*, double*, double*);
extern __global__ void orcu_dotblksum_1e7(int, double*);
#endif

extern __global__ void kernInfNorm(double*,int,int,int,double*);
extern __global__ void kernRedInfNorm(int,double*);
#if(VMANNRM)
extern __global__ void kernNorm2(double*,int,int,int,double*);
extern __global__ void kernRedNorm(int,double*);

#else
extern __global__ void orcu_norm2kernel_1e5(const int , double*, double*);
extern __global__ void orcu_norm2blksum_1e5(int, double*);
extern __global__ void orcu_norm2kernel_1e6(const int, double*, double*);
extern __global__ void orcu_norm2blksum_1e6(int, double*);
extern __global__ void orcu_norm2kernel_1e7(const int, double*, double*);
extern __global__ void orcu_norm2blksum_1e7(int, double*);
#endif

#if(VMANXPY)

#else
__global__ void orcu_axpykernel_1e5(const int , double , double* , double*);
__global__ void orcu_axpykernel_1e6(const int , double , double* , double*);
__global__ void orcu_axpykernel_1e7(const int , double , double* , double*);
#endif

#if(VMANWXP)
extern __global__ void kernWAXPY(double*, double*,double, int, double*);
extern __global__ void kernWXPY(double*, double*, int, double*);
#else
__global__ void orcu_waxpykernel_1e5(const int, double, double*, double*, double*);
__global__ void orcu_waxpykernel_1e6(const int, double, double*, double*, double*);
__global__ void orcu_waxpykernel_1e7(const int, double, double*, double*, double*);
#endif


__device__ void warpReduce(volatile double*, int);
__device__ void warpDotReduce(volatile double*, int);
__device__ void warpMaxReduce(volatile double*, int);

__device__ void mt19937si(uint);
__device__ void mt19937sai(uint*,uint);
__device__ uint mt19937s(void);
__device__ uint mt19937sl(void);


extern PetscBool  synchronizeGPU;


#define CHKERRGPU(err) if (((int)err) != (int)CUBLAS_STATUS_SUCCESS) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"CUDA error %s",cudaGetErrorString(err))

/* #define WaitForGPU() synchronizeGPU ? cudaDeviceSynchronize() : 0 */

/* Flag used to indicate the state of the Vecotr type with regards to the 
   location of the current uptodate data, if the memory has been allocated,
   and if the arrays have been initialized.                               */
typedef enum {VEC_UNALLOC,VEC_ALLOC,VEC_GPU, VEC_CPU,VEC_SYNCHED} VecGPUFlag;
typedef struct{
  VECHEADER
  VecGPUFlag          syncState; /* memory state */
  PetscInt*           length;    /* local length of the vector */
  PetscScalar*        scalar;    /* misc. on device scalar place holder */
  PetscScalar*        cpuptr;    /* pointer to array on device */
  PetscScalar*        devptr;    /* pointer to array on host   */
  PetscScalar*        zval;      /* future usage for keeping reduction result on device */
  cudaStream_t        streamid;  /* streamID created when Vector created */
  PetscInt*           offset;    /* global array offset */
  PetscInt*           segment;   /* something or other */
 }Vec_SeqGPU;






/*
  Below are legacy functions kept around to check compatability
  All of them have been stripped of functionality, and only printf
  to state they are called.
*/


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



