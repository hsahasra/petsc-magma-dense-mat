/* 
    Defines the vector component of PETSc. Vectors generally represent 
  degrees of freedom for finite element/finite difference functions
  on a grid. They have more mathematical structure then simple arrays.
*/

#ifndef __PETSCVEC_H 
#define __PETSCVEC_H
#include "petscis.h"

PETSC_EXTERN_CXX_BEGIN

/*S
     Vec - Abstract PETSc vector object

   Level: beginner

  Concepts: field variables, unknowns, arrays

.seealso:  VecCreate(), VecType, VecSetType()
S*/
typedef struct _p_Vec*         Vec;

/*S
     VecScatter - Object used to manage communication of data
       between vectors in parallel. Manages both scatters and gathers

   Level: beginner

  Concepts: scatter

.seealso:  VecScatterCreate(), VecScatterBegin(), VecScatterEnd()
S*/
typedef struct _p_VecScatter*  VecScatter;

/*E
  ScatterMode - Determines the direction of a scatter

  Level: beginner

.seealso: VecScatter, VecScatterBegin(), VecScatterEnd()
E*/
typedef enum {SCATTER_FORWARD=0, SCATTER_REVERSE=1, SCATTER_FORWARD_LOCAL=2, SCATTER_REVERSE_LOCAL=3, SCATTER_LOCAL=2} ScatterMode;

/*MC
    SCATTER_FORWARD - Scatters the values as dictated by the VecScatterCreate() call

    Level: beginner

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_REVERSE, SCATTER_FORWARD_LOCAL,
          SCATTER_REVERSE_LOCAL

M*/

/*MC
    SCATTER_REVERSE - Moves the values in the opposite direction then the directions indicated in
         in the VecScatterCreate()

    Level: beginner

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_FORWARD, SCATTER_FORWARD_LOCAL,
          SCATTER_REVERSE_LOCAL

M*/

/*MC
    SCATTER_FORWARD_LOCAL - Scatters the values as dictated by the VecScatterCreate() call except NO parallel communication
       is done. Any variables that have be moved between processes are ignored

    Level: developer

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_REVERSE, SCATTER_FORWARD,
          SCATTER_REVERSE_LOCAL

M*/

/*MC
    SCATTER_REVERSE_LOCAL - Moves the values in the opposite direction then the directions indicated in
         in the VecScatterCreate()  except NO parallel communication
       is done. Any variables that have be moved between processes are ignored

    Level: developer

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_FORWARD, SCATTER_FORWARD_LOCAL,
          SCATTER_REVERSE

M*/

/*J
    VecType - String with the name of a PETSc vector or the creation function
       with an optional dynamic library name, for example
       http://www.mcs.anl.gov/petsc/lib.a:myveccreate()

   Level: beginner

.seealso: VecSetType(), Vec
J*/
#define VecType char*
#define VECSEQ         "seq"
#define VECMPI         "mpi"
#define VECSTANDARD    "standard"   /* seq on one process and mpi on several */
#define VECSHARED      "shared"
#define VECSIEVE       "sieve"
#define VECSEQGPU      "seqgpu"
#define VECSEQCUSP     "seqcusp"
#define VECMPICUSP     "mpicusp"
#define VECCUSP        "cusp"       /* seqcusp on one process and mpicusp on several */
#define VECNEST        "nest"
#define VECSEQPTHREAD  "seqpthread"
#define VECMPIPTHREAD  "mpipthread"
#define VECPTHREAD     "pthread"    /* seqpthread on one process and mpipthread on several */


/* Logging support */
#define    VEC_FILE_CLASSID 1211214
extern  PetscClassId VEC_CLASSID;
extern  PetscClassId VEC_SCATTER_CLASSID;


extern PetscErrorCode  VecInitializePackage(const char[]);
extern PetscErrorCode  VecFinalizePackage(void);

extern PetscErrorCode  VecCreate(MPI_Comm,Vec*);
extern PetscErrorCode  VecCreateSeq(MPI_Comm,PetscInt,Vec*);
extern PetscErrorCode  VecCreateMPI(MPI_Comm,PetscInt,PetscInt,Vec*);
extern PetscErrorCode  VecCreateSeqWithArray(MPI_Comm,PetscInt,PetscInt,const PetscScalar[],Vec*);
extern PetscErrorCode  VecCreateMPIWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscScalar[],Vec*);
extern PetscErrorCode  VecCreateShared(MPI_Comm,PetscInt,PetscInt,Vec*);
extern PetscErrorCode  VecSetFromOptions(Vec);
extern PetscErrorCode  VecSetUp(Vec);
extern PetscErrorCode  VecDestroy(Vec*);
extern PetscErrorCode  VecZeroEntries(Vec);
extern PetscErrorCode  VecSetOptionsPrefix(Vec,const char[]);
extern PetscErrorCode  VecAppendOptionsPrefix(Vec,const char[]);
extern PetscErrorCode  VecGetOptionsPrefix(Vec,const char*[]);

extern PetscErrorCode  VecSetSizes(Vec,PetscInt,PetscInt);

extern PetscErrorCode  VecDotNorm2(Vec,Vec,PetscScalar*,PetscScalar*);
extern PetscErrorCode  VecDot(Vec,Vec,PetscScalar*);
extern PetscErrorCode  VecTDot(Vec,Vec,PetscScalar*);  
extern PetscErrorCode  VecMDot(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode  VecMTDot(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode  VecGetSubVector(Vec,IS,Vec*);
extern PetscErrorCode  VecRestoreSubVector(Vec,IS,Vec*);

/*E
    NormType - determines what type of norm to compute

    Level: beginner

.seealso: VecNorm(), VecNormBegin(), VecNormEnd(), MatNorm()
E*/
typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
extern const char *NormTypes[];
#define NORM_MAX NORM_INFINITY

/*MC
     NORM_1 - the one norm, ||v|| = sum_i | v_i |. ||A|| = max_j || v_*j ||, maximum column sum

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_2, NORM_FROBENIUS, 
           NORM_INFINITY, NORM_1_AND_2

M*/

/*MC
     NORM_2 - the two norm, ||v|| = sqrt(sum_i (v_i)^2) (vectors only)

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_FROBENIUS, 
           NORM_INFINITY, NORM_1_AND_2

M*/

/*MC
     NORM_FROBENIUS - ||A|| = sqrt(sum_ij (A_ij)^2), same as NORM_2 for vectors

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_2, 
           NORM_INFINITY, NORM_1_AND_2

M*/

/*MC
     NORM_INFINITY - ||v|| = max_i |v_i|. ||A|| = max_i || v_i* ||, maximum row sum

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_2, 
           NORM_FROBINIUS, NORM_1_AND_2

M*/

/*MC
     NORM_1_AND_2 - computes both the 1 and 2 norm of a vector

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_2, 
           NORM_FROBINIUS, NORM_INFINITY

M*/

/*MC
     NORM_MAX - see NORM_INFINITY

   Level: beginner

M*/

extern PetscErrorCode  VecNorm(Vec,NormType,PetscReal *);
extern PetscErrorCode  VecNormAvailable(Vec,NormType,PetscBool *,PetscReal *);
extern PetscErrorCode  VecNormalize(Vec,PetscReal *);
extern PetscErrorCode  VecSum(Vec,PetscScalar*);
extern PetscErrorCode  VecMax(Vec,PetscInt*,PetscReal *);
extern PetscErrorCode  VecMin(Vec,PetscInt*,PetscReal *);
extern PetscErrorCode  VecScale(Vec,PetscScalar);
extern PetscErrorCode  VecCopy(Vec,Vec);        
extern PetscErrorCode  VecSetRandom(Vec,PetscRandom);
extern PetscErrorCode  VecSet(Vec,PetscScalar);
extern PetscErrorCode  VecSwap(Vec,Vec);
extern PetscErrorCode  VecAXPY(Vec,PetscScalar,Vec);  
extern PetscErrorCode  VecAXPBY(Vec,PetscScalar,PetscScalar,Vec);  
extern PetscErrorCode  VecMAXPY(Vec,PetscInt,const PetscScalar[],Vec[]);
extern PetscErrorCode  VecAYPX(Vec,PetscScalar,Vec);
extern PetscErrorCode  VecWAXPY(Vec,PetscScalar,Vec,Vec);
extern PetscErrorCode  VecAXPBYPCZ(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
extern PetscErrorCode  VecPointwiseMax(Vec,Vec,Vec);    
extern PetscErrorCode  VecPointwiseMaxAbs(Vec,Vec,Vec);    
extern PetscErrorCode  VecPointwiseMin(Vec,Vec,Vec);    
extern PetscErrorCode  VecPointwiseMult(Vec,Vec,Vec);    
extern PetscErrorCode  VecPointwiseDivide(Vec,Vec,Vec);    
extern PetscErrorCode  VecMaxPointwiseDivide(Vec,Vec,PetscReal*);    
extern PetscErrorCode  VecShift(Vec,PetscScalar);
extern PetscErrorCode  VecReciprocal(Vec);
extern PetscErrorCode  VecPermute(Vec, IS, PetscBool );
extern PetscErrorCode  VecSqrtAbs(Vec);
extern PetscErrorCode  VecLog(Vec);
extern PetscErrorCode  VecExp(Vec);
extern PetscErrorCode  VecAbs(Vec);
extern PetscErrorCode  VecDuplicate(Vec,Vec*);          
extern PetscErrorCode  VecDuplicateVecs(Vec,PetscInt,Vec*[]);         
extern PetscErrorCode  VecDestroyVecs(PetscInt, Vec*[]); 
extern PetscErrorCode  VecStrideNormAll(Vec,NormType,PetscReal[]);
extern PetscErrorCode  VecStrideMaxAll(Vec,PetscInt [],PetscReal []);
extern PetscErrorCode  VecStrideMinAll(Vec,PetscInt [],PetscReal []);
extern PetscErrorCode  VecStrideScaleAll(Vec,const PetscScalar[]);

extern PetscErrorCode  VecStrideNorm(Vec,PetscInt,NormType,PetscReal*);
extern PetscErrorCode  VecStrideMax(Vec,PetscInt,PetscInt *,PetscReal *);
extern PetscErrorCode  VecStrideMin(Vec,PetscInt,PetscInt *,PetscReal *);
extern PetscErrorCode  VecStrideScale(Vec,PetscInt,PetscScalar);
extern PetscErrorCode  VecStrideSet(Vec,PetscInt,PetscScalar);


extern PetscErrorCode  VecStrideGather(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode  VecStrideScatter(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode  VecStrideGatherAll(Vec,Vec[],InsertMode);
extern PetscErrorCode  VecStrideScatterAll(Vec[],Vec,InsertMode);

extern PetscErrorCode  VecSetValues(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode  VecGetValues(Vec,PetscInt,const PetscInt[],PetscScalar[]);
extern PetscErrorCode  VecAssemblyBegin(Vec);
extern PetscErrorCode  VecAssemblyEnd(Vec);
extern PetscErrorCode  VecStashSetInitialSize(Vec,PetscInt,PetscInt);
extern PetscErrorCode  VecStashView(Vec,PetscViewer);
extern PetscErrorCode  VecStashGetInfo(Vec,PetscInt*,PetscInt*,PetscInt*,PetscInt*);

/*MC
   VecSetValue - Set a single entry into a vector.

   Synopsis:
   PetscErrorCode VecSetValue(Vec v,int row,PetscScalar value, InsertMode mode);

   Not Collective

   Input Parameters:
+  v - the vector
.  row - the row location of the entry
.  value - the value to insert
-  mode - either INSERT_VALUES or ADD_VALUES

   Notes:
   For efficiency one should use VecSetValues() and set several or 
   many values simultaneously if possible.

   These values may be cached, so VecAssemblyBegin() and VecAssemblyEnd() 
   MUST be called after all calls to VecSetValues() have been completed.

   VecSetValues() uses 0-based indices in Fortran as well as in C.

   Level: beginner

.seealso: VecSetValues(), VecAssemblyBegin(), VecAssemblyEnd(), VecSetValuesBlockedLocal(), VecSetValueLocal()
M*/
PETSC_STATIC_INLINE PetscErrorCode VecSetValue(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValues(v,1,&i,&va,mode);}


extern PetscErrorCode  VecSetBlockSize(Vec,PetscInt);
extern PetscErrorCode  VecGetBlockSize(Vec,PetscInt*);
extern PetscErrorCode  VecSetValuesBlocked(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

/* Dynamic creation and loading functions */
extern PetscFList VecList;
extern PetscBool  VecRegisterAllCalled;
extern PetscErrorCode  VecSetType(Vec, const VecType);
extern PetscErrorCode  VecGetType(Vec, const VecType *);
extern PetscErrorCode  VecRegister(const char[],const char[],const char[],PetscErrorCode (*)(Vec));
extern PetscErrorCode  VecRegisterAll(const char []);
extern PetscErrorCode  VecRegisterDestroy(void);

/*MC
  VecRegisterDynamic - Adds a new vector component implementation

  Synopsis:
  PetscErrorCode VecRegisterDynamic(const char *name, const char *path, const char *func_name, PetscErrorCode (*create_func)(Vec))

  Not Collective

  Input Parameters:
+ name        - The name of a new user-defined creation routine
. path        - The path (either absolute or relative) of the library containing this routine
. func_name   - The name of routine to create method context
- create_func - The creation routine itself

  Notes:
  VecRegisterDynamic() may be called multiple times to add several user-defined vectors

  If dynamic libraries are used, then the fourth input argument (routine_create) is ignored.

  Sample usage:
.vb
    VecRegisterDynamic("my_vec","/home/username/my_lib/lib/libO/solaris/libmy.a", "MyVectorCreate", MyVectorCreate);
.ve

  Then, your vector type can be chosen with the procedural interface via
.vb
    VecCreate(MPI_Comm, Vec *);
    VecSetType(Vec,"my_vector_name");
.ve
   or at runtime via the option
.vb
    -vec_type my_vector_name
.ve

  Notes: $PETSC_ARCH occuring in pathname will be replaced with appropriate values.
         If your function is not being put into a shared library then use VecRegister() instead
        
  Level: advanced

.keywords: Vec, register
.seealso: VecRegisterAll(), VecRegisterDestroy(), VecRegister()
M*/
#if defined(PETSC_USE_DYNAMIC_LIBRARIES)
#define VecRegisterDynamic(a,b,c,d) VecRegister(a,b,c,0)
#else
#define VecRegisterDynamic(a,b,c,d) VecRegister(a,b,c,d)
#endif


extern PetscErrorCode  VecScatterCreate(Vec,IS,Vec,IS,VecScatter *);
extern PetscErrorCode  VecScatterCreateEmpty(MPI_Comm,VecScatter *);
extern PetscErrorCode  VecScatterCreateLocal(VecScatter,PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],PetscInt);
extern PetscErrorCode  VecScatterBegin(VecScatter,Vec,Vec,InsertMode,ScatterMode);
extern PetscErrorCode  VecScatterEnd(VecScatter,Vec,Vec,InsertMode,ScatterMode); 
extern PetscErrorCode  VecScatterDestroy(VecScatter*);
extern PetscErrorCode  VecScatterCopy(VecScatter,VecScatter *);
extern PetscErrorCode  VecScatterView(VecScatter,PetscViewer);
extern PetscErrorCode  VecScatterRemap(VecScatter,PetscInt *,PetscInt*);
extern PetscErrorCode  VecScatterGetMerged(VecScatter,PetscBool *);

extern PetscErrorCode  VecGetArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode  VecRestoreArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode  VecGetArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode  VecRestoreArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode  VecGetArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode  VecRestoreArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode  VecGetArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);
extern PetscErrorCode  VecRestoreArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);

extern PetscErrorCode  VecPlaceArray(Vec,const PetscScalar[]);
extern PetscErrorCode  VecResetArray(Vec);
extern PetscErrorCode  VecReplaceArray(Vec,const PetscScalar[]);
extern PetscErrorCode  VecGetArrays(const Vec[],PetscInt,PetscScalar**[]);
extern PetscErrorCode  VecRestoreArrays(const Vec[],PetscInt,PetscScalar**[]);

extern PetscErrorCode  VecView(Vec,PetscViewer);
extern PetscErrorCode  VecViewFromOptions(Vec, const char *);
extern PetscErrorCode  VecEqual(Vec,Vec,PetscBool *);
extern PetscErrorCode  VecLoad(Vec, PetscViewer);

extern PetscErrorCode  VecGetSize(Vec,PetscInt*);
extern PetscErrorCode  VecGetLocalSize(Vec,PetscInt*);
extern PetscErrorCode  VecGetOwnershipRange(Vec,PetscInt*,PetscInt*);
extern PetscErrorCode  VecGetOwnershipRanges(Vec,const PetscInt *[]);

extern PetscErrorCode  VecSetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping);
extern PetscErrorCode  VecSetValuesLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

/*MC
   VecSetValueLocal - Set a single entry into a vector using the local numbering

   Synopsis:
   PetscErrorCode VecSetValueLocal(Vec v,int row,PetscScalar value, InsertMode mode);

   Not Collective

   Input Parameters:
+  v - the vector
.  row - the row location of the entry
.  value - the value to insert
-  mode - either INSERT_VALUES or ADD_VALUES

   Notes:
   For efficiency one should use VecSetValues() and set several or 
   many values simultaneously if possible.

   These values may be cached, so VecAssemblyBegin() and VecAssemblyEnd() 
   MUST be called after all calls to VecSetValues() have been completed.

   VecSetValues() uses 0-based indices in Fortran as well as in C.

   Level: beginner

.seealso: VecSetValues(), VecAssemblyBegin(), VecAssemblyEnd(), VecSetValuesBlockedLocal(), VecSetValue()
M*/
PETSC_STATIC_INLINE PetscErrorCode VecSetValueLocal(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValuesLocal(v,1,&i,&va,mode);}

extern PetscErrorCode  VecSetLocalToGlobalMappingBlock(Vec,ISLocalToGlobalMapping);
extern PetscErrorCode  VecSetValuesBlockedLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode  VecGetLocalToGlobalMappingBlock(Vec,ISLocalToGlobalMapping*);
extern PetscErrorCode  VecGetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping*);

extern PetscErrorCode  VecDotBegin(Vec,Vec,PetscScalar *);
extern PetscErrorCode  VecDotEnd(Vec,Vec,PetscScalar *);
extern PetscErrorCode  VecTDotBegin(Vec,Vec,PetscScalar *);
extern PetscErrorCode  VecTDotEnd(Vec,Vec,PetscScalar *);
extern PetscErrorCode  VecNormBegin(Vec,NormType,PetscReal *);
extern PetscErrorCode  VecNormEnd(Vec,NormType,PetscReal *);

extern PetscErrorCode  VecMDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode  VecMDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode  VecMTDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode  VecMTDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode PetscCommSplitReductionBegin(MPI_Comm);


#if defined(PETSC_USE_DEBUG)
#define VecValidValues(vec,argnum,input) do {                           \
    PetscErrorCode     _ierr;                                           \
    PetscInt          _n,_i;                                            \
    const PetscScalar *_x;                                              \
                                                                        \
    if (vec->petscnative || vec->ops->getarray) {                       \
      _ierr = VecGetLocalSize(vec,&_n);CHKERRQ(_ierr);                  \
      _ierr = VecGetArrayRead(vec,&_x);CHKERRQ(_ierr);                  \
      for (_i=0; _i<_n; _i++) {                                         \
        if (input) {                                                    \
          if (PetscIsInfOrNanScalar(_x[_i])) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_FP,"Vec entry at local location %D is not-a-number or infinite at beginning of function: Parameter number %d",_i,argnum); \
        } else {                                                        \
          if (PetscIsInfOrNanScalar(_x[_i])) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_FP,"Vec entry at local location %D is not-a-number or infinite at end of function: Parameter number %d",_i,argnum); \
        }                                                               \
      }                                                                 \
      _ierr = VecRestoreArrayRead(vec,&_x);CHKERRQ(_ierr);              \
    }                                                                   \
  } while (0)
#else
#define VecValidValues(vec,argnum,input)
#endif


typedef enum {VEC_IGNORE_OFF_PROC_ENTRIES,VEC_IGNORE_NEGATIVE_INDICES} VecOption;
extern PetscErrorCode  VecSetOption(Vec,VecOption,PetscBool );

/*
   Expose VecGetArray()/VecRestoreArray() to users. Allows this to work without any function
   call overhead on any 'native' Vecs.
*/

#include "petsc-private/vecimpl.h"

extern PetscErrorCode  VecContourScale(Vec,PetscReal,PetscReal);

/*
    These numbers need to match the entries in 
  the function table in vecimpl.h
*/
typedef enum { VECOP_VIEW = 33, VECOP_LOAD = 41, VECOP_DUPLICATE = 0} VecOperation;
extern PetscErrorCode  VecSetOperation(Vec,VecOperation,void(*)(void));

/*
     Routines for dealing with ghosted vectors:
  vectors with ghost elements at the end of the array.
*/
extern PetscErrorCode  VecMPISetGhost(Vec,PetscInt,const PetscInt[]);
extern PetscErrorCode  VecCreateGhost(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);  
extern PetscErrorCode  VecCreateGhostWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);  
extern PetscErrorCode  VecCreateGhostBlock(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);  
extern PetscErrorCode  VecCreateGhostBlockWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);  
extern PetscErrorCode  VecGhostGetLocalForm(Vec,Vec*);
extern PetscErrorCode  VecGhostRestoreLocalForm(Vec,Vec*);
extern PetscErrorCode  VecGhostUpdateBegin(Vec,InsertMode,ScatterMode);
extern PetscErrorCode  VecGhostUpdateEnd(Vec,InsertMode,ScatterMode);

extern PetscErrorCode  VecConjugate(Vec);

extern PetscErrorCode  VecScatterCreateToAll(Vec,VecScatter*,Vec*);
extern PetscErrorCode  VecScatterCreateToZero(Vec,VecScatter*,Vec*);

extern PetscErrorCode  PetscViewerMathematicaGetVector(PetscViewer, Vec);
extern PetscErrorCode  PetscViewerMathematicaPutVector(PetscViewer, Vec);

/*S
     Vecs - Collection of vectors where the data for the vectors is stored in 
            one contiguous memory

   Level: advanced

   Notes:
    Temporary construct for handling multiply right hand side solves

    This is faked by storing a single vector that has enough array space for 
    n vectors

  Concepts: parallel decomposition

S*/
        struct _n_Vecs  {PetscInt n; Vec v;};
typedef struct _n_Vecs* Vecs;
#define VecsDestroy(x)            (VecDestroy(&(x)->v)         || PetscFree(x))
#define VecsCreateSeq(comm,p,m,x) (PetscNew(struct _n_Vecs,x) || VecCreateSeq(comm,p*m,&(*(x))->v) || (-1 == ((*(x))->n = (m))))
#define VecsCreateSeqWithArray(comm,p,m,a,x) (PetscNew(struct _n_Vecs,x) || VecCreateSeqWithArray(comm,1,p*m,a,&(*(x))->v) || (-1 == ((*(x))->n = (m))))
#define VecsDuplicate(x,y)        (PetscNew(struct _n_Vecs,y) || VecDuplicate(x->v,&(*(y))->v) || (-1 == ((*(y))->n = (x)->n)))



#if defined(PETSC_HAVE_CUSP)
typedef struct _p_PetscCUSPIndices* PetscCUSPIndices;
extern PetscErrorCode PetscCUSPIndicesCreate(PetscInt, PetscInt*,PetscInt, PetscInt*,PetscCUSPIndices*);
extern PetscErrorCode PetscCUSPIndicesDestroy(PetscCUSPIndices*);
extern PetscErrorCode VecCUSPCopyToGPUSome_Public(Vec,PetscCUSPIndices);
extern PetscErrorCode VecCUSPCopyFromGPUSome_Public(Vec,PetscCUSPIndices);

#if defined(PETSC_HAVE_TXPETSCGPU)
extern PetscErrorCode VecCUSPResetIndexBuffersFlagsGPU_Public(PetscCUSPIndices);
extern PetscErrorCode VecCUSPCopySomeToContiguousBufferGPU_Public(Vec,PetscCUSPIndices);
extern PetscErrorCode VecCUSPCopySomeFromContiguousBufferGPU_Public(Vec,PetscCUSPIndices);
extern PetscErrorCode VecScatterInitializeForGPU(VecScatter,Vec,ScatterMode);
extern PetscErrorCode VecScatterFinalizeForGPU(VecScatter);
#endif

extern PetscErrorCode  VecCreateSeqCUSP(MPI_Comm,PetscInt,Vec*);
extern PetscErrorCode  VecCreateMPICUSP(MPI_Comm,PetscInt,PetscInt,Vec*);
#endif

#if defined(PETSC_HAVE_PTHREADCLASSES)
extern PetscErrorCode VecSetNThreads(Vec,PetscInt);
extern PetscErrorCode VecGetNThreads(Vec,PetscInt*);
extern PetscErrorCode VecGetThreadOwnershipRange(Vec,PetscInt,PetscInt*,PetscInt*);
extern PetscErrorCode VecSetThreadAffinities(Vec,const PetscInt[]);
extern PetscErrorCode VecCreateSeqPThread(MPI_Comm,PetscInt,PetscInt,PetscInt[],Vec*);
extern PetscErrorCode VecCreateMPIPThread(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt[],Vec*);
#endif

extern PetscErrorCode  VecNestGetSubVecs(Vec,PetscInt*,Vec**);
extern PetscErrorCode  VecNestGetSubVec(Vec,PetscInt,Vec*);
extern PetscErrorCode  VecNestSetSubVecs(Vec,PetscInt,PetscInt*,Vec*);
extern PetscErrorCode  VecNestSetSubVec(Vec,PetscInt,Vec);
extern PetscErrorCode  VecCreateNest(MPI_Comm,PetscInt,IS*,Vec*,Vec*);
extern PetscErrorCode  VecNestGetSize(Vec,PetscInt*);

PETSC_EXTERN_CXX_END
#endif

