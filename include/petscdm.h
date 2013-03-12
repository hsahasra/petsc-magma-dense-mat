/*
      Objects to manage the interactions between the mesh data structures and the algebraic objects
*/
#if !defined(__PETSCDM_H)
#define __PETSCDM_H
#include <petscmat.h>
#include <petscdmtypes.h>

PETSC_EXTERN PetscErrorCode DMInitializePackage(const char[]);

PETSC_EXTERN PetscClassId DM_CLASSID;

/*J
    DMType - String with the name of a PETSc DM or the creation function
       with an optional dynamic library name, for example
       http://www.mcs.anl.gov/petsc/lib.a:mydmcreate()

   Level: beginner

.seealso: DMSetType(), DM
J*/
typedef const char* DMType;
#define DMDA        "da"
#define DMADDA      "adda"
#define DMCOMPOSITE "composite"
#define DMSLICED    "sliced"
#define DMSHELL     "shell"
#define DMMESH      "mesh"
#define DMPLEX      "plex"
#define DMCARTESIAN "cartesian"
#define DMREDUNDANT "redundant"
#define DMAKKT      "akkt"
#define DMPATCH     "patch"
#define DMSWARM     "swarm"

PETSC_EXTERN PetscFunctionList DMList;
PETSC_EXTERN PetscBool         DMRegisterAllCalled;
PETSC_EXTERN PetscErrorCode DMCreate(MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMSetType(DM, DMType);
PETSC_EXTERN PetscErrorCode DMGetType(DM, DMType *);
PETSC_EXTERN PetscErrorCode DMRegister(const char[],const char[],const char[],PetscErrorCode (*)(DM));
PETSC_EXTERN PetscErrorCode DMRegisterAll(const char []);
PETSC_EXTERN PetscErrorCode DMRegisterDestroy(void);


/*MC
  DMRegisterDynamic - Adds a new DM component implementation

  Synopsis:
  #include "petscdm.h"
  PetscErrorCode DMRegisterDynamic(const char *name,const char *path,const char *func_name, PetscErrorCode (*create_func)(DM))

  Not Collective

  Input Parameters:
+ name        - The name of a new user-defined creation routine
. path        - The path (either absolute or relative) of the library containing this routine
. func_name   - The name of routine to create method context
- create_func - The creation routine itself

  Notes:
  DMRegisterDynamic() may be called multiple times to add several user-defined DMs

  If dynamic libraries are used, then the fourth input argument (routine_create) is ignored.

  Sample usage:
.vb
    DMRegisterDynamic("my_da","/home/username/my_lib/lib/libO/solaris/libmy.a", "MyDMCreate", MyDMCreate);
.ve

  Then, your DM type can be chosen with the procedural interface via
.vb
    DMCreate(MPI_Comm, DM *);
    DMSetType(DM,"my_da_name");
.ve
   or at runtime via the option
.vb
    -da_type my_da_name
.ve

  Notes: $PETSC_ARCH occuring in pathname will be replaced with appropriate values.
         If your function is not being put into a shared library then use DMRegister() instead

  Level: advanced

.keywords: DM, register
.seealso: DMRegisterAll(), DMRegisterDestroy(), DMRegister()
M*/
#if defined(PETSC_USE_DYNAMIC_LIBRARIES)
#define DMRegisterDynamic(a,b,c,d) DMRegister(a,b,c,0)
#else
#define DMRegisterDynamic(a,b,c,d) DMRegister(a,b,c,d)
#endif

PETSC_EXTERN PetscErrorCode DMView(DM,PetscViewer);
PETSC_EXTERN PetscErrorCode DMLoad(DM,PetscViewer);
PETSC_EXTERN PetscErrorCode DMDestroy(DM*);
PETSC_EXTERN PetscErrorCode DMCreateGlobalVector(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMCreateLocalVector(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMGetLocalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMRestoreLocalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMGetGlobalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMRestoreGlobalVector(DM,Vec *);
PETSC_EXTERN PetscErrorCode DMClearGlobalVectors(DM);
PETSC_EXTERN PetscErrorCode DMGetNamedGlobalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMRestoreNamedGlobalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMGetNamedLocalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMRestoreNamedLocalVector(DM,const char*,Vec*);
PETSC_EXTERN PetscErrorCode DMGetLocalToGlobalMapping(DM,ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode DMGetLocalToGlobalMappingBlock(DM,ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode DMCreateFieldIS(DM,PetscInt*,char***,IS**);
PETSC_EXTERN PetscErrorCode DMGetBlockSize(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMCreateColoring(DM,ISColoringType,MatType,ISColoring*);
PETSC_EXTERN PetscErrorCode DMCreateMatrix(DM,MatType,Mat*);
PETSC_EXTERN PetscErrorCode DMSetMatrixPreallocateOnly(DM,PetscBool);
PETSC_EXTERN PetscErrorCode DMCreateInterpolation(DM,DM,Mat*,Vec*);
PETSC_EXTERN PetscErrorCode DMCreateInjection(DM,DM,VecScatter*);
PETSC_EXTERN PetscErrorCode DMGetWorkArray(DM,PetscInt,PetscDataType,void*);
PETSC_EXTERN PetscErrorCode DMRestoreWorkArray(DM,PetscInt,PetscDataType,void*);
PETSC_EXTERN PetscErrorCode DMRefine(DM,MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMCoarsen(DM,MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMRefineHierarchy(DM,PetscInt,DM[]);
PETSC_EXTERN PetscErrorCode DMCoarsenHierarchy(DM,PetscInt,DM[]);
PETSC_EXTERN PetscErrorCode DMCoarsenHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,Vec,Mat,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMRefineHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMRestrict(DM,Mat,Vec,Mat,DM);
PETSC_EXTERN PetscErrorCode DMInterpolate(DM,Mat,DM);
PETSC_EXTERN PetscErrorCode DMSetFromOptions(DM);
PETSC_EXTERN PetscErrorCode DMViewFromOptions(DM,const char[]);

PETSC_EXTERN PetscErrorCode DMSetUp(DM);
PETSC_EXTERN PetscErrorCode DMCreateInterpolationScale(DM,DM,Mat,Vec*);
PETSC_EXTERN PetscErrorCode DMCreateAggregates(DM,DM,Mat*);
PETSC_EXTERN PetscErrorCode DMGlobalToLocalHookAdd(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode DMGlobalToLocalBegin(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMGlobalToLocalEnd(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToGlobalBegin(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMLocalToGlobalEnd(DM,Vec,InsertMode,Vec);
PETSC_EXTERN PetscErrorCode DMConvert(DM,DMType,DM*);

PETSC_EXTERN PetscErrorCode DMGetCoordinateDM(DM,DM*);
PETSC_EXTERN PetscErrorCode DMGetCoordinates(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMSetCoordinates(DM,Vec);
PETSC_EXTERN PetscErrorCode DMGetCoordinatesLocal(DM,Vec*);
PETSC_EXTERN PetscErrorCode DMSetCoordinatesLocal(DM,Vec);
PETSC_EXTERN PetscErrorCode DMLocatePoints(DM,Vec,IS*);

/* block hook interface */
PETSC_EXTERN PetscErrorCode DMSubDomainHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,VecScatter,VecScatter,DM,void*),void*);
PETSC_EXTERN PetscErrorCode DMSubDomainRestrict(DM,VecScatter,VecScatter,DM);

PETSC_EXTERN PetscErrorCode DMSetOptionsPrefix(DM,const char []);
PETSC_EXTERN PetscErrorCode DMSetVecType(DM,VecType);
PETSC_EXTERN PetscErrorCode DMSetMatType(DM,MatType);
PETSC_EXTERN PetscErrorCode DMSetApplicationContext(DM,void*);
PETSC_EXTERN PetscErrorCode DMSetApplicationContextDestroy(DM,PetscErrorCode (*)(void**));
PETSC_EXTERN PetscErrorCode DMGetApplicationContext(DM,void*);
PETSC_EXTERN PetscErrorCode DMSetVariableBounds(DM,PetscErrorCode (*)(DM,Vec,Vec));
PETSC_EXTERN PetscErrorCode DMHasVariableBounds(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMHasColoring(DM,PetscBool *);
PETSC_EXTERN PetscErrorCode DMComputeVariableBounds(DM,Vec,Vec);

PETSC_EXTERN PetscErrorCode DMCreateSubDM(DM, PetscInt, PetscInt[], IS *, DM *);
PETSC_EXTERN PetscErrorCode DMCreateFieldDecomposition(DM,PetscInt*,char***,IS**,DM**);
PETSC_EXTERN PetscErrorCode DMCreateDomainDecomposition(DM,PetscInt*,char***,IS**,IS**,DM**);
PETSC_EXTERN PetscErrorCode DMCreateDomainDecompositionScatters(DM,PetscInt,DM*,VecScatter**,VecScatter**,VecScatter**);

PETSC_EXTERN PetscErrorCode DMGetRefineLevel(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMGetCoarsenLevel(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMFinalizePackage(void);

PETSC_EXTERN PetscErrorCode VecGetDM(Vec, DM*);
PETSC_EXTERN PetscErrorCode VecSetDM(Vec, DM);
PETSC_EXTERN PetscErrorCode MatGetDM(Mat, DM*);
PETSC_EXTERN PetscErrorCode MatSetDM(Mat, DM);

typedef struct NLF_DAAD* NLF;

#define DM_FILE_CLASSID 1211221

/* FEM support */
PETSC_EXTERN PetscErrorCode DMPrintCellVector(PetscInt, const char [], PetscInt, const PetscScalar []);
PETSC_EXTERN PetscErrorCode DMPrintCellMatrix(PetscInt, const char [], PetscInt, PetscInt, const PetscScalar []);

PETSC_EXTERN PetscErrorCode DMGetDefaultSection(DM, PetscSection *);
PETSC_EXTERN PetscErrorCode DMSetDefaultSection(DM, PetscSection);
PETSC_EXTERN PetscErrorCode DMGetDefaultGlobalSection(DM, PetscSection *);
PETSC_EXTERN PetscErrorCode DMSetDefaultGlobalSection(DM, PetscSection);
PETSC_EXTERN PetscErrorCode DMGetDefaultSF(DM, PetscSF *);
PETSC_EXTERN PetscErrorCode DMSetDefaultSF(DM, PetscSF);
PETSC_EXTERN PetscErrorCode DMCreateDefaultSF(DM, PetscSection, PetscSection);
PETSC_EXTERN PetscErrorCode DMGetPointSF(DM, PetscSF *);
PETSC_EXTERN PetscErrorCode DMSetPointSF(DM, PetscSF);

PETSC_EXTERN PetscErrorCode DMGetNumFields(DM, PetscInt *);
PETSC_EXTERN PetscErrorCode DMSetNumFields(DM, PetscInt);
PETSC_EXTERN PetscErrorCode DMGetField(DM, PetscInt, PetscObject *);

typedef struct {
  PetscInt         numQuadPoints; /* The number of quadrature points on an element */
  const PetscReal *quadPoints;    /* The quadrature point coordinates */
  const PetscReal *quadWeights;   /* The quadrature weights */
  PetscInt         numBasisFuncs; /* The number of finite element basis functions on an element */
  PetscInt         numComponents; /* The number of components for each basis function */
  const PetscReal *basis;         /* The basis functions tabulated at the quadrature points */
  const PetscReal *basisDer;      /* The basis function derivatives tabulated at the quadrature points */
} PetscQuadrature;

typedef struct {
  PetscQuadrature *quad;
  void (**f0Funcs)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[]); /* The f_0 functions for each field */
  void (**f1Funcs)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[]); /* The f_1 functions for each field */
  void (**g0Funcs)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[]); /* The g_0 functions for each field pair */
  void (**g1Funcs)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[]); /* The g_1 functions for each field pair */
  void (**g2Funcs)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[]); /* The g_2 functions for each field pair */
  void (**g3Funcs)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[]); /* The g_3 functions for each field pair */
  PetscScalar (**bcFuncs)(const PetscReal x[]); /* The boundary condition function for each field component */
} PetscFEM;

typedef enum {PETSC_UNIT_LENGTH, PETSC_UNIT_MASS, PETSC_UNIT_TIME, PETSC_UNIT_CURRENT, PETSC_UNIT_TEMPERATURE, PETSC_UNIT_AMOUNT, PETSC_UNIT_LUMINOSITY, NUM_PETSC_UNITS} PetscUnit;

struct _DMInterpolationInfo {
  MPI_Comm   comm;
  PetscInt   dim;    /*1 The spatial dimension of points */
  PetscInt   nInput; /* The number of input points */
  PetscReal *points; /* The input point coordinates */
  PetscInt  *cells;  /* The cell containing each point */
  PetscInt   n;      /* The number of local points */
  Vec        coords; /* The point coordinates */
  PetscInt   dof;    /* The number of components to interpolate */
};
typedef struct _DMInterpolationInfo *DMInterpolationInfo;

PetscErrorCode DMInterpolationCreate(MPI_Comm, DMInterpolationInfo *);
PetscErrorCode DMInterpolationSetDim(DMInterpolationInfo, PetscInt);
PetscErrorCode DMInterpolationGetDim(DMInterpolationInfo, PetscInt *);
PetscErrorCode DMInterpolationSetDof(DMInterpolationInfo, PetscInt);
PetscErrorCode DMInterpolationGetDof(DMInterpolationInfo, PetscInt *);
PetscErrorCode DMInterpolationAddPoints(DMInterpolationInfo, PetscInt, PetscReal[]);
PetscErrorCode DMInterpolationSetUp(DMInterpolationInfo, DM, PetscBool);
PetscErrorCode DMInterpolationGetCoordinates(DMInterpolationInfo, Vec *);
PetscErrorCode DMInterpolationGetVector(DMInterpolationInfo, Vec *);
PetscErrorCode DMInterpolationRestoreVector(DMInterpolationInfo, Vec *);
PetscErrorCode DMInterpolationEvaluate(DMInterpolationInfo, DM, Vec, Vec);
PetscErrorCode DMInterpolationDestroy(DMInterpolationInfo *);
#endif
