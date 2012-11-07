
#ifndef __TCOPENMPIMPLH
#define __TCOPENMPIMPLH

#include <petsc-private/threadcommimpl.h>
EXTERN_C_BEGIN
extern PetscErrorCode PetscThreadCommCreate_OpenMP(PetscThreadComm);
EXTERN_C_END

extern PetscErrorCode PetscThreadCommRunKernel_OpenMP(MPI_Comm,PetscThreadCommJobCtx);

#endif
