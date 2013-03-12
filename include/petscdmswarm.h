/*
  DMSwarm, for parallel particle problems.
*/
#if !defined(__PETSCDMSWARM_H)
#define __PETSCDMSWARM_H

#include <petscdm.h>

/*S
  DMSWARM - DM object that encapsulates a set of particles. These are often used as Lagrangian traces and in the Material Point Method.

  Level: intermediate

  Concepts: particles

.seealso:  DM, DMSwarmCreate()
S*/
PETSC_EXTERN PetscErrorCode DMSwarmCreate(MPI_Comm, DM*);
PETSC_EXTERN PetscErrorCode DMSwarmClone(DM, DM*);

#endif
