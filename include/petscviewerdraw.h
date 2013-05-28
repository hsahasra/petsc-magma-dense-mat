/*
     PetscViewers are objects where other objects can be looked at or stored.
*/

#if !defined(__PETSCVIEWERDRAW_H)
#define __PETSCVIEWERDRAW_H

#include <petscdraw.h>
#include <petscviewer.h>

PETSC_EXTERN PetscErrorCode PetscViewerDrawSetDrawType(PetscViewer,PetscDrawType);

#endif
