#include <petsc-private/dmswarmimpl.h>   /*I      "petscdmswarm.h"   I*/

/* Logging support */
PetscLogEvent DMSWARM_Advect;

#undef __FUNCT__
#define __FUNCT__ "DMView_Swarm"
PetscErrorCode DMView_Swarm(DM dm, PetscViewer viewer)
{
  PetscBool      iascii, isbinary;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidHeaderSpecific(viewer, PETSC_VIEWER_CLASSID, 2);
  ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERASCII, &iascii);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERBINARY, &isbinary);CHKERRQ(ierr);
#if 0
  if (iascii) {
    ierr = DMSwarmView_Ascii(dm, viewer);CHKERRQ(ierr);
  } else if (isbinary) {
    ierr = DMSwarmView_Binary(dm, viewer);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDestroy_Swarm"
PetscErrorCode DMDestroy_Swarm(DM dm)
{
  DM_Swarm      *swarm = (DM_Swarm*) dm->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (--swarm->refct > 0) PetscFunctionReturn(0);
  ierr = DMDestroy(&user->vdm);CHKERRQ(ierr);
  ierr = DataExDestroy(swarm->ex);CHKERRQ(ierr);
  ierr = DataBucketDestroy(&swarm->db);CHKERRQ(ierr);
  /* This was originally freed in DMDestroy(), but that prevents reference counting of backend objects */
  ierr = PetscFree(swarm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMSetUp_Swarm"
PetscErrorCode DMSetUp_Swarm(DM dm)
{
  DM_Swarm      *swarm = (DM_Swarm*) dm->data;
  PetscInt       pointSize = 2*sizeof(PetscScalar), initSize, maxSize;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  /* Create DataBucket */
  ierr = DataBucketCreate(&swarm->db);CHKERRQ(ierr);
  ierr = DataBucketRegisterField(swarm->db, "default", pointSize, PETSC_NULL);CHKERRQ(ierr);
  ierr = DataBucketFinalize(db);CHKERRQ(ierr);
  initSize = 1000; /* Could initialize based upon velocity DM */
  maxSize  = initSize;
  ierr = DataBucketSetInitialSizes(swarm->db, initSize, maxSize);CHKERRQ(ierr);
  /* Set point coordinates */
  if (swarm->pointPlacement == SWARM_PLACEMENT_LATTICE) {
    PetscInt  Nxp[]   = {2,2}; /* change with -lattice_layout_N{x,y,z} */
    PetscReal perturb = 0.1;   /* change with -lattice_layout_perturb */

    ierr = SwarmMPntStd_CoordAssignment_LatticeLayout2d(swarm->vdm, Nxp, perturb, swarm->db);CHKERRQ(ierr);
  } else if (swarm->pointPlacement == SWARM_PLACEMENT_RANDOM) {
    PetscInt nPerCell = 9; /* change with -random_layout_Np */

    ierr = SwarmMPntStd_CoordAssignment_RandomLayout2d(swarm->vdm, nPerCell, swarm->db);CHKERRQ(ierr);
  }
  /* Create the data exchanger needed for parallel particle movement */
  ierr = SwarmDMDA2dDataExchangerCreate(swarm->vdm, &swarm->ex);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateCoordinateDM_Swarm"
PetscErrorCode DMCreateCoordinateDM_Swarm(DM dm, DM *cdm)
{
  PetscSection   section;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMSwarmClone(dm, cdm);CHKERRQ(ierr);
  ierr = PetscSectionCreate(((PetscObject) dm)->comm, &section);CHKERRQ(ierr);
  ierr = DMSetDefaultSection(*cdm, section);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
