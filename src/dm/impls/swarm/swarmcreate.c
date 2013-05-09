#define PETSCDM_DLL
#include <petsc-private/dmswarmimpl.h>    /*I   "petscdmswarm.h"   I*/
#include <petscdmda.h>

const char *const DMSwarmPlacements[] = {"LATTICE", "RANDOM", "DMSwarmPlacement", "DM_SWARM_PLACEMENT_", 0};

#undef __FUNCT__
#define __FUNCT__ "DMSetFromOptions_Swarm"
PetscErrorCode  DMSetFromOptions_Swarm(DM dm)
{
  DM_Swarm      *swarm = (DM_Swarm *) dm->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  ierr = PetscOptionsHead("DMSwarm Options");CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-dm_swarm_point_placement", "Method for laying out points, e.g. lattice, random", "DMSwarmSetPlacement", DMSwarmPlacements, (PetscEnum) swarm->pointPlacement, (PetscEnum *) &swarm->pointPlacement, NULL);CHKERRQ(ierr);
  /* Handle DMSwarm refinement */
  /* Handle viewing */
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* External function declarations here */
extern PetscErrorCode DMCreateCoordinateDM_Swarm(DM dm, DM *cdm);
extern PetscErrorCode DMSetUp_Swarm(DM dm);
extern PetscErrorCode DMDestroy_Swarm(DM dm);
extern PetscErrorCode DMView_Swarm(DM dm, PetscViewer viewer);

#undef __FUNCT__
#define __FUNCT__ "DMCreateGlobalVector_Swarm"
static PetscErrorCode DMCreateGlobalVector_Swarm(DM dm,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreateGlobalVector_Section_Private(dm,vec);CHKERRQ(ierr);
  /* ierr = VecSetOperation(*vec, VECOP_DUPLICATE, (void(*)(void)) VecDuplicate_MPI_DM);CHKERRQ(ierr); */
  /* ierr = VecSetOperation(*vec, VECOP_VIEW, (void (*)(void))VecView_Swarm);CHKERRQ(ierr); */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateLocalVector_Swarm"
static PetscErrorCode DMCreateLocalVector_Swarm(DM dm,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreateLocalVector_Section_Private(dm,vec);CHKERRQ(ierr);
  /* ierr = VecSetOperation(*vec, VECOP_VIEW, (void(*)(void)) VecView_Swarm_Local);CHKERRQ(ierr); */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMInitialize_Swarm"
PetscErrorCode DMInitialize_Swarm(DM dm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscStrallocpy(VECSTANDARD, (char**)&dm->vectype);CHKERRQ(ierr);

  dm->ops->view                            = DMView_Swarm;
  dm->ops->setfromoptions                  = DMSetFromOptions_Swarm;
  dm->ops->setup                           = DMSetUp_Swarm;
  dm->ops->createglobalvector              = DMCreateGlobalVector_Swarm;
  dm->ops->createlocalvector               = DMCreateLocalVector_Swarm;
  dm->ops->createlocaltoglobalmapping      = NULL;
  dm->ops->createlocaltoglobalmappingblock = NULL;
  dm->ops->createfieldis                   = NULL;
  dm->ops->createcoordinatedm              = DMCreateCoordinateDM_Swarm;
  dm->ops->getcoloring                     = 0;
  dm->ops->creatematrix                    = NULL;
  dm->ops->createinterpolation             = 0;
  dm->ops->getaggregates                   = 0;
  dm->ops->getinjection                    = 0;
  dm->ops->refine                          = NULL;
  dm->ops->coarsen                         = 0;
  dm->ops->refinehierarchy                 = 0;
  dm->ops->coarsenhierarchy                = 0;
  dm->ops->globaltolocalbegin              = NULL;
  dm->ops->globaltolocalend                = NULL;
  dm->ops->localtoglobalbegin              = NULL;
  dm->ops->localtoglobalend                = NULL;
  dm->ops->destroy                         = DMDestroy_Swarm;
  dm->ops->createsubdm                     = NULL;
  dm->ops->locatepoints                    = NULL;
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "DMCreate_Swarm"
PetscErrorCode DMCreate_Swarm(DM dm)
{
  DM_Swarm      *swarm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  ierr     = PetscNewLog(dm, DM_Swarm, &swarm);CHKERRQ(ierr);
  dm->data = swarm;

  swarm->refct = 1;
  swarm->db    = NULL;
  swarm->vdm   = NULL;
  swarm->pointPlacement = DM_SWARM_PLACEMENT_LATTICE;

  ierr = DMInitialize_Swarm(dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "DMSwarmCreate"
/*@
  DMSwarmCreate - Creates a DMSwarm object, which encapsulates a set of particles. These are often used as Lagrangian traces and in the Material Point Method.

  Collective on MPI_Comm

  Input Parameter:
. comm - The communicator for the DMSwarm object

  Output Parameter:
. mesh  - The DMSwarm object

  Level: beginner

.keywords: DMSwarm, create
@*/
PetscErrorCode DMSwarmCreate(MPI_Comm comm, DM *swamr)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(swamr,2);
  ierr = DMCreate(comm, swamr);CHKERRQ(ierr);
  ierr = DMSetType(*swamr, DMSWARM);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMSwarmClone"
/*@
  DMSwarmClone - Creates a DMSwarm object with the same particles as the original.

  Collective on MPI_Comm

  Input Parameter:
. dm - The original DMSwarm object

  Output Parameter:
. newdm  - The new DMSwarm object

  Level: beginner

.keywords: DMSwarm, create
@*/
PetscErrorCode DMSwarmClone(DM dm, DM *newdm)
{
  DM_Swarm      *swarm;
  void          *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidPointer(newdm,2);
  ierr         = DMCreate(PetscObjectComm((PetscObject)dm), newdm);CHKERRQ(ierr);
  ierr         = PetscSFDestroy(&(*newdm)->sf);CHKERRQ(ierr);
  ierr         = PetscObjectReference((PetscObject) dm->sf);CHKERRQ(ierr);
  (*newdm)->sf = dm->sf;
  swarm         = (DM_Swarm*) dm->data;
  swarm->refct++;
  (*newdm)->data = swarm;
  ierr           = PetscObjectChangeTypeName((PetscObject) *newdm, DMSWARM);CHKERRQ(ierr);
  ierr           = DMInitialize_Swarm(*newdm);CHKERRQ(ierr);
  ierr           = DMGetApplicationContext(dm, &ctx);CHKERRQ(ierr);
  ierr           = DMSetApplicationContext(*newdm, ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
