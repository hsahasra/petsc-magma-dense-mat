/* 
   Private context for a Newton line search method for solving
   systems of nonlinear equations
 */

#ifndef __SNES_LS_H
#define __SNES_LS_H
#include <private/snesimpl.h>

typedef struct {
  PetscErrorCode           (*LineSearch)(SNES,void*,Vec,Vec,Vec,PetscReal,PetscReal,Vec,Vec,PetscReal*,PetscReal*,PetscBool *);
  void                     *lsP;                              /* user-defined line-search context (optional) */
  /* --------------- Parameters used by line search method ----------------- */
  PetscReal                alpha;		                                                   /* used to determine sufficient reduction */
  PetscReal                maxstep;                                                          /* maximum step size */
  PetscReal                minlambda;                                                        /* determines smallest line search lambda used */
  PetscErrorCode           (*precheckstep)(SNES,Vec,Vec,void*,PetscBool *);                  /* step-checking routine (optional) */
  void                     *precheck;                                                        /* user-defined step-checking context (optional) */
  PetscErrorCode           (*postcheckstep)(SNES,Vec,Vec,Vec,void*,PetscBool *,PetscBool *); /* step-checking routine (optional) */
  void                     *postcheck;                                                       /* user-defined step-checking context (optional) */
  PetscViewer              monitor;
  PetscReal                damping;
} SNES_LS;

#endif

