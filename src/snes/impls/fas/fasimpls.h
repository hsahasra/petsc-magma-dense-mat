#ifndef _SNES_FASIMPLS
#define _SNES_FASIMPLS

#include <private/snesimpl.h>
#include <private/dmimpl.h>

typedef struct {

  /* flags for knowing the global place of this FAS object */
  PetscInt       level;                        /* level = 0 coarsest level */
  PetscInt       levels;                       /* if level + 1 = levels; we're the last turtle */


  /* smoothing objects */
  SNES           upsmooth;                     /* the SNES for presmoothing */
  SNES           downsmooth;                   /* the SNES for postsmoothing */

  /* coarse grid correction objects */
  SNES           next;                         /* the SNES instance for the next level in the hierarchy */
  Mat            interpolate;                  /* interpolation */
  Mat            restrct;                      /* restriction operator */
  Vec            rscale;                       /* the pointwise scaling of the restriction operator */

  /* method parameters */
  PetscInt       n_cycles;                     /* number of cycles on this level */
  PetscInt       max_up_it;                    /* number of pre-smooths */
  PetscInt       max_down_it;                  /* number of post-smooth cycles */

} SNES_FAS;

#endif
