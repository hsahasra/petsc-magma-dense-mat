/* Defines the basic SNES object */
#include <../src/snes/impls/fas/fasimpls.h>

/*MC
Full Approximation Scheme nonlinear multigrid solver.

The nonlinear problem is solved via the repeated application of nonlinear preconditioners and coarse-grid corrections

.seealso: SNESCreate(), SNES, SNESSetType(), SNESType (for list of available types)
M*/

extern PetscErrorCode SNESDestroy_FAS(SNES snes);
extern PetscErrorCode SNESSetUp_FAS(SNES snes);
extern PetscErrorCode SNESSetFromOptions_FAS(SNES snes);
extern PetscErrorCode SNESView_FAS(SNES snes, PetscViewer viewer);
extern PetscErrorCode SNESSolve_FAS(SNES snes);
extern PetscErrorCode SNESReset_FAS(SNES snes);

EXTERN_C_BEGIN

#undef __FUNCT__
#define __FUNCT__ "SNESCreate_FAS"
PetscErrorCode SNESCreate_FAS(SNES snes)
{
  SNES_FAS * fas;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  snes->ops->destroy        = SNESDestroy_FAS;
  snes->ops->setup          = SNESSetUp_FAS;
  snes->ops->setfromoptions = SNESSetFromOptions_FAS;
  snes->ops->view           = SNESView_FAS;
  snes->ops->solve          = SNESSolve_FAS;
  snes->ops->reset          = SNESReset_FAS;

  snes->usesksp             = PETSC_FALSE;
  snes->usespc              = PETSC_FALSE;

  ierr = PetscNewLog(snes, SNES_FAS, &fas);CHKERRQ(ierr);
  snes->data                = (void*) fas;
  fas->level                = 0;
  fas->levels               = 1;
  fas->n_cycles             = 1;
  fas->max_up_it            = 1;
  fas->max_down_it          = 1;
  fas->upsmooth             = PETSC_NULL;
  fas->downsmooth           = PETSC_NULL;
  fas->next                 = PETSC_NULL;
  fas->interpolate          = PETSC_NULL;
  fas->restrct              = PETSC_NULL;

  PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "SNESFASGetLevels"
PetscErrorCode SNESFASGetLevels(SNES snes, PetscInt * levels) {
  SNES_FAS * fas = (SNES_FAS *)snes->data;
  PetscFunctionBegin;
  *levels = fas->levels;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESFASGetSNES"
PetscErrorCode SNESFASGetSNES(SNES snes, PetscInt level, SNES * lsnes) {
  SNES_FAS * fas = (SNES_FAS *)snes->data;
  PetscInt levels = fas->level;
  PetscInt i;
  PetscFunctionBegin;
  *lsnes = snes;
  if (fas->level < level) {
    SETERRQ(((PetscObject)snes)->comm, PETSC_ERR_ARG_OUTOFRANGE, "SNESFASGetSNES should only be called on a finer SNESFAS instance than the level.");
  }
  if (level > levels - 1) {
    SETERRQ(((PetscObject)snes)->comm, PETSC_ERR_ARG_OUTOFRANGE, "Level %d doesn't exist in the SNESFAS.");
  }
  for (i = fas->level; i > level; i--) {
    *lsnes = fas->next;
    fas = (SNES_FAS *)(*lsnes)->data;
  }
  if (fas->level != level) SETERRQ(((PetscObject)snes)->comm, PETSC_ERR_ARG_OUTOFRANGE, "SNESFASGetSNES didn't return the right level!");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESFASSetLevels"
PetscErrorCode SNESFASSetLevels(SNES snes, PetscInt levels, MPI_Comm * comms) {
  PetscErrorCode ierr;
  PetscInt i;
  SNES_FAS * fas = (SNES_FAS *)snes->data;
  MPI_Comm comm;
  PetscFunctionBegin;
  comm = ((PetscObject)snes)->comm;
  if (levels == fas->levels) {
    if (!comms)
      PetscFunctionReturn(0);
  }
  /* user has changed the number of levels; reset */
  ierr = SNESReset(snes);CHKERRQ(ierr);
  /* destroy any coarser levels if necessary */
  if (fas->next) SNESDestroy(&fas->next);CHKERRQ(ierr);
  fas->next = PETSC_NULL;
  /* setup the finest level */
  for (i = levels - 1; i >= 0; i--) {
    if (comms) comm = comms[i];
    fas->level = i;
    fas->levels = levels;
    fas->next = PETSC_NULL;
    if (i > 0) {
      ierr = SNESCreate(comm, &fas->next);CHKERRQ(ierr);
      ierr = PetscObjectIncrementTabLevel((PetscObject)fas->next, (PetscObject)snes, levels - i);CHKERRQ(ierr);
      ierr = SNESSetType(fas->next, SNESFAS);CHKERRQ(ierr);
      fas = (SNES_FAS *)fas->next->data;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESFASSetInterpolation"
PetscErrorCode SNESFASSetInterpolation(SNES snes, PetscInt level, Mat mat) {
  SNES_FAS * fas =  (SNES_FAS *)snes->data;
  PetscInt top_level = fas->level,i;

  PetscFunctionBegin;
  if (level > top_level)
    SETERRQ1(((PetscObject)snes)->comm, PETSC_ERR_ARG_OUTOFRANGE, "Bad level number %d in SNESFASSetInterpolation", level);
  /* get to the correct level */
  for (i = fas->level; i > level; i--) {
    fas = (SNES_FAS *)fas->next->data;
  }
  if (fas->level != level)
    SETERRQ(((PetscObject)snes)->comm, PETSC_ERR_ARG_WRONG, "Inconsistent level labelling in SNESFASSetInterpolation");
  fas->interpolate = mat;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESFASSetRestriction"
PetscErrorCode SNESFASSetRestriction(SNES snes, PetscInt level, Mat mat) {
  SNES_FAS * fas =  (SNES_FAS *)snes->data;
  PetscInt top_level = fas->level,i;

  PetscFunctionBegin;
  if (level > top_level)
    SETERRQ1(((PetscObject)snes)->comm, PETSC_ERR_ARG_OUTOFRANGE, "Bad level number %d in SNESFASSetRestriction", level);
  /* get to the correct level */
  for (i = fas->level; i > level; i--) {
    fas = (SNES_FAS *)fas->next->data;
  }
  if (fas->level != level)
    SETERRQ(((PetscObject)snes)->comm, PETSC_ERR_ARG_WRONG, "Inconsistent level labelling in SNESFASSetRestriction");
  fas->restrct = mat;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESFASSetRScale"
PetscErrorCode SNESFASSetRScale(SNES snes, PetscInt level, Vec rscale) {
  SNES_FAS * fas =  (SNES_FAS *)snes->data;
  PetscInt top_level = fas->level,i;

  PetscFunctionBegin;
  if (level > top_level)
    SETERRQ1(((PetscObject)snes)->comm, PETSC_ERR_ARG_OUTOFRANGE, "Bad level number %d in SNESFASSetRestriction", level);
  /* get to the correct level */
  for (i = fas->level; i > level; i--) {
    fas = (SNES_FAS *)fas->next->data;
  }
  if (fas->level != level)
    SETERRQ(((PetscObject)snes)->comm, PETSC_ERR_ARG_WRONG, "Inconsistent level labelling in SNESFASSetRestriction");
  fas->rscale = rscale;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESReset_FAS"
PetscErrorCode SNESReset_FAS(SNES snes)
{
  PetscErrorCode ierr = 0;
  SNES_FAS * fas = (SNES_FAS *)snes->data;

  PetscFunctionBegin;
  /* destroy local data created in SNESSetup_FAS */
#if 0
  /* recurse -- reset should destroy the structures -- destroy should destroy the structures recursively */
#endif
  if (fas->next) ierr = SNESReset_FAS(fas->next);CHKERRQ(ierr);
#if 0
#endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESDestroy_FAS"
PetscErrorCode SNESDestroy_FAS(SNES snes)
{
  SNES_FAS * fas = (SNES_FAS *)snes->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* recursively resets and then destroys */
  ierr = SNESReset_FAS(snes);CHKERRQ(ierr);
  if (fas->upsmooth)     ierr = SNESDestroy(&fas->upsmooth);CHKERRQ(ierr);
  if (fas->downsmooth)   ierr = SNESDestroy(&fas->downsmooth);CHKERRQ(ierr);
  if (fas->interpolate == fas->restrct) {
    if (fas->interpolate)  ierr = MatDestroy(&fas->interpolate);CHKERRQ(ierr);
    fas->restrct = PETSC_NULL;
  } else {
    if (fas->interpolate)  ierr = MatDestroy(&fas->interpolate);CHKERRQ(ierr);
    if (fas->restrct)      ierr = MatDestroy(&fas->restrct);CHKERRQ(ierr);
  }
  if (fas->rscale)       ierr = VecDestroy(&fas->rscale);CHKERRQ(ierr);
  if (snes->work)        ierr = VecDestroyVecs(snes->nwork,&snes->work);CHKERRQ(ierr);
  if (fas->next)         ierr = SNESDestroy(&fas->next);CHKERRQ(ierr);
  ierr = PetscFree(fas);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESSetUp_FAS"
PetscErrorCode SNESSetUp_FAS(SNES snes)
{
  SNES_FAS       *fas = (SNES_FAS *) snes->data,*tmp;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* should call the SNESSetFromOptions() only when approriate */
  tmp = fas;
  while (tmp) {
    if (tmp->upsmooth) {ierr = SNESSetFromOptions(tmp->upsmooth);}
    if (tmp->downsmooth) {ierr = SNESSetFromOptions(tmp->upsmooth);}
    tmp = tmp->next ? (SNES_FAS*) tmp->next->data : 0;
  }

  if (!snes->work || snes->nwork != 2) {ierr = SNESDefaultGetWork(snes, 2);CHKERRQ(ierr);} /* work vectors used for intergrid transfers */
  /* gets the solver ready for solution */
  if (snes->dm) {
    /* construct EVERYTHING from the DM -- including the progressive set of smoothers */
    if (fas->next) {
      /* for now -- assume the DM and the evaluation functions have been set externally */
      if (!fas->next->dm) {
        ierr = DMCoarsen(snes->dm, ((PetscObject)fas->next)->comm, &fas->next->dm);CHKERRQ(ierr);
        ierr = SNESSetDM(fas->next, fas->next->dm);CHKERRQ(ierr);
      }
      /* set the interpolation and restriction from the DM */
      if (!fas->interpolate) {
        ierr = DMGetInterpolation(fas->next->dm, snes->dm, &fas->interpolate, &fas->rscale);CHKERRQ(ierr);
        fas->restrct = fas->interpolate;
      }
    }
    /* set the DMs of the pre and post-smoothers here */
    if (fas->upsmooth)  SNESSetDM(fas->upsmooth,   snes->dm);CHKERRQ(ierr);
    if (fas->downsmooth)SNESSetDM(fas->downsmooth, snes->dm);CHKERRQ(ierr);
  }
  if (fas->next) {
    /* gotta set up the solution vector for this to work */
    if (!fas->next->vec_sol)ierr = VecDuplicate(fas->rscale, &fas->next->vec_sol);CHKERRQ(ierr);
    ierr = SNESSetUp(fas->next);CHKERRQ(ierr);
  }
  /* got to set them all up at once */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESSetFromOptions_FAS"
PetscErrorCode SNESSetFromOptions_FAS(SNES snes)
{
  SNES_FAS   *fas = (SNES_FAS *) snes->data;
  PetscInt levels = 1;
  PetscBool flg, monflg;
  PetscErrorCode ierr;
  const char * def_smooth = SNESNRICHARDSON;
  char pre_type[256];
  char post_type[256];
  char                    monfilename[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  ierr = PetscOptionsHead("SNESFAS Options-----------------------------------");CHKERRQ(ierr);

  /* number of levels -- only process on the finest level */
  if (fas->levels - 1 == fas->level) {
    ierr = PetscOptionsInt("-snes_fas_levels", "Number of Levels", "SNESFASSetLevels", levels, &levels, &flg);CHKERRQ(ierr);
    if (!flg && snes->dm) {
      ierr = DMGetRefineLevel(snes->dm,&levels);CHKERRQ(ierr);
      levels++;
    }
    ierr = SNESFASSetLevels(snes, levels, PETSC_NULL);CHKERRQ(ierr);
  }

  /* type of pre and/or post smoothers -- set both at once */
  ierr = PetscMemcpy(post_type, def_smooth, 256);CHKERRQ(ierr);
  ierr = PetscMemcpy(pre_type, def_smooth, 256);CHKERRQ(ierr);
  ierr = PetscOptionsList("-snes_fas_smoother_type","Nonlinear smoother method","SNESSetType",SNESList,def_smooth,pre_type,256,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscMemcpy(post_type, pre_type, 256);CHKERRQ(ierr);
  } else {
    ierr = PetscOptionsList("-snes_fas_smoothup_type",  "Nonlinear smoother method","SNESSetType",SNESList,def_smooth,pre_type, 256,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsList("-snes_fas_smoothdown_type","Nonlinear smoother method","SNESSetType",SNESList,def_smooth,post_type,256,&flg);CHKERRQ(ierr);
  }

  /* options for the number of preconditioning cycles and cycle type */
  ierr = PetscOptionsInt("-snes_fas_up_it","Number of upsmoother iterations","PCMGSetNumberSmoothUp",fas->max_up_it,&fas->max_up_it,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-snes_fas_down_it","Number of downsmoother iterations","PCMGSetNumberSmoothUp",fas->max_down_it,&fas->max_down_it,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-snes_fas_cycles","Number of cycles","PCMGSetNumberSmoothUp",fas->n_cycles,&fas->n_cycles,&flg);CHKERRQ(ierr);

  ierr = PetscOptionsString("-snes_fas_monitor","Monitor for smoothers","SNESMonitorSet","stdout",monfilename,PETSC_MAX_PATH_LEN,&monflg);CHKERRQ(ierr);

  /* other options for the coarsest level */
  if (fas->level == 0) {
    ierr = PetscOptionsList("-snes_fas_coarse_smoother_type","Coarsest smoother method","SNESSetType",SNESList,def_smooth,pre_type,256,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscMemcpy(post_type, pre_type, 256);CHKERRQ(ierr);
    } else {
      ierr = PetscOptionsList("-snes_fas_coarse_smoothup_type",  "Nonlinear smoother method","SNESSetType",SNESList,def_smooth,pre_type, 256,&flg);CHKERRQ(ierr);
      ierr = PetscOptionsList("-snes_fas_coarse_smoothdown_type","Nonlinear smoother method","SNESSetType",SNESList,def_smooth,post_type,256,&flg);CHKERRQ(ierr);
    }
  }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  /* setup from the determined types if the smoothers don't exist */
  if (!fas->upsmooth) {
    const char     *prefix;
    ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
    ierr = SNESCreate(((PetscObject)snes)->comm, &fas->upsmooth);CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(fas->upsmooth,prefix);CHKERRQ(ierr);
    if (fas->level || (fas->levels == 1)) {
      ierr = SNESAppendOptionsPrefix(fas->upsmooth,"fas_levels_");CHKERRQ(ierr);
    } else {
      ierr = SNESAppendOptionsPrefix(fas->upsmooth,"fas_coarse_");CHKERRQ(ierr);
    }
    ierr = PetscObjectIncrementTabLevel((PetscObject)fas->upsmooth, (PetscObject)snes, 1);CHKERRQ(ierr);
    ierr = SNESSetType(fas->upsmooth, pre_type);CHKERRQ(ierr);
    if (snes->ops->computefunction) {
      ierr = SNESSetFunction(fas->upsmooth,PETSC_NULL,snes->ops->computefunction,snes->funP);CHKERRQ(ierr);
    }
  }
  if (fas->upsmooth) {
    ierr = SNESSetTolerances(fas->upsmooth, 0.0, 0.0, 0.0, fas->max_up_it, 1000);CHKERRQ(ierr);
  }

  if (!fas->downsmooth && fas->level != 0) {
    const char     *prefix;
    ierr = SNESGetOptionsPrefix(snes,&prefix);CHKERRQ(ierr);
    ierr = SNESCreate(((PetscObject)snes)->comm, &fas->downsmooth);CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(fas->downsmooth,prefix);CHKERRQ(ierr);
    ierr = SNESAppendOptionsPrefix(fas->downsmooth,"fas_levels_");CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)fas->downsmooth, (PetscObject)snes, 1);CHKERRQ(ierr); 
    ierr = SNESSetType(fas->downsmooth, pre_type);CHKERRQ(ierr);
    if (snes->ops->computefunction) {
      ierr = SNESSetFunction(fas->downsmooth,PETSC_NULL,snes->ops->computefunction,snes->funP);CHKERRQ(ierr);
    }
  }
  if (fas->downsmooth) {
    ierr = SNESSetTolerances(fas->downsmooth, 0.0, 0.0, 0.0, fas->max_down_it, 1000);CHKERRQ(ierr);
  }

  if (monflg) {
    if (fas->upsmooth)   ierr = SNESMonitorSet(fas->upsmooth,SNESMonitorDefault,PETSC_NULL,(PetscErrorCode (*)(void**))PetscViewerDestroy);CHKERRQ(ierr);
    if (fas->downsmooth) ierr = SNESMonitorSet(fas->downsmooth,SNESMonitorDefault,PETSC_NULL,(PetscErrorCode (*)(void**))PetscViewerDestroy);CHKERRQ(ierr);
  }

  /* recursive option setting for the smoothers */
  if (fas->next) {ierr = SNESSetFromOptions_FAS(fas->next);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESView_FAS"
PetscErrorCode SNESView_FAS(SNES snes, PetscViewer viewer)
{
  SNES_FAS   *fas = (SNES_FAS *) snes->data;
  PetscBool      iascii;
  PetscErrorCode ierr;
  PetscInt levels = fas->levels;
  PetscInt i;

  PetscFunctionBegin;
  ierr = PetscTypeCompare((PetscObject) viewer, PETSCVIEWERASCII, &iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer, "FAS, levels = %d\n",  fas->levels);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);
    for (i = levels - 1; i >= 0; i--) {
      ierr = PetscViewerASCIIPrintf(viewer, "level: %d\n",  fas->level);CHKERRQ(ierr);
      if (fas->upsmooth) {
        ierr = PetscViewerASCIIPrintf(viewer, "up-smoother on level %D\n",  fas->level);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPushTab(viewer);
        ierr = SNESView(fas->upsmooth, viewer);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPopTab(viewer);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer, "no up-smoother on level %D\n",  fas->level);CHKERRQ(ierr);
      }
      if (fas->downsmooth) {
        ierr = PetscViewerASCIIPrintf(viewer, "down-smoother on level %D\n",  fas->level);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPushTab(viewer);
        ierr = SNESView(fas->downsmooth, viewer);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPopTab(viewer);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer, "no down-smoother on level %D\n",  fas->level);CHKERRQ(ierr);
      }
      if (fas->next) fas = (SNES_FAS *)fas->next->data;
    }
    ierr = PetscViewerASCIIPopTab(viewer);
  } else {
    SETERRQ1(((PetscObject)snes)->comm,PETSC_ERR_SUP,"Viewer type %s not supported for SNESFAS",((PetscObject)viewer)->type_name);
  }
  PetscFunctionReturn(0);
}

/*

Defines the FAS cycle as:

fine problem: F(x) = 0
coarse problem: F^c(x) = b^c

b^c = F^c(I^c_fx^f - I^c_fF(x))

correction:

x = x + I(x^c - Rx)

 */

#undef __FUNCT__
#define __FUNCT__ "FASCycle_Private"
PetscErrorCode FASCycle_Private(SNES snes, Vec B, Vec X, Vec F) {

  PetscErrorCode ierr;
  Vec X_c, Xo_c, F_c, B_c;
  SNES_FAS * fas = (SNES_FAS *)snes->data;
  PetscInt i;

  PetscFunctionBegin;
  /* pre-smooth -- just update using the pre-smoother */
  if (fas->upsmooth) {
    ierr = SNESSolve(fas->upsmooth, B, X);CHKERRQ(ierr);
  } else if (snes->pc) {
    ierr = SNESSolve(snes->pc, B, X);CHKERRQ(ierr);
  }
  if (fas->next) {
    for (i = 0; i < fas->n_cycles; i++) {
      ierr = SNESComputeFunction(snes, X, F);CHKERRQ(ierr);
      X_c  = fas->next->vec_sol;
      Xo_c = fas->next->work[0];
      F_c  = fas->next->vec_func;
      B_c  = fas->next->work[1];
      /* inject the solution to coarse */
      ierr = MatRestrict(fas->restrct, X, Xo_c);CHKERRQ(ierr);
      ierr = VecPointwiseMult(Xo_c, fas->rscale, Xo_c);CHKERRQ(ierr);
      if (B) {
        ierr = VecAYPX(F, -1.0, B);CHKERRQ(ierr); /* F = B - F (defect) */
      } else {
        ierr = VecScale(F, -1.0);CHKERRQ(ierr);
      }

      /* restrict the defect */
      ierr = MatRestrict(fas->restrct, F, B_c);CHKERRQ(ierr);
      ierr = VecPointwiseMult(B_c,  fas->rscale, B_c);CHKERRQ(ierr);

      /* solve the coarse problem corresponding to F^c(x^c) = b^c = Rb + F^c(Rx) - RF(x) */
      fas->next->vec_rhs = PETSC_NULL;                                           /*unset the RHS to evaluate function instead of residual*/
      ierr = SNESComputeFunction(fas->next, Xo_c, F_c);CHKERRQ(ierr);
      ierr = VecAXPY(B_c, 1.0, F_c);CHKERRQ(ierr);                               /* add F_c(X) to the RHS */

      /* set initial guess of the coarse problem to the projected fine solution */
      ierr = VecCopy(Xo_c, X_c);CHKERRQ(ierr);

      /* recurse to the next level */
      ierr = FASCycle_Private(fas->next, B_c, X_c, F_c);CHKERRQ(ierr);

      /* correct as x <- x + I(x^c - Rx)*/
      ierr = VecAXPY(X_c, -1.0, Xo_c);CHKERRQ(ierr);
      ierr = MatInterpolateAdd(fas->interpolate, X_c, X, X);CHKERRQ(ierr);
    }
  }
    /* down-smooth -- just update using the down-smoother */
  if (fas->level != 0) {
    if (fas->downsmooth) {
      ierr = SNESSolve(fas->downsmooth, B, X);CHKERRQ(ierr);
    } else if (snes->pc) {
      ierr = SNESSolve(snes->pc, B, X);CHKERRQ(ierr);
    }
  }
  ierr = SNESComputeFunction(snes, X, F);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FASInitialGuess_Private"
PetscErrorCode FASInitialGuess_Private(SNES snes, Vec B, Vec X) {

  PetscErrorCode ierr;
  Vec X_c, B_c;
  SNES_FAS * fas;

  PetscFunctionBegin;
  fas = (SNES_FAS *)snes->data;
  /* pre-smooth -- just update using the pre-smoother */
  if (fas->level == 0) {
    if (fas->upsmooth) {
      ierr = SNESSolve(fas->upsmooth, B, X);CHKERRQ(ierr);
    } else if (snes->pc) {
      ierr = SNESSolve(snes->pc, B, X);CHKERRQ(ierr);
    }
  }
  if (fas->next) {
    X_c  = fas->next->vec_sol;
    B_c  = fas->next->work[0];
    /* inject the solution to coarse */
    ierr = MatRestrict(fas->restrct, X, X_c);CHKERRQ(ierr);
    ierr = VecPointwiseMult(X_c, fas->rscale, X_c);CHKERRQ(ierr);
    if (B) {
      ierr = MatRestrict(fas->restrct, B, B_c);CHKERRQ(ierr);
      ierr = VecPointwiseMult(B_c, fas->rscale, B_c);CHKERRQ(ierr);
    } else {
      B_c = PETSC_NULL;
    }
    /* recurse to the next level */
    ierr = FASInitialGuess_Private(fas->next, B_c, X_c);CHKERRQ(ierr);
    ierr = MatInterpolate(fas->interpolate, X_c, X);CHKERRQ(ierr);
  }
  /* down-smooth -- just update using the down-smoother */
  if (fas->level != 0) {
    if (fas->downsmooth) {
      ierr = SNESSolve(fas->downsmooth, B, X);CHKERRQ(ierr);
    } else if (snes->pc) {
      ierr = SNESSolve(snes->pc, B, X);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SNESSolve_FAS"

PetscErrorCode SNESSolve_FAS(SNES snes)
{
  PetscErrorCode ierr;
  PetscInt i, maxits;
  Vec X, F, B;
  PetscReal fnorm;
  PetscFunctionBegin;
  maxits = snes->max_its;            /* maximum number of iterations */
  snes->reason = SNES_CONVERGED_ITERATING;
  X = snes->vec_sol;
  F = snes->vec_func;
  B = snes->vec_rhs;

  /*norm setup */
  ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
  snes->iter = 0;
  snes->norm = 0.;
  ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
  if (snes->domainerror) {
    snes->reason = SNES_DIVERGED_FUNCTION_DOMAIN;
    PetscFunctionReturn(0);
  }
  ierr = VecNorm(F, NORM_2, &fnorm);CHKERRQ(ierr); /* fnorm <- ||F||  */
  if (PetscIsInfOrNanReal(fnorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"Infinite or not-a-number generated in norm");
  ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
  snes->norm = fnorm;
  ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
  SNESLogConvHistory(snes,fnorm,0);
  ierr = SNESMonitor(snes,0,fnorm);CHKERRQ(ierr);

  /* set parameter for default relative tolerance convergence test */
  snes->ttol = fnorm*snes->rtol;
  /* test convergence */
  ierr = (*snes->ops->converged)(snes,0,0.0,0.0,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
  if (snes->reason) PetscFunctionReturn(0);
  for (i = 0; i < maxits; i++) {
    /* Call general purpose update function */
    if (snes->ops->update) {
      ierr = (*snes->ops->update)(snes, snes->iter);CHKERRQ(ierr);
    }
    ierr = FASCycle_Private(snes, B, X, F);CHKERRQ(ierr);
    ierr = VecNorm(F, NORM_2, &fnorm);CHKERRQ(ierr); /* fnorm <- ||F||  */
    /* Monitor convergence */
    ierr = PetscObjectTakeAccess(snes);CHKERRQ(ierr);
    snes->iter = i+1;
    snes->norm = fnorm;
    ierr = PetscObjectGrantAccess(snes);CHKERRQ(ierr);
    SNESLogConvHistory(snes,snes->norm,0);
    ierr = SNESMonitor(snes,snes->iter,snes->norm);CHKERRQ(ierr);
    /* Test for convergence */
    ierr = (*snes->ops->converged)(snes,snes->iter,0.0,0.0,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
    if (snes->reason) break;
  }
  if (i == maxits) {
    ierr = PetscInfo1(snes, "Maximum number of iterations has been reached: %D\n", maxits);CHKERRQ(ierr);
    if (!snes->reason) snes->reason = SNES_DIVERGED_MAX_IT;
  }
  PetscFunctionReturn(0);
}
