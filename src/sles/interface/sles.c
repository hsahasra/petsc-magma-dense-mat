#ifndef lint
static char vcid[] = "$Id: sles.c,v 1.16 1995/04/13 19:55:21 curfman Exp $";
#endif

#include "slesimpl.h"
#include "options.h"

/*@
   SLESPrintHelp - Prints SLES options.

   Input Parameter:
.  sles - the SLES context

   Notes:
   Keywords:  SLES, help
@*/
int SLESPrintHelp(SLES sles)
{
  VALIDHEADER(sles,SLES_COOKIE);
  fprintf(stderr,"SLES options:\n");
  KSPPrintHelp(sles->ksp);
  PCPrintHelp(sles->pc);
  return 0;
}

/*@
   SLESSetOptionsPrefix - Sets the prefix used for searching for all 
   SLES options in the database.

   Input Parameter:
.  sles - the SLES context
.  prefix - the prefix to prepend to all option names

   Notes:
   Keywords:  SLES, options, prefix, database
@*/
int SLESSetOptionsPrefix(SLES sles,char *prefix)
{
  VALIDHEADER(sles,SLES_COOKIE);
  KSPSetOptionsPrefix(sles->ksp,prefix);
  PCSetOptionsPrefix(sles->pc,prefix);
  return 0;
}

/*@
   SLESSetFromOptions - Sets various SLES parameters from user options.
   Also takes all KSP and PC options.

   Input Parameter:
.  sles - the SLES context

   Notes:
   Keywords:  SLES, options, database
@*/
int SLESSetFromOptions(SLES sles)
{
  VALIDHEADER(sles,SLES_COOKIE);
  KSPSetFromOptions(sles->ksp);
  PCSetFromOptions(sles->pc);
  return 0;
}
/*@
   SLESCreate - Creates a linear equation solver context.

   Output Parameter:
.  sles - the newly created SLES context

   Notes:
   Keywords:  SLES, create, context
@*/
int SLESCreate(SLES *outsles)
{
  int ierr;
  SLES sles;
  *outsles = 0;
  PETSCHEADERCREATE(sles,_SLES,SLES_COOKIE,0,MPI_COMM_WORLD);
  PLogObjectCreate(sles);
  if ((ierr = KSPCreate(&sles->ksp))) SETERR(ierr,0);
  if ((ierr = PCCreate(&sles->pc))) SETERR(ierr,0);
  PLogObjectParent(sles,sles->ksp);
  PLogObjectParent(sles,sles->pc);
  sles->setupcalled = 0;
  *outsles = sles;
  return 0;
}

/*@
   SLESDestroy - Destroys the SLES context.

   Input Parameters:
.  sles - the SLES context

   Notes:
   Keywords:  SLES, destroy, context
@*/
int SLESDestroy(SLES sles)
{
  int ierr;
  VALIDHEADER(sles,SLES_COOKIE);
  ierr = KSPDestroy(sles->ksp); CHKERR(ierr);
  ierr = PCDestroy(sles->pc); CHKERR(ierr);
  PLogObjectDestroy(sles);
  PETSCHEADERDESTROY(sles);
  return 0;
}
extern int PCPreSolve(PC,KSP),PCPostSolve(PC,KSP);
/*@
   SLESSolve - Solves a linear system.

   Input Parameters:
.  sles - the SLES context
.  b - the right hand side

   Output Parameters:
.  x - the approximate solution
.  its - the number of iterations used

   Notes:
   Keywords:  SLES, solve, linear system
@*/
int SLESSolve(SLES sles,Vec b,Vec x,int *its)
{
  int ierr;
  KSP ksp;
  PC  pc;
  VALIDHEADER(sles,SLES_COOKIE);
  ksp = sles->ksp; pc = sles->pc;
  KSPSetRhs(ksp,b);
  KSPSetSolution(ksp,x);
  KSPSetBinv(ksp,pc);
  if (!sles->setupcalled) {
    if ((ierr = PCSetVector(pc,b))) SETERR(ierr,0);
    if ((ierr = KSPSetUp(sles->ksp))) SETERR(ierr,0);
    if ((ierr = PCSetUp(sles->pc))) SETERR(ierr,0);
    sles->setupcalled = 1;
  }
  ierr = PCPreSolve(pc,ksp); CHKERR(ierr);
  ierr = KSPSolve(ksp,its); CHKERR(ierr);
  ierr = PCPostSolve(pc,ksp); CHKERR(ierr);
  return 0;
}

/*@
   SLESGetKSP - Returns the KSP context for a SLES solver.

   Input Parameter:
.  sles - the SLES context

   Output Parameter:
.  ksp - the Krylov space context

   Notes:  
   The user can then directly manipulate the KSP context to set various 
   options, etc.
   
   Keywords:  SLES, KSP, context
@*/
int SLESGetKSP(SLES sles,KSP *ksp)
{
  VALIDHEADER(sles,SLES_COOKIE);
  *ksp = sles->ksp;
  return 0;
}
/*@
   SLESGetPC - Returns the preconditioner (PC) context for a SLES solver.

   Input Parameter:
.  sles - the SLES context

   Output Parameter:
.  pc - the preconditioner context

   Notes:  
   The user can then directly manipulate the PC context to set various 
   options, etc.

   Keywords:  SLES, PC, context
@*/
int SLESGetPC(SLES sles,PC *pc)
{
  VALIDHEADER(sles,SLES_COOKIE);
  *pc = sles->pc;
  return 0;
}

#include "mat/matimpl.h"
/*@
   SLESSetOperators - Sets the matrix associated with the linear system
   and a (possibly) different one associated with the preconditioner. 

   Input Parameters:
.  sles - the sles context
.  mat - the matrix to use
.  pmat - matrix to be used in constructing the preconditioner, usually
          the same as mat.  If pmat is 0, the old preconditioner is reused.  
.  flag - use 0 or MAT_SAME_NONZERO_PATTERN

   Notes:
   The flag can be used to eliminate unnecessary repeated work in the 
   repeated solution of linear systems of the same size using the same 
   preconditioner.  The user can set flag = MAT_SAME_NONZERO_PATTERN to 
   indicate that the preconditioning matrix has the same nonzero pattern 
   during successive linear solves.

   Keywords:  SLES, operators, matrix, preconditioner
@*/
int SLESSetOperators(SLES sles,Mat mat,Mat pmat,int flag)
{
  VALIDHEADER(sles,SLES_COOKIE);
  VALIDHEADER(mat,MAT_COOKIE);
  if (pmat) {VALIDHEADER(pmat,MAT_COOKIE);}
  PCSetOperators(sles->pc,mat,pmat,flag);
  sles->setupcalled = 0;  /* so that next solve call will call setup */
  return 0;
}
