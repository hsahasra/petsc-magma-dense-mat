#include <private/fortranimpl.h>
#include <petscsnes.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define sneslinesearchsetpostcheck_      SNESLINESEARCHSETPOSTCHECK
#define sneslinesearchset_               SNESLINESEARCHSET
#define sneslinesearchcubic_             SNESLINESEARCHCUBIC
#define sneslinesearchquadratic_         SNESLINESEARCHQUADRATIC
#define sneslinesearchno_                SNESLINESEARCHNO
#define sneslinesearchnonorms_           SNESLINESEARCHNONORMS
#define sneslinesearchsetprecheck_       SNESLINESEARCHSETPRECHECK
#define snessetupdate_                   SNESSETUPDATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define sneslinesearchsetpostcheck_      sneslinesearchsetpostcheck
#define sneslinesearchset_               sneslinesearchset
#define sneslinesearchcubic_             sneslinesearchcubic     
#define sneslinesearchquadratic_         sneslinesearchquadratic    
#define sneslinesearchno_                sneslinesearchno    
#define sneslinesearchnonorms_           sneslinesearchnonorms    
#define sneslinesearchsetprecheck_       sneslinesearchsetprecheck
#define snessetupdate_                   snessetupdate
#endif


/* These are not extern C because they are passed into non-extern C user level functions */
PetscErrorCode OurSNESLineSearch(SNES snes,void *ctx,Vec x,Vec f,Vec y,PetscReal fnorm,PetscReal xnorm,Vec g,Vec w,PetscReal*ynorm,PetscReal*gnorm,PetscBool  *flag)
{
  PetscErrorCode ierr = 0;
  (*(void (PETSC_STDCALL *)(SNES*,void*,Vec*,Vec*,Vec*,PetscReal*,PetscReal*,Vec*,Vec*,PetscReal*,PetscReal*,PetscBool *,PetscErrorCode*))(((PetscObject)snes)->fortran_func_pointers[6]))(&snes,ctx,&x,&f,&y,&fnorm,&xnorm,&g,&w,ynorm,gnorm,flag,&ierr);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode OurSNESLineSearchPostCheck(SNES snes,Vec x,Vec y,Vec z,void *checkCtx,PetscBool  *flag1,PetscBool  *flag2)
{
  PetscErrorCode ierr = 0;
  (*(void (PETSC_STDCALL *)(SNES*,Vec*,Vec*,Vec*,void*,PetscBool *,PetscBool *,PetscErrorCode*))(((PetscObject)snes)->fortran_func_pointers[7]))(&snes,&x,&y,&z,checkCtx,flag1,flag2,&ierr);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode OurSNESLineSearchPreCheck(SNES snes,Vec x,Vec y,void *checkCtx,PetscBool  *flag1)
{
  PetscErrorCode ierr = 0;
  (*(void (PETSC_STDCALL *)(SNES*,Vec*,Vec*,void*,PetscBool *,PetscErrorCode*))(((PetscObject)snes)->fortran_func_pointers[8]))(&snes,&x,&y,checkCtx,flag1,&ierr);CHKERRQ(ierr);
  return 0;
}
PetscErrorCode OurSNESSetUpdate(SNES snes,PetscInt b)
{
  PetscErrorCode ierr = 0;
  (*(void (PETSC_STDCALL *)(SNES*,PetscInt*,PetscErrorCode*))(((PetscObject)snes)->fortran_func_pointers[9]))(&snes,&b,&ierr);CHKERRQ(ierr);
  return 0;
}

EXTERN_C_BEGIN
void PETSC_STDCALL sneslinesearchsetpostcheck_(SNES *snes,void (PETSC_STDCALL *f)(SNES*,Vec*,Vec *,Vec *,void *,PetscBool *,PetscBool *,PetscErrorCode*),void *ctx,PetscErrorCode *ierr)
{
  PetscObjectAllocateFortranPointers(*snes,12);
  ((PetscObject)*snes)->fortran_func_pointers[7] = (PetscVoidFunction)f;
  *ierr = SNESLineSearchSetPostCheck(*snes,OurSNESLineSearchPostCheck,ctx);
}  

void PETSC_STDCALL sneslinesearchsetprecheck_(SNES *snes,void (PETSC_STDCALL *f)(SNES*,Vec*,Vec *,void *,PetscBool *,PetscErrorCode*),void *ctx,PetscErrorCode *ierr)
{
  PetscObjectAllocateFortranPointers(*snes,12);
  ((PetscObject)*snes)->fortran_func_pointers[8] = (PetscVoidFunction)f;
  *ierr = SNESLineSearchSetPreCheck(*snes,OurSNESLineSearchPreCheck,ctx);
}  

void PETSC_STDCALL snessetupdate_(SNES *snes,void (PETSC_STDCALL *f)(SNES*,PetscInt*,PetscErrorCode*),PetscErrorCode *ierr)
{
  PetscObjectAllocateFortranPointers(*snes,12);
  ((PetscObject)*snes)->fortran_func_pointers[9] = (PetscVoidFunction)f;
  *ierr = SNESSetUpdate(*snes,OurSNESSetUpdate);
}  
/* -----------------------------------------------------------------------------------------------------*/
void sneslinesearchcubic_(SNES *snes,void *lsctx,Vec *x,Vec *f,Vec *y,PetscReal*fnorm,PetscReal *xnorm,Vec *g,Vec *w,PetscReal *ynorm,PetscReal *gnorm,PetscBool  *flag,PetscErrorCode *ierr)
{
  *ierr = SNESLineSearchCubic(*snes,lsctx,*x,*f,*y,*fnorm,*xnorm,*g,*w,ynorm,gnorm,flag);
}
void sneslinesearchquadratic_(SNES *snes,void *lsctx,Vec *x,Vec *f,Vec *y,PetscReal*fnorm,PetscReal *xnorm,Vec *g,Vec *w,PetscReal *ynorm,PetscReal *gnorm,PetscBool  *flag,PetscErrorCode *ierr)
{
  *ierr = SNESLineSearchQuadratic(*snes,lsctx,*x,*f,*y,*fnorm,*xnorm,*g,*w,ynorm,gnorm,flag);
}
void sneslinesearchno_(SNES *snes,void *lsctx,Vec *x,Vec *f,Vec *y,PetscReal*fnorm,PetscReal *xnorm,Vec *g,Vec *w,PetscReal *ynorm,PetscReal *gnorm,PetscBool  *flag,PetscErrorCode *ierr)
{
  *ierr = SNESLineSearchNo(*snes,lsctx,*x,*f,*y,*fnorm,*xnorm,*g,*w,ynorm,gnorm,flag);
}
void sneslinesearchnonorms_(SNES *snes,void *lsctx,Vec *x,Vec *f,Vec *y,PetscReal*fnorm,PetscReal *xnorm,Vec *g,Vec *w,PetscReal *ynorm,PetscReal *gnorm,PetscBool  *flag,PetscErrorCode *ierr)
{
  *ierr = SNESLineSearchNoNorms(*snes,lsctx,*x,*f,*y,*fnorm,*xnorm,*g,*w,ynorm,gnorm,flag);
}

void PETSC_STDCALL sneslinesearchset_(SNES *snes,void (PETSC_STDCALL *f)(SNES*,void *,Vec*,Vec*,Vec*,PetscReal*,PetscReal*,Vec*,Vec*,PetscReal*,PetscReal*,PetscBool *,PetscErrorCode*),void *ctx,PetscErrorCode *ierr)
{
  PetscObjectAllocateFortranPointers(*snes,12);
  if ((PetscVoidFunction)f == (PetscVoidFunction)sneslinesearchcubic_) {
    *ierr = SNESLineSearchSet(*snes,SNESLineSearchCubic,ctx);
  } else if ((PetscVoidFunction)f == (PetscVoidFunction)sneslinesearchquadratic_) {
    *ierr = SNESLineSearchSet(*snes,SNESLineSearchQuadratic,ctx);
  } else if ((PetscVoidFunction)f == (PetscVoidFunction)sneslinesearchno_) {
    *ierr = SNESLineSearchSet(*snes,SNESLineSearchNo,ctx);
  } else if ((PetscVoidFunction)f == (PetscVoidFunction)sneslinesearchnonorms_) {
    *ierr = SNESLineSearchSet(*snes,SNESLineSearchNoNorms,ctx);
  } else {
    ((PetscObject)*snes)->fortran_func_pointers[6] = (PetscVoidFunction)f;
    *ierr = SNESLineSearchSet(*snes,OurSNESLineSearch,ctx);
  }
}

/* -----------------------------------------------------------------------------------------------------*/
EXTERN_C_END
