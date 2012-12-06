#include <petsc-private/pcimpl.h>   /*I "petscpc.h" I*/

#define PCPOLYNOMIAL_MAXORDER 32

#define PCPOLYNOMIAL_NEUMANN      0
#define PCPOLYNOMIAL_LEASTSQUARES 1
#define PCPOLYNOMIAL_CHEBYSHEV    2
#define PCPOLYNOMIAL_NTYPES       3
static const char *PCPOLYNOMIAL_TYPE[64] = {
  "neumann","leastsquares","chebyshev"
};
static PetscErrorCode PCPolynomialEstimateExtremeEigenvalues(PC pc);

typedef struct {
  Vec       xplusy;
  Vec       tempvec;
  Vec       z[PCPOLYNOMIAL_MAXORDER+1];
  PetscBool known_max_eig,known_min_eig;/* if using given extreme eigenvalues */
  PetscReal max_eig,min_eig; /* estimates of extreme eigenvalues */
  PetscReal omega;
  PetscReal kcoeff[PCPOLYNOMIAL_MAXORDER+2];
  PetscInt  order;
  KSP       ksp;
  PetscInt  polytype;
} PC_Polynomial;



#undef __FUNCT__  
#define __FUNCT__ "PCSetUp_Polynomial"
static PetscErrorCode PCSetUp_Polynomial(PC pc)
{
  PC_Polynomial  *pcpoly = (PC_Polynomial*)pc->data;
  PetscInt       i;
  PetscBool      isgpu1,isgpu2;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!pc->setupcalled) {
    ierr = MatGetVecs(pc->pmat,&pcpoly->xplusy,&pcpoly->tempvec);CHKERRQ(ierr);
    /*
    ierr = PetscObjectTypeCompare((PetscObject)pc->pmat,MATSEQDIA,&issame);CHKERRQ(ierr);
    if (issame) {
      ierr = VecSetType(pcpoly->xplusy, VECSEQGPU);CHKERRQ(ierr);
      ierr = VecSetType(pcpoly->tempvec, VECSEQGPU);CHKERRQ(ierr);
    }
     */
    /*
    ierr = VecSetFromOptions(pcpoly->xplusy); CHKERRQ(ierr);
    ierr = VecSetFromOptions(pcpoly->tempvec); CHKERRQ(ierr);
     */
    ierr = PetscObjectTypeCompare((PetscObject)pc->pmat,MATAIJCUSP,&isgpu1);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)pc->pmat,MATSEQSGGPU,&isgpu2);CHKERRQ(ierr);
    if (isgpu1 || isgpu2) {
      ierr = VecSetType(pcpoly->xplusy,VECCUSP); CHKERRQ(ierr);
      ierr = VecSetType(pcpoly->tempvec,VECCUSP); CHKERRQ(ierr);
    }
    if (pcpoly->polytype == PCPOLYNOMIAL_CHEBYSHEV) {
      for (i=0;i<=pcpoly->order;i++) {
        ierr = VecDuplicate(pcpoly->tempvec,&pcpoly->z[i]);CHKERRQ(ierr);
      }
    }
    PetscLogObjectParent(pc,pcpoly->xplusy);
    PetscLogObjectParent(pc,pcpoly->tempvec);
    if (pcpoly->order < 0) {
      SETERRQ(PETSC_COMM_WORLD,1,"polynomial order must be greater than 0");
    }

    if (pcpoly->polytype == PCPOLYNOMIAL_LEASTSQUARES && pcpoly->order > 8) {
      SETERRQ(PETSC_COMM_WORLD,1,"polynomial leastsquares order must be in interval [0, 8]\n");
    }


    if (pcpoly->polytype == PCPOLYNOMIAL_CHEBYSHEV && 
        pcpoly->order > PCPOLYNOMIAL_MAXORDER) {
      SETERRQ1(PETSC_COMM_WORLD,1,"polynomial order must be in interval [0, %d]\n",PCPOLYNOMIAL_MAXORDER);
    }

  }
  ierr = PCPolynomialEstimateExtremeEigenvalues(pc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* -------------------------------------------------------------------------- */
/*
   PCApply_Polynomial - Applies the Polynomial preconditioner to a vector.

   Input Parameters:
.  pc - the preconditioner context
.  x - input vector

   Output Parameter:
.  y - output vector

   Application Interface Routine: PCApply()
 */
#undef __FUNCT__  
#define __FUNCT__ "PCApply_Polynomial"
static PetscErrorCode PCApply_Polynomial(PC pc,Vec x,Vec y)
{
  PC_Polynomial      *pcpoly = (PC_Polynomial*)pc->data;
  PetscInt i;
  PetscScalar scalar[9]; 
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Need to approximate omega=1/max(eigv(A)) */
  pcpoly->omega = 1.0/pcpoly->max_eig;

  /* Compute P(A)*x */
  if (pcpoly->polytype == PCPOLYNOMIAL_NEUMANN) {
    ierr = MatMult(pc->pmat,x,pcpoly->tempvec);CHKERRQ(ierr);
    ierr = VecScale(pcpoly->tempvec,pcpoly->omega);CHKERRQ(ierr);
    ierr = VecWAXPY(y,-1.0,pcpoly->tempvec,x);CHKERRQ(ierr);
    for (i=1;i<=pcpoly->order;i++) {
      ierr = VecWAXPY(pcpoly->xplusy,1.0,x,y);CHKERRQ(ierr);
      ierr = MatMult(pc->pmat,pcpoly->xplusy,pcpoly->tempvec);CHKERRQ(ierr);
      ierr = VecWAXPY(y,-pcpoly->omega,pcpoly->tempvec,pcpoly->xplusy);CHKERRQ(ierr);
    }
    ierr = VecAXPY(y,1.0,x);CHKERRQ(ierr);
    ierr = VecScale(y,pcpoly->omega);CHKERRQ(ierr);
  } else if (pcpoly->polytype == PCPOLYNOMIAL_LEASTSQUARES) {
    scalar[1] = 4/pcpoly->max_eig;
    for (i=2;i<pcpoly->order;i++) {
      scalar[i] = scalar[1]*scalar[i-1];
    }

    /* Reconstruct leastsquares coefficients */
    switch (pcpoly->order) {
    case 0:
      pcpoly->kcoeff[0] = 1.0;
      break;
    case 1:
      pcpoly->kcoeff[0] = 5.0;
      pcpoly->kcoeff[1] = -scalar[1];
      break;
    case 2:
      pcpoly->kcoeff[0] = 14.0;
      pcpoly->kcoeff[1] = -7*scalar[1];
      pcpoly->kcoeff[2] = scalar[2];
      break;
    case 3:
      pcpoly->kcoeff[0] = 30.0;
      pcpoly->kcoeff[1] = -27*scalar[1];
      pcpoly->kcoeff[2] = 9*scalar[2];
      pcpoly->kcoeff[3] = -scalar[3];
      break;
    case 4:
      pcpoly->kcoeff[0] = 55.0;
      pcpoly->kcoeff[1] = -77.0*scalar[1];
      pcpoly->kcoeff[2] = 44.0*scalar[2];
      pcpoly->kcoeff[3] = -11.0*scalar[3];
      pcpoly->kcoeff[4] = scalar[4];
      break;
    case 5:
      pcpoly->kcoeff[0] = 91.0;
      pcpoly->kcoeff[1] = -182.0*scalar[1];
      pcpoly->kcoeff[2] = 156.0*scalar[2];
      pcpoly->kcoeff[3] = -65.0*scalar[3];
      pcpoly->kcoeff[4] = 13.0*scalar[4];
      pcpoly->kcoeff[5] = -1.0*scalar[5];
      break;
    case 6:
      pcpoly->kcoeff[0] =  140.0;
      pcpoly->kcoeff[1] = -378.0*scalar[1];
      pcpoly->kcoeff[2] =  450.0*scalar[2];
      pcpoly->kcoeff[3] = -275.0*scalar[3];
      pcpoly->kcoeff[4] =   90.0*scalar[4];
      pcpoly->kcoeff[5] =  -15.0*scalar[5];
      pcpoly->kcoeff[6] =        scalar[6];
      break;

    case 7:
      pcpoly->kcoeff[0] =  204.0;
      pcpoly->kcoeff[1] = -714.0*scalar[1];
      pcpoly->kcoeff[2] = 1122.0*scalar[2];
      pcpoly->kcoeff[3] = -935.0*scalar[3];
      pcpoly->kcoeff[4] =  442.0*scalar[4];
      pcpoly->kcoeff[5] = -119.0*scalar[5];
      pcpoly->kcoeff[6] =   17.0*scalar[6];
      pcpoly->kcoeff[7] =       -scalar[7];
      break;

    case 8:
      pcpoly->kcoeff[0] =  285.0;
      pcpoly->kcoeff[1] =-1254.0*scalar[1];
      pcpoly->kcoeff[2] = 2508.0*scalar[2];
      pcpoly->kcoeff[3] =-2717.0*scalar[3];
      pcpoly->kcoeff[4] = 1729.0*scalar[4];
      pcpoly->kcoeff[5] = -665.0*scalar[5];
      pcpoly->kcoeff[6] =  152.0*scalar[6];
      pcpoly->kcoeff[7] =  -19.0*scalar[7];
      pcpoly->kcoeff[8] =        scalar[8];
      break;
    }
    ierr = VecCopy(x,y);CHKERRQ(ierr);
    ierr = VecScale(y,pcpoly->kcoeff[pcpoly->order]);CHKERRQ(ierr);
    for (i=pcpoly->order-1;i>=0;i--) {
      ierr = MatMult(pc->pmat,y,pcpoly->tempvec);CHKERRQ(ierr);
      ierr = VecWAXPY(y,pcpoly->kcoeff[i],x,pcpoly->tempvec);CHKERRQ(ierr);
    }
  } else if (pcpoly->polytype == PCPOLYNOMIAL_CHEBYSHEV) {
    PetscReal theta;
    PetscReal delta;
    PetscReal tmp;
    theta = (pcpoly->min_eig + pcpoly->max_eig)/2.0;
    delta = (pcpoly->max_eig - pcpoly->min_eig)/2.0;
    pcpoly->kcoeff[0] = 1.0;
    pcpoly->kcoeff[1] = theta/delta;
    pcpoly->kcoeff[2] = 2*pcpoly->kcoeff[1]*theta/delta
                        - pcpoly->kcoeff[0];
    ierr = VecCopy(x,pcpoly->z[0]);CHKERRQ(ierr);
    ierr = VecScale(pcpoly->z[0],1.0/theta);CHKERRQ(ierr);

    ierr = MatMult(pc->pmat,x,pcpoly->z[1]);CHKERRQ(ierr);
    tmp = theta*theta - delta*delta/2.0;
    if (tmp == 0.0) {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,"Bad eigenvalue estimates; max == min\n");
    }
    ierr = VecScale(pcpoly->z[1],-1.0/tmp);CHKERRQ(ierr);
    ierr = VecAXPY(pcpoly->z[1],2*theta/tmp,x);CHKERRQ(ierr);

    for (i=2;i<=pcpoly->order;i++) {
      pcpoly->kcoeff[i+1] = 2*pcpoly->kcoeff[i]*theta/delta
                            - pcpoly->kcoeff[i-1];
      tmp = 2.0*pcpoly->kcoeff[i]/pcpoly->kcoeff[i+1]/delta;
      ierr = MatMult(pc->pmat,pcpoly->z[i-1],pcpoly->z[i]);CHKERRQ(ierr);
      ierr = VecAXPY(pcpoly->z[i],-theta,pcpoly->z[i-1]);CHKERRQ(ierr);
      /* Note: this was wrong in Liang text, he omitted a delta */
      ierr = VecScale(pcpoly->z[i],-tmp);CHKERRQ(ierr);
      ierr = VecAXPY(pcpoly->z[i],tmp,x);CHKERRQ(ierr);
      ierr = VecAXPY(pcpoly->z[i],-pcpoly->kcoeff[i-1]/pcpoly->kcoeff[i+1],pcpoly->z[i-2]);CHKERRQ(ierr);
    }
    ierr = VecCopy(pcpoly->z[pcpoly->order],y);CHKERRQ(ierr);

  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCReset_Polynomial"
static PetscErrorCode PCReset_Polynomial(PC pc)
{
  PC_Polynomial      *pcpoly = (PC_Polynomial*)pc->data;
  PetscInt i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecDestroy(&pcpoly->xplusy);CHKERRQ(ierr);
  ierr = VecDestroy(&pcpoly->tempvec);CHKERRQ(ierr);
  if (pcpoly->polytype == PCPOLYNOMIAL_CHEBYSHEV) {
    for (i=0;i<pcpoly->order+1;i++) {
      ierr = VecDestroy(&pcpoly->z[i]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*
   PCDestroy_Polynomial - Destroys the private context for the Polynomial preconditioner
   that was created with PCCreate_Polynomial().

   Input Parameter:
.  pc - the preconditioner context

   Application Interface Routine: PCDestroy()
*/
#undef __FUNCT__  
#define __FUNCT__ "PCDestroy_Polynomial"
static PetscErrorCode PCDestroy_Polynomial(PC pc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCReset_Polynomial(pc);CHKERRQ(ierr);
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PCPolynomialEstimateExtremeEigenvalues"
PetscErrorCode PCPolynomialEstimateExtremeEigenvalues(PC pc)
{
  PC_Polynomial      *pcpoly = (PC_Polynomial*)pc->data;
  PetscBool      issame;
  PetscScalar le,ue;
  KSP ksp;
  PC subpc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if ((pcpoly->known_min_eig && pcpoly->known_max_eig) ||
      (pcpoly->known_max_eig && (pcpoly->polytype==PCPOLYNOMIAL_NEUMANN ||
                                 pcpoly->polytype==PCPOLYNOMIAL_LEASTSQUARES))){

    ierr = PetscInfo2(pc,"Extreme Eigenvalue estimates: (%G,%G)\n",pcpoly->min_eig,pcpoly->max_eig);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = PetscObjectTypeCompare((PetscObject)pc,PCPOLYNOMIAL,&issame);
 CHKERRQ(ierr);
  if (!issame) {
    PetscFunctionReturn(0);
  }
  ierr = KSPCreate(((PetscObject)pc)->comm,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,pc->pmat,pc->pmat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = VecSet(pcpoly->tempvec,1.0);CHKERRQ(ierr);
  ierr = KSPSetComputeSingularValues(ksp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,10);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&subpc);CHKERRQ(ierr);
  ierr = PCSetType(subpc,PCNONE);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,pcpoly->tempvec,pcpoly->xplusy);CHKERRQ(ierr);
  ierr = KSPComputeExtremeSingularValues(ksp,&ue,&le);CHKERRQ(ierr);
  if (!pcpoly->known_max_eig) {
    pcpoly->max_eig=ue;
  }
  if (!pcpoly->known_min_eig) {
    pcpoly->min_eig=le;
  }

  ierr = PetscInfo2(pc,"Extreme Eigenvalue estimates: (%G,%G)\n",pcpoly->min_eig,pcpoly->max_eig);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCSetFromOptions_Polynomial"
static PetscErrorCode PCSetFromOptions_Polynomial(PC pc)
{
  PC_Polynomial  *pcpoly = (PC_Polynomial*)pc->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsHead("Polynomial options");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-pc_polynomial_order","Order of polynomial preconditioner",
                         PETSC_NULL,pcpoly->order,&pcpoly->order,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEList("-pc_polynomial_type","neumann, chebyshev, or leastsquares",PETSC_NULL,PCPOLYNOMIAL_TYPE,
                           PCPOLYNOMIAL_NTYPES, PCPOLYNOMIAL_TYPE[pcpoly->polytype], &pcpoly->polytype,PETSC_NULL);CHKERRQ(ierr);

  ierr = PetscOptionsReal("-pc_polynomial_maxeig","use given value for maximum eigenvalue",PETSC_NULL,pcpoly->max_eig,&pcpoly->max_eig,&pcpoly->known_max_eig);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-pc_polynomial_mineig","use given value for minimum eigenvalue",PETSC_NULL,pcpoly->min_eig,&pcpoly->min_eig,&pcpoly->known_min_eig);CHKERRQ(ierr);

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------- */
/*
   PCCreate_Polynomial - Creates a Polynomial preconditioner context, PC_Polynomial, 
   and sets this as the private data within the generic preconditioning 
   context, PC, that was created within PCCreate().

   Input Parameter:
.  pc - the preconditioner context

   Application Interface Routine: PCCreate()
*/

/*MC
     PCPOLYNOMIAL - Polynomial (i.e. diagonal scaling preconditioning)

   Options Database Key:
+  -pc_polynomial_type     neumann, chebyshev, or leastsquares
-  -pc_polynomial_order    degree polynomial

   Level: beginner

  Concepts: Polynomial, diagonal scaling, preconditioners
M*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCCreate_Polynomial"
PetscErrorCode  PCCreate_Polynomial(PC pc)
{
  PC_Polynomial  *pcpoly;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /*
     Creates the private data structure for this preconditioner and
     attach it to the PC object.
  */
  ierr      = PetscNewLog(pc,PC_Polynomial,&pcpoly);CHKERRQ(ierr);
  pc->data  = (void*)pcpoly;
  pcpoly->order = 5;
  pcpoly->min_eig=0.0;
  pcpoly->max_eig=1.0;
  pcpoly->known_max_eig=pcpoly->known_min_eig=PETSC_FALSE;
  pcpoly->polytype = PCPOLYNOMIAL_CHEBYSHEV;
  

  /*
      Set the pointers for the functions that are provided above.
      Now when the user-level routines (such as PCApply(), PCDestroy(), etc.)
      are called, they will automatically call these functions.  Note we
      choose not to provide a couple of these functions since they are
      not needed.
  */
  pc->ops->apply               = PCApply_Polynomial;
  pc->ops->applytranspose      = PCApply_Polynomial;
  pc->ops->setup               = PCSetUp_Polynomial;
  pc->ops->reset               = PCReset_Polynomial;
  pc->ops->destroy             = PCDestroy_Polynomial;
  pc->ops->setfromoptions      = PCSetFromOptions_Polynomial;
  pc->ops->view                = 0;
  pc->ops->applyrichardson     = 0;
  pc->ops->applysymmetricleft  = PCApply_Polynomial;
  pc->ops->applysymmetricright = PCApply_Polynomial;
  PetscFunctionReturn(0);
}
EXTERN_C_END
