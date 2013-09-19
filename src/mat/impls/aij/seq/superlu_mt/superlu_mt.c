

/*
     Defines the data structure for the base matrix type (SeqAIJ)
*/
#include <../src/mat/impls/aij/seq/aij.h>    /*I "petscmat.h" I*/

EXTERN_C_BEGIN

#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#include <pcsp_defs.h>
#else
#include <pzsp_defs.h>
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
#include <pssp_defs.h>
#else
#include <pdsp_defs.h>
#endif
#endif
EXTERN_C_END

/*
     This is the data we are "ADDING" to the SeqAIJ matrix type to get the SuperLU matrix type
*/
typedef struct {
  SuperMatrix         A,L,U,B,X,AC;
  superlumt_options_t options;
  PetscInt            *perm_r2;
  PetscReal           *R, *C;
  equed_t             equed;
  PetscReal           rpg, rcond;
  MatStructure        flg;
  PetscScalar         *rhs_dup;
  superlu_memusage_t  mem_usage;
  Gstat_t             gstat;

  /* Flag to clean up (non-global) SuperLU objects during Destroy */
  PetscBool CleanUpSuperLU;
} Mat_SuperLU;

extern PetscErrorCode MatFactorInfo_SuperLU(Mat,PetscViewer);
extern PetscErrorCode MatLUFactorNumeric_SuperLU(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatDestroy_SuperLU(Mat);
extern PetscErrorCode MatView_SuperLU(Mat,PetscViewer);
extern PetscErrorCode MatAssemblyEnd_SuperLU(Mat,MatAssemblyType);
extern PetscErrorCode MatSolve_SuperLU(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SuperLU(Mat,Vec,Vec);
extern PetscErrorCode MatLUFactorSymbolic_SuperLU(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatDuplicate_SuperLU(Mat, MatDuplicateOption, Mat*);

/*
    Utility function
*/
#undef __FUNCT__
#define __FUNCT__ "MatFactorInfo_SuperLU"
PetscErrorCode MatFactorInfo_SuperLU(Mat A,PetscViewer viewer)
{
  PetscErrorCode      ierr;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"SuperLU-MT run parameters:\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    These are the methods provided to REPLACE the corresponding methods of the
   SeqAIJ matrix class. Other methods of SeqAIJ are not replaced
*/
#undef __FUNCT__
#define __FUNCT__ "MatLUFactorNumeric_SuperLU"
PetscErrorCode MatLUFactorNumeric_SuperLU(Mat F,Mat A,const MatFactorInfo *info)
{
  Mat_SuperLU    *lu = (Mat_SuperLU*)F->spptr;
  Mat_SeqAIJ     *aa;
  PetscInt       sinfo;

  PetscFunctionBegin;
  if (lu->flg == SAME_NONZERO_PATTERN) { /* successing numerical factorization */
    lu->options.fact = DOFACT;
    Destroy_SuperMatrix_Store(&lu->A);
    /* need some way of knowing there was previously a factorization that needs freeing */
    if (lu->mem_usage.total_needed > 0) {
      PetscStackCall("SuperLU:Destroy_SuperNode_Matrix",Destroy_SuperNode_Matrix(&lu->L));
      PetscStackCall("SuperLU:Destroy_CompCol_Matrix",Destroy_CompCol_Matrix(&lu->U));
    }
  }

  /* Create the SuperMatrix for lu->A=A^T:
       Since SuperLU likes column-oriented matrices,we pass it the transpose,
       and then solve A^T X = B in MatSolve(). */
  aa = (Mat_SeqAIJ*)(A)->data;
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
  PetscStackCall("SuperLU:cCreate_CompCol_Matrix",cCreate_CompCol_Matrix(&lu->A,A->cmap->n,A->rmap->n,aa->nz,(singlecomplex*)aa->a,aa->j,aa->i,SLU_NC,SLU_C,SLU_GE));
#else
  PetscStackCall("SuperLU:zCreate_CompCol_Matrix",zCreate_CompCol_Matrix(&lu->A,A->cmap->n,A->rmap->n,aa->nz,(doublecomplex*)aa->a,aa->j,aa->i,SLU_NC,SLU_Z,SLU_GE));
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
  PetscStackCall("SuperLU:sCreate_CompCol_Matrix",sCreate_CompCol_Matrix(&lu->A,A->cmap->n,A->rmap->n,aa->nz,aa->a,aa->j,aa->i,SLU_NC,SLU_S,SLU_GE));
#else
  PetscStackCall("SuperLU:dCreate_CompCol_Matrix",dCreate_CompCol_Matrix(&lu->A,A->cmap->n,A->rmap->n,aa->nz,aa->a,aa->j,aa->i,SLU_NC,SLU_D,SLU_GE));
#endif
#endif

  /* Numerical factorization */
  if (F->factortype == MAT_FACTOR_LU) {
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#else
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
#else
    pdgstrf_init(lu->options.nprocs, lu->options.fact, lu->options.trans, lu->options.refact, lu->options.panel_size, lu->options.relax,
		 lu->options.diag_pivot_thresh, lu->options.usepr, lu->options.drop_tol, lu->options.perm_c, lu->options.perm_r,
		 NULL,0, &lu->A, &lu->AC, &lu->options, &lu->gstat);
    PetscStackCall("SuperLU:pdgstrf",pdgstrf(&lu->options, &lu->AC, lu->perm_r2,&lu->L, &lu->U,&lu->gstat,&sinfo));
#endif
#endif
  } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Factor type not supported");

  if (sinfo > 0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_MAT_LU_ZRPVT,"Zero pivot in row %D",sinfo);
  else if (sinfo < 0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB, "info = %D, the %D-th argument in pcgstrf() had an illegal value", sinfo,-sinfo);

  lu->flg                = SAME_NONZERO_PATTERN;
  F->ops->solve          = MatSolve_SuperLU;
  F->ops->solvetranspose = MatSolveTranspose_SuperLU;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatDestroy_SuperLU"
PetscErrorCode MatDestroy_SuperLU(Mat A)
{
  PetscErrorCode ierr;
  Mat_SuperLU    *lu=(Mat_SuperLU*)A->spptr;

  PetscFunctionBegin;
  if (lu && lu->CleanUpSuperLU) { /* Free the SuperLU datastructures */
    PetscStackCall("SuperLU:Destroy_SuperMatrix_Store",Destroy_SuperMatrix_Store(&lu->A));
    PetscStackCall("SuperLU:Destroy_SuperMatrix_Store",Destroy_SuperMatrix_Store(&lu->B));
    PetscStackCall("SuperLU:Destroy_SuperMatrix_Store",Destroy_SuperMatrix_Store(&lu->X));
    if (lu->mem_usage.total_needed > 0) {
      PetscStackCall("SuperLU:Destroy_SuperNode_Matrix",Destroy_SuperNode_Matrix(&lu->L));
      PetscStackCall("SuperLU:Destroy_CompCol_Matrix",Destroy_CompCol_Matrix(&lu->U));
    }
  }
  if (lu) {
    ierr = PetscFree(lu->options.perm_r);CHKERRQ(ierr);
    ierr = PetscFree(lu->perm_r2);CHKERRQ(ierr);
    ierr = PetscFree(lu->options.perm_c);CHKERRQ(ierr);
    ierr = PetscFree(lu->R);CHKERRQ(ierr);
    ierr = PetscFree(lu->C);CHKERRQ(ierr);
    ierr = PetscFree(lu->rhs_dup);CHKERRQ(ierr);
      }
  ierr = PetscFree(A->spptr);CHKERRQ(ierr);

  /* clear composed functions */
  ierr = PetscObjectComposeFunction((PetscObject)A,"MatFactorGetSolverPackage_C",NULL);CHKERRQ(ierr);
  ierr = MatDestroy_SeqAIJ(A);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatView_SuperLU"
PetscErrorCode MatView_SuperLU(Mat A,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  PetscBool         iascii;
  PetscViewerFormat format;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_INFO) {
      ierr = MatFactorInfo_SuperLU(A,viewer);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatSolve_SuperLU_Private"
PetscErrorCode MatSolve_SuperLU_Private(Mat A,Vec b,Vec x)
{
  Mat_SuperLU      *lu = (Mat_SuperLU*)A->spptr;
  PetscScalar      *barray,*xarray;
  PetscErrorCode   ierr;
  PetscInt         info;
  PetscReal        ferr,berr;
  static PetscBool cite = PETSC_FALSE;
  int              nthreads = 1;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister("@article{superlu99,\n  author  = {James W. Demmel and Stanley C. Eisenstat and\n             John R. Gilbert and Xiaoye S. Li and Joseph W. H. Liu},\n  title = {A supernodal approach to sparse partial pivoting},\n  journal = {SIAM J. Matrix Analysis and Applications},\n  year = {1999},\n  volume  = {20},\n  number = {3},\n  pages = {720-755}\n}\n",&cite);CHKERRQ(ierr);

  lu->B.ncol = 1;   /* Set the number of right-hand side */
  ierr = VecGetArray(b,&barray);CHKERRQ(ierr);
  ierr = VecGetArray(x,&xarray);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
  ((DNformat*)lu->B.Store)->nzval = (singlecomplex*)barray;
  ((DNformat*)lu->X.Store)->nzval = (singlecomplex*)xarray;
#else
  ((DNformat*)lu->B.Store)->nzval = (doublecomplex*)barray;
  ((DNformat*)lu->X.Store)->nzval = (doublecomplex*)xarray;
#endif
#else
  ((DNformat*)lu->B.Store)->nzval = barray;
  ((DNformat*)lu->X.Store)->nzval = xarray;
#endif

  lu->options.fact = FACTORED; /* Indicate the factored form of A is supplied. */
  if (A->factortype == MAT_FACTOR_LU) {
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
    PetscStackCall("SuperLU:cgssvx",pcgssvx(nthreads,&lu->options, &lu->A, lu->options.perm_c, lu->options.perm_r,  &lu->equed, lu->R, lu->C,
                                     &lu->L, &lu->U, &lu->B, &lu->X, &lu->rpg, &lu->rcond, &ferr, &berr,
                                     &lu->mem_usage, &info));
#else
    PetscStackCall("SuperLU:zgssvx",pzgssvx(nthreads,&lu->options, &lu->A, lu->options.perm_c, lu->options.perm_r,  &lu->equed, lu->R, lu->C,
                                     &lu->L, &lu->U, &lu->B, &lu->X, &lu->rpg, &lu->rcond, &ferr, &berr,
                                     &lu->mem_usage,  &info));
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
    PetscStackCall("SuperLU:sgssvx",psgssvx(nthreads,&lu->options, &lu->A, lu->options.perm_c, lu->options.perm_r,  &lu->equed, lu->R, lu->C,
                                     &lu->L, &lu->U, lu->work, lu->lwork, &lu->B, &lu->X, &lu->rpg, &lu->rcond, &ferr, &berr,
                                     &lu->mem_usage, &info));
#else
    PetscStackCall("SuperLU:dgssvx",pdgssvx(nthreads,&lu->options, &lu->A, lu->options.perm_c, lu->options.perm_r,  &lu->equed, lu->R, lu->C,
                                     &lu->L, &lu->U, &lu->B, &lu->X, &lu->rpg, &lu->rcond, &ferr, &berr,
                                     &lu->mem_usage, &info));
#endif
#endif
  } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Factor type not supported");
  ierr = VecRestoreArray(b,&barray);CHKERRQ(ierr);
  ierr = VecRestoreArray(x,&xarray);CHKERRQ(ierr);

  if (info > 0) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"  Warning: gssvx() returns info %D\n",info);
  } else if (info < 0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB, "info = %D, the %D-th argument in gssvx() had an illegal value", info,-info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatSolve_SuperLU"
PetscErrorCode MatSolve_SuperLU(Mat A,Vec b,Vec x)
{
  Mat_SuperLU    *lu = (Mat_SuperLU*)A->spptr;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  lu->options.trans = TRANS;
  ierr = MatSolve_SuperLU_Private(A,b,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatSolveTranspose_SuperLU"
PetscErrorCode MatSolveTranspose_SuperLU(Mat A,Vec b,Vec x)
{
  Mat_SuperLU    *lu = (Mat_SuperLU*)A->spptr;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  lu->options.trans = NOTRANS;
  ierr = MatSolve_SuperLU_Private(A,b,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Note the r permutation is ignored
*/
#undef __FUNCT__
#define __FUNCT__ "MatLUFactorSymbolic_SuperLU"
PetscErrorCode MatLUFactorSymbolic_SuperLU(Mat F,Mat A,IS r,IS c,const MatFactorInfo *info)
{
  Mat_SuperLU *lu = (Mat_SuperLU*)(F->spptr);

  PetscFunctionBegin;
  lu->flg                 = DIFFERENT_NONZERO_PATTERN;
  lu->CleanUpSuperLU      = PETSC_TRUE;
  F->ops->lufactornumeric = MatLUFactorNumeric_SuperLU;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatFactorGetSolverPackage_seqaij_superlu"
PetscErrorCode MatFactorGetSolverPackage_seqaij_superlu(Mat A,const MatSolverPackage *type)
{
  PetscFunctionBegin;
  *type = MATSOLVERSUPERLU;
  PetscFunctionReturn(0);
}


/*MC
  MATSOLVERSUPERLU_MT = "superlu" - A solver package providing solvers LU for OpenMP sequential matrices
  via the external package SuperLU_MT.

  Use ./configure --download-superlu_mt to have PETSc installed with SuperLU-MT

   Notes: Do not confuse this with MATSOLVERSUPERLU_DIST which is for MPI parallel sparse solves

   Level: beginner

.seealso: PCLU, PCILU, MATSOLVERSUPERLU_DIST, MATSOLVERMUMPS, PCFactorSetMatSolverPackage(), MatSolverPackage
M*/

#undef __FUNCT__
#define __FUNCT__ "MatGetFactor_seqaij_superlu_mt"
PETSC_EXTERN PetscErrorCode MatGetFactor_seqaij_superlu_mt(Mat A,MatFactorType ftype,Mat *F)
{
  Mat            B;
  Mat_SuperLU    *lu;
  PetscErrorCode ierr;
  PetscInt       m=A->rmap->n,n=A->cmap->n,i;

  PetscFunctionBegin;
  ierr = MatCreate(PetscObjectComm((PetscObject)A),&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,A->rmap->n,A->cmap->n,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetType(B,((PetscObject)A)->type_name);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(B,0,NULL);CHKERRQ(ierr);

  if (ftype == MAT_FACTOR_LU || ftype == MAT_FACTOR_ILU) {
    B->ops->lufactorsymbolic  = MatLUFactorSymbolic_SuperLU;
    B->ops->ilufactorsymbolic = MatLUFactorSymbolic_SuperLU;
  } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Factor type not supported");

  B->ops->destroy = MatDestroy_SuperLU;
  B->ops->view    = MatView_SuperLU;
  B->factortype   = ftype;
  B->assembled    = PETSC_TRUE;           /* required by -ksp_view */
  B->preallocated = PETSC_TRUE;

  ierr = PetscNewLog(B,Mat_SuperLU,&lu);CHKERRQ(ierr);
  lu->options.nprocs            = 1;
  lu->options.fact              = DOFACT;
  lu->options.trans             = TRANS;
  lu->options.refact            = NO;
  lu->options.panel_size        = sp_ienv(1);
  lu->options.relax             = sp_ienv(2);
  lu->options.diag_pivot_thresh = 1.e-3;
  lu->options.usepr             = NO;
  lu->options.drop_tol          = 0.0;
  lu->options.SymmetricMode     = NO;
  lu->options.PrintStat         = NO;
  lu->options.perm_c            = NULL;
  lu->options.perm_r            = NULL;
  lu->options.work              = NULL;
  lu->options.lwork             = 0;

  StatAlloc(n, lu->options.nprocs, lu->options.panel_size, lu->options.relax, &lu->gstat);
  StatInit(n, lu->options.nprocs, &lu->gstat);

  ierr = PetscOptionsBegin(PetscObjectComm((PetscObject)A),((PetscObject)A)->prefix,"SuperLU Options","Mat");CHKERRQ(ierr);
  PetscOptionsEnd();

  /* Allocate spaces (notice sizes are for the transpose) */
  ierr = PetscMalloc(n*sizeof(PetscInt),&lu->options.perm_r);CHKERRQ(ierr);
  for (i=0; i<n; i++) lu->options.perm_r[i] = i;
  ierr = PetscMalloc(n*sizeof(PetscInt),&lu->perm_r2);CHKERRQ(ierr);
  ierr = PetscMalloc(m*sizeof(PetscInt),&lu->options.perm_c);CHKERRQ(ierr);
  for (i=0; i<n; i++) lu->options.perm_c[i] = i;
  ierr = PetscMalloc(n*sizeof(PetscScalar),&lu->R);CHKERRQ(ierr);
  ierr = PetscMalloc(m*sizeof(PetscScalar),&lu->C);CHKERRQ(ierr);

  /* create rhs and solution x without allocate space for .Store */
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
  PetscStackCall("SuperLU:cCreate_Dense_Matrix(",cCreate_Dense_Matrix(&lu->B, m, 1, NULL, m, SLU_DN, SLU_C, SLU_GE));
  PetscStackCall("SuperLU:cCreate_Dense_Matrix(",cCreate_Dense_Matrix(&lu->X, m, 1, NULL, m, SLU_DN, SLU_C, SLU_GE));
#else
  PetscStackCall("SuperLU:zCreate_Dense_Matrix",zCreate_Dense_Matrix(&lu->B, m, 1, NULL, m, SLU_DN, SLU_Z, SLU_GE));
  PetscStackCall("SuperLU:zCreate_Dense_Matrix",zCreate_Dense_Matrix(&lu->X, m, 1, NULL, m, SLU_DN, SLU_Z, SLU_GE));
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
  PetscStackCall("SuperLU:sCreate_Dense_Matrix",sCreate_Dense_Matrix(&lu->B, m, 1, NULL, m, SLU_DN, SLU_S, SLU_GE));
  PetscStackCall("SuperLU:sCreate_Dense_Matrix",sCreate_Dense_Matrix(&lu->X, m, 1, NULL, m, SLU_DN, SLU_S, SLU_GE));
#else
  PetscStackCall("SuperLU:dCreate_Dense_Matrix",dCreate_Dense_Matrix(&lu->B, m, 1, NULL, m, SLU_DN, SLU_D, SLU_GE));
  PetscStackCall("SuperLU:dCreate_Dense_Matrix",dCreate_Dense_Matrix(&lu->X, m, 1, NULL, m, SLU_DN, SLU_D, SLU_GE));
#endif
#endif

  ierr     = PetscObjectComposeFunction((PetscObject)B,"MatFactorGetSolverPackage_C",MatFactorGetSolverPackage_seqaij_superlu);CHKERRQ(ierr);
  B->spptr = lu;
  *F       = B;
  PetscFunctionReturn(0);
}

