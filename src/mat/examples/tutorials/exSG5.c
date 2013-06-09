
static char help[] = "Reads a PETSc matrix and vector from a file and solves a linear system.\n\
It is copied and intended to move dirty codes from ksp/examples/tutorials/ex10.c and simplify ex10.c.\n\
  Input parameters include\n\
  -f0 <input_file> : first file to load (small system)\n\
  -f1 <input_file> : second file to load (larger system)\n\n\
  -trans  : solve transpose system instead\n\n";
/*
  This code  can be used to test PETSc interface to other packages.\n\
  Examples of command line options:       \n\
   ex30 -f0 <datafile> -ksp_type preonly  \n\
        -help -ksp_view                  \n\
        -num_numfac <num_numfac> -num_rhs <num_rhs> \n\
        -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package spooles or superlu or superlu_dist or mumps \n\
        -ksp_type preonly -pc_type cholesky -pc_factor_mat_solver_package spooles or mumps \n\
   mpiexec -n <np> ex30 -f0 <datafile> -ksp_type cg -pc_type asm -pc_asm_type basic -sub_pc_type icc -mat_type sbaij

   ./ex30 -f0 $D/small -mat_sigma -3.999999999999999 -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu -mat_superlu_conditionnumber -ckerror -mat_superlu_diagpivotthresh 0
   ./ex30 -f0 $D/small -mat_sigma -3.999999999999999 -ksp_type fgmres -pc_type hypre -pc_hypre_type boomeramg -ksp_type fgmres -ckError
   ./ex30 -f0 $D/small -mat_sigma -3.999999999999999 -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package petsc -pc_factor_shift_type NONZERO -pc_factor_shift_amount 1.e-5 -ckerror
 \n\n";
*/
/*T
   Concepts: KSP solving a linear system
   Processors: n
T*/

#include <petscksp.h>
#include <sys/time.h>
#include "papi.h"

#define MAX_CONV_RUNTIME (1000000*2*60)
#define NUM_EVENTS 5
//#define USE_PAPI 1
//not portable, even a little


long getTimeInMicroSeconds(){
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_sec*1000000L + (now.tv_usec);
}

long startTime = 0;

PetscErrorCode timedConvergenceTest(KSP ksp ,PetscInt iterNum, PetscReal rnorm, KSPConvergedReason* reason, void* cctx){
  PetscErrorCode ierr;
  long now = getTimeInMicroSeconds();
  //static int startTime = 0;


  ierr = KSPDefaultConverged(ksp, iterNum, rnorm, reason, cctx); CHKERRQ(ierr);

  //if(iterNum<=0)
  //  startTime = now;
  //printf("iter %d, start %ld, now %ld, time %ld, max %ld, reason %ld\n", iterNum, startTime, now, now-startTime, MAX_CONV_RUNTIME, *reason); 
  //if(*reason == KSP_CONVERGED_ITERATING && iterNum > 2) //limit num iters
  //  *reason = KSP_DIVERGED_NULL; //not sure what to use for the reason, NULL seems best?

  if(*reason == KSP_CONVERGED_ITERATING && now-startTime > MAX_CONV_RUNTIME) //do not know if we have converged
    *reason = KSP_DIVERGED_NULL; //not sure what to use for the reason, NULL seems best?

  return ierr;

}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  KSP            ksp, kspaij, kspdiag;
  Mat            A,B,C,Adiag;
  Vec            x,xaij,xdiag,bNatural,b,u,b2,b2aij,bdiag,b2diag,xlex;        /* approx solution, RHS, exact solution */
  PetscViewer    fd;              /* viewer */
  char           file[4][PETSC_MAX_PATH_LEN];     /* input file name */
  PetscBool      table = PETSC_FALSE,flg,flgB=PETSC_FALSE,trans=PETSC_FALSE,partition=PETSC_FALSE,initialguess = PETSC_FALSE;
  PetscBool      outputSoln=PETSC_FALSE;
  PetscErrorCode ierr;
  PetscInt       its,num_numfac,n,M;
  PetscReal      my_rtol;
  PetscReal      rnorm,enorm,myRtol;
  PetscLogDouble tsetup,tsetup1,tsetup2,tsolve,tsolve1,tsolve2;
  PetscBool      preload=PETSC_TRUE,diagonalscale,isSymmetric,ckrnorm=PETSC_TRUE,Test_MatDuplicate=PETSC_FALSE,ckerror=PETSC_FALSE;
  PetscBool      diagcomp=PETSC_FALSE, aijcomp=PETSC_FALSE, aijonly=PETSC_FALSE, useMyRtol=PETSC_FALSE;
  PetscMPIInt    rank;
  PetscScalar    sigma;
  PetscInt       m;
  PetscInt       bSize;

#if USE_PAPI
  int papi_err;
  uint32_t papi_events[NUM_EVENTS];

  /*{PAPI_TOT_CYC, PAPI_FP_OPS, sse_packed_code, FP_COMP_OPS_EXE_X87,
     FP_COMP_OPS_EXE_SSE_FP, FP_COMP_OPS_EXE_SSE_FP_SCALAR,
     FP_COMP_OPS_EXE_SSE_DOUBLE_PRECISION, FP_COMP_OPS_EXE_SSE_SINGLE_PRECISION,
     FP_COMP_OPS_EXE_SSE2_INTEGER*/

  long long papi_values[NUM_EVENTS];
  int papi_num_hwcntrs;
  int retval, sse_packed_code;

  // initialize papi and get some native event codes
  retval = PAPI_library_init(PAPI_VER_CURRENT);
  if (retval != PAPI_VER_CURRENT && retval > 0) {
    fprintf(stderr,"PAPI library version mismatch!\en");
    exit(1); }
  if ( PAPI_event_name_to_code( "FP_COMP_OPS_EXE:SSE_FP_PACKED", &sse_packed_code ) != PAPI_OK )
    printf("Could not get code for sse fp packed\n");

  papi_events[0] = PAPI_LD_INS;
  papi_events[1] = PAPI_SR_INS;
  papi_events[2] = PAPI_FP_OPS;
#endif

  //printf("Initialize PETSc\n");
  PetscInitialize(&argc,&args,(char *)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-table",&table,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-trans",&trans,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-partition",&partition,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-initialguess",&initialguess,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-output_solution",&outputSoln,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-ckrnorm",&ckrnorm,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-ckerror",&ckerror,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-diagcomp",&diagcomp,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-aijcomp",&aijcomp,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-aijonly",&aijonly,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-myrtol",&myRtol,&useMyRtol);
  void* conv_ctx;
  KSPDefaultConvergedCreate(&conv_ctx);//malloc(sizeof(KSPDefaultConvergedCtx));

#if USE_PAPI
  if((papi_err = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT ){
      printf("Library initialization error! \n");
      exit(1);
  }
  if ((papi_num_hwcntrs = PAPI_num_counters()) < PAPI_OK){
    printf("There are no counters available. \n");
    exit(1);
  }
  printf("There are %d counters in this system\n",papi_num_hwcntrs);
#endif

  /* 
     Determine files from which we read the two linear systems
     (matrix and right-hand-side vector).
  */
  ierr = PetscOptionsGetString(PETSC_NULL,"-f",file[0],PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscStrcpy(file[1],file[0]);CHKERRQ(ierr);
    preload = PETSC_FALSE;
  } else {
    ierr = PetscOptionsGetString(PETSC_NULL,"-f0",file[0],PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary file with the -f0 or -f option");
    ierr = PetscOptionsGetString(PETSC_NULL,"-f1",file[1],PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg) {preload = PETSC_FALSE;} /* don't bother with second system */
  }

  /* -----------------------------------------------------------
                  Beginning of linear solver loop
     ----------------------------------------------------------- */
  /* 
     Loop through the linear solve 2 times.  
      - The intention here is to preload and solve a small system;
        then load another (larger) system and solve it as well.
        This process preloads the instructions with the smaller
        system so that more accurate performance monitoring (via
        -log_summary) can be done with the larger one (that actually
        is the system of interest). 
  */
  PetscPreLoadBegin(preload,"Load system");

    /* - - - - - - - - - - - New Stage - - - - - - - - - - - - -
                           Load system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /* 
       Open binary file.  Note that we use FILE_MODE_READ to indicate
       reading from this file.
    */
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[PetscPreLoadIt],FILE_MODE_READ,&fd);CHKERRQ(ierr);
    
    /*
       Load the matrix and vector; then destroy the viewer.
    */
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetType(A,MATSEQSGGPU); CHKERRQ(ierr);
    ierr = MatLoad(A,fd);CHKERRQ(ierr);
  
    if (!preload){
      flg = PETSC_FALSE;
      ierr = PetscOptionsGetString(PETSC_NULL,"-rhs",file[2],PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
      //printf("rhs is %s\n", file[2]);
      if (flg){ /* rhs is stored in a separate file */
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr); 
        ierr = PetscInfo(0,"Loading RHS from file\n");CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[2],FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = MatGetLocalSize(A,&m,PETSC_NULL);CHKERRQ(ierr);
        ierr = VecCreate(PETSC_COMM_WORLD,&bNatural);CHKERRQ(ierr);
	ierr = VecLoad(bNatural, fd); CHKERRQ(ierr);
	ierr = VecGetSize(bNatural,&bSize);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
	ierr = VecSetSizes(b,PETSC_DECIDE,bSize);CHKERRQ(ierr);
	ierr = VecSetType(b,VECSEQGPU);CHKERRQ(ierr);
	ierr = VecCopy(bNatural,b);CHKERRQ(ierr);

        ierr = PetscObjectSetName((PetscObject)b, "Rhs vector");CHKERRQ(ierr);
      } else {
        /* if file contains no RHS, then use a vector of all ones */
        ierr = PetscInfo(0,"Using vector of ones for RHS\n");CHKERRQ(ierr);
        ierr = MatGetLocalSize(A,&m,PETSC_NULL);CHKERRQ(ierr);
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,m,PETSC_DECIDE);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecSet(b,1.0);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)b, "Rhs vector");CHKERRQ(ierr);
      }
    }
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr); 

    /* Test MatDuplicate() */
    if (Test_MatDuplicate){
      ierr = MatDuplicate(A,MAT_COPY_VALUES,&B);CHKERRQ(ierr);
      ierr = MatEqual(A,B,&flg);CHKERRQ(ierr);
      if (!flg){
        PetscPrintf(PETSC_COMM_WORLD,"  A != B \n");CHKERRQ(ierr);
      } 
      ierr = MatDestroy(&B);CHKERRQ(ierr); 
    }

    /* Add a shift to A */
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-mat_sigma",&sigma,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscOptionsGetString(PETSC_NULL,"-fB",file[2],PETSC_MAX_PATH_LEN,&flgB);CHKERRQ(ierr);
      if (flgB){
        /* load B to get A = A + sigma*B */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[2],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
        ierr = MatSetOptionsPrefix(B,"B_");CHKERRQ(ierr); /* e.g., ./ex30 -f0 <A> -fB <B> -mat_sigma 1.0 -B_mat_view_draw */
	ierr = MatLoad(B,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        ierr = MatAXPY(A,sigma,B,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr); /* A <- sigma*B + A */  
      } else {
        ierr = MatShift(A,sigma);CHKERRQ(ierr); 
      }
    }

    /* Make A singular for testing zero-pivot of ilu factorization        */
    /* Example: ./ex30 -f0 <datafile> -test_zeropivot -set_row_zero -pc_factor_shift_nonzero */
    flg  = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL, "-test_zeropivot", &flg,PETSC_NULL);CHKERRQ(ierr);
    if (flg) {
      PetscInt          row,ncols;
      const PetscInt    *cols;
      const PetscScalar *vals;
      PetscBool         flg1=PETSC_FALSE;
      PetscScalar       *zeros;
      row = 0;      
      ierr = MatGetRow(A,row,&ncols,&cols,&vals);CHKERRQ(ierr);     
      ierr = PetscMalloc(sizeof(PetscScalar)*(ncols+1),&zeros);
      ierr = PetscMemzero(zeros,(ncols+1)*sizeof(PetscScalar));CHKERRQ(ierr);
      flg1 = PETSC_FALSE;
      ierr = PetscOptionsGetBool(PETSC_NULL, "-set_row_zero", &flg1,PETSC_NULL);CHKERRQ(ierr);
      if (flg1){ /* set entire row as zero */
        ierr = MatSetValues(A,1,&row,ncols,cols,zeros,INSERT_VALUES);CHKERRQ(ierr);
      } else { /* only set (row,row) entry as zero */
        ierr = MatSetValues(A,1,&row,1,&row,zeros,INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

    /* 
       If the loaded matrix is larger than the vector (due to being padded 
       to match the block size of the system), then create a new padded vector.
    */
    
    ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
    if (m != n) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ, "This example is not intended for rectangular matrices (%d, %d)", m, n);
    ierr = MatGetSize(A,&M,PETSC_NULL);CHKERRQ(ierr);
    ierr = VecGetSize(b,&m);CHKERRQ(ierr);
    if (M != m) { /* Create a new vector b by padding the old one */
      PetscInt    j,mvec,start,end,indx;
      Vec         tmp;
      PetscScalar *bold;

      ierr = VecCreate(PETSC_COMM_WORLD,&tmp);CHKERRQ(ierr);
      ierr = VecSetSizes(tmp,n,PETSC_DECIDE);CHKERRQ(ierr);
      ierr = VecSetFromOptions(tmp);CHKERRQ(ierr);
      ierr = VecGetOwnershipRange(b,&start,&end);CHKERRQ(ierr);
      ierr = VecGetLocalSize(b,&mvec);CHKERRQ(ierr);
      ierr = VecGetArray(b,&bold);CHKERRQ(ierr);
      for (j=0; j<mvec; j++) {
        indx = start+j;
        ierr  = VecSetValues(tmp,1,&indx,bold+j,INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = VecRestoreArray(b,&bold);CHKERRQ(ierr);
      ierr = VecDestroy(&b);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(tmp);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(tmp);CHKERRQ(ierr);
      b = tmp;
    }
    ierr = VecDuplicate(b,&b2);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x, "Solution vector");CHKERRQ(ierr);
    ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)u, "True Solution vector");CHKERRQ(ierr);
    ierr = VecSet(x,0.0);CHKERRQ(ierr);

    if (ckerror){ /* Set true solution */
      ierr = VecSet(u,1.0);CHKERRQ(ierr);
      ierr = MatMult(A,u,b);CHKERRQ(ierr);
    }

    /* - - - - - - - - - - - New Stage - - - - - - - - - - - - -
                      Setup solve for system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    if (partition) {
      MatPartitioning mpart;
      IS              mis,nis,is;
      PetscInt        *count;
      PetscMPIInt     size;
      Mat             BB;
      ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
      ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
      ierr = PetscMalloc(size*sizeof(PetscInt),&count);CHKERRQ(ierr);
      ierr = MatPartitioningCreate(PETSC_COMM_WORLD, &mpart);CHKERRQ(ierr);
      ierr = MatPartitioningSetAdjacency(mpart, A);CHKERRQ(ierr);
      /* ierr = MatPartitioningSetVertexWeights(mpart, weight);CHKERRQ(ierr); */
      ierr = MatPartitioningSetFromOptions(mpart);CHKERRQ(ierr);
      ierr = MatPartitioningApply(mpart, &mis);CHKERRQ(ierr);
      ierr = MatPartitioningDestroy(&mpart);CHKERRQ(ierr);
      ierr = ISPartitioningToNumbering(mis,&nis);CHKERRQ(ierr);
      ierr = ISPartitioningCount(mis,size,count);CHKERRQ(ierr);
      ierr = ISDestroy(&mis);CHKERRQ(ierr);
      ierr = ISInvertPermutation(nis, count[rank], &is);CHKERRQ(ierr);
      ierr = PetscFree(count);CHKERRQ(ierr);
      ierr = ISDestroy(&nis);CHKERRQ(ierr);
      ierr = ISSort(is);CHKERRQ(ierr);
      ierr = MatGetSubMatrix(A,is,is,MAT_INITIAL_MATRIX,&BB);CHKERRQ(ierr);

      /* need to move the vector also */
      ierr = ISDestroy(&is);CHKERRQ(ierr);
      ierr = MatDestroy(&A);CHKERRQ(ierr);
      A    = BB;
    }

    /*
       We also explicitly time this stage via PetscGetTime()
    */
    ierr = PetscGetTime(&tsetup1);CHKERRQ(ierr);

    /*
       Create linear solver; set operators; set runtime options.
    */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    //KSPSetConvergenceTest(ksp, &timedConvergenceTest, conv_ctx, NULL);
    ierr = KSPSetInitialGuessNonzero(ksp,initialguess);CHKERRQ(ierr);
    num_numfac = 1;
    //ierr = PetscOptionsGetInt(PETSC_NULL,"-num_numfac",&num_numfac,PETSC_NULL);CHKERRQ(ierr);
    //while ( num_numfac-- ){
     
    ierr = KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);

    if ( diagcomp ) {
      //printf("Converting to diag numbered sggpu matrix for comparison\n");
      ierr = KSPCreate(PETSC_COMM_WORLD,&kspdiag);CHKERRQ(ierr);
      //KSPSetConvergenceTest(kspdiag, &timedConvergenceTest, conv_ctx, NULL);
      ierr = KSPSetInitialGuessNonzero(kspdiag,initialguess);CHKERRQ(ierr);

      ierr = MatConvertMatToDiag_SeqSGGPU(A,&Adiag);
      ierr = MatConvertVecLexDiag_SeqSGGPU(A,b,0,&bdiag);
      ierr = VecDuplicate(bdiag,&xdiag);CHKERRQ(ierr);
      ierr = VecDuplicate(bdiag,&b2diag);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)xdiag, "DIAG Solution vector");CHKERRQ(ierr);

      ierr = KSPSetOperators(kspdiag,Adiag,Adiag,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      ierr = KSPSetFromOptions(kspdiag);CHKERRQ(ierr);
    }

    if ( aijcomp || aijonly ) {
      ierr = KSPCreate(PETSC_COMM_WORLD,&kspaij);CHKERRQ(ierr);
      //KSPSetConvergenceTest(kspaij, &timedConvergenceTest, conv_ctx, NULL);
      ierr = KSPSetInitialGuessNonzero(kspaij,initialguess);CHKERRQ(ierr);
      ierr = MatConvert(A,MATSEQAIJ,MAT_INITIAL_MATRIX,&C);CHKERRQ(ierr);
      ierr = KSPSetOperators(kspaij,C,C,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      ierr = VecDuplicate(bNatural,&xaij);CHKERRQ(ierr);
      ierr = VecDuplicate(bNatural,&b2aij);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)xaij, "AIJ Solution vector");CHKERRQ(ierr);
      ierr = KSPSetFromOptions(kspaij);CHKERRQ(ierr);
      //printf("Converting to SeqAIJ matrix type for comparison\n");

      /* Check whether C is symmetric */
      flg  = PETSC_FALSE;
      ierr = PetscOptionsGetBool(PETSC_NULL, "-check_symmetry", &flg,PETSC_NULL);CHKERRQ(ierr);
      if (flg) {
	Mat Ctrans;
	ierr = MatTranspose(C, MAT_INITIAL_MATRIX,&Ctrans);
	ierr = MatEqual(C, Ctrans, &isSymmetric);
	if (isSymmetric) {
	  PetscPrintf(PETSC_COMM_WORLD,"input matrix is symmetric \n");CHKERRQ(ierr);
	} else {
	  PetscPrintf(PETSC_COMM_WORLD,"input matrix is non-symmetric \n");CHKERRQ(ierr);
	}
	ierr = MatDestroy(&Ctrans);CHKERRQ(ierr);
      }

    }

    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    if ( useMyRtol ) {
      printf("Setting ksp rtol to %f\n", myRtol);
      fflush(stdout);
      ierr = KSPSetTolerances(ksp,myRtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    }

      /* 
       Here we explicitly call KSPSetUp() and KSPSetUpOnBlocks() to
       enable more precise profiling of setting up the preconditioner.
       These calls are optional, since both will be called within
       KSPSolve() if they haven't been called already.
      */
      //ierr = KSPSetUp(ksp);CHKERRQ(ierr);
      //ierr = KSPSetUpOnBlocks(ksp);CHKERRQ(ierr);
      //ierr = PetscGetTime(&tsetup2);CHKERRQ(ierr);
      //tsetup = tsetup2 - tsetup1;

      /*
       Tests "diagonal-scaling of preconditioned residual norm" as used 
       by many ODE integrator codes including SUNDIALS. Note this is different
       than diagonally scaling the matrix before computing the preconditioner
      */
      //diagonalscale = PETSC_FALSE;
      //ierr = PetscOptionsGetBool(PETSC_NULL,"-diagonal_scale",&diagonalscale,PETSC_NULL);CHKERRQ(ierr);
      //if (diagonalscale) {
      //  PC       pc;
      //  PetscInt j,start,end,n;
      //  Vec      scale;
      //  
      //  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
      //  ierr = VecGetSize(x,&n);CHKERRQ(ierr);
      //  ierr = VecDuplicate(x,&scale);CHKERRQ(ierr);
      //  ierr = VecGetOwnershipRange(scale,&start,&end);CHKERRQ(ierr);
      //  for (j=start; j<end; j++) {
      //    ierr = VecSetValue(scale,j,((PetscReal)(j+1))/((PetscReal)n),INSERT_VALUES);CHKERRQ(ierr);
      //  }
      //  ierr = VecAssemblyBegin(scale);CHKERRQ(ierr);
      //  ierr = VecAssemblyEnd(scale);CHKERRQ(ierr);
      //  ierr = PCSetDiagonalScale(pc,scale);CHKERRQ(ierr);
      //  ierr = VecDestroy(&scale);CHKERRQ(ierr);
      //}

      /* - - - - - - - - - - - New Stage - - - - - - - - - - - - -
                           Solve system
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      /*
       Solve linear system; we also explicitly time this stage.
      */
      startTime = getTimeInMicroSeconds();

      ierr = PetscGetTime(&tsolve1);CHKERRQ(ierr);

#if USE_PAPI
      /*      int e;
      for(e=1;e<NUM_EVENTS;e++){
	printf("Count = %d (%ld)\n", e, papi_events[e-1]);
	if ( (papi_err = PAPI_start_counters(papi_events, e)) != PAPI_OK){
	  printf("Error starting counters! (%d)\n", papi_err); fflush(NULL);
	  return papi_err;
	}
	if ((papi_err=PAPI_stop_counters(papi_values, e)) != PAPI_OK){
	  printf("Error stoping counters! (%d)\n", papi_err); fflush(NULL);
          return papi_err;
	}
	}*/
      if ( (papi_err = PAPI_start_counters(papi_events, 3)) != PAPI_OK){
        printf("Error starting counters! (%d)\n", papi_err); fflush(NULL);
	return papi_err;
      }
#endif

      //if (trans) {
      //  ierr = KSPSolveTranspose(ksp,b,x);CHKERRQ(ierr);
      //  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
      //} else {
      //PetscInt  num_rhs=1;
      //ierr = PetscOptionsGetInt(PETSC_NULL,"-num_rhs",&num_rhs,PETSC_NULL);CHKERRQ(ierr);
        //while ( num_rhs-- ) {
      if ( !aijonly ) {
	//printf("Calling KSPSolve\n");
	//fflush(stdout);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	//printf("Done with KSPSolve\n");
	//fflush(stdout);
	//}
#if USE_PAPI
	if ((papi_err=PAPI_read_counters(papi_values, NUM_EVENTS)) != PAPI_OK){
	  printf("Error reading counters! (%d)\n", papi_err); fflush(NULL);
	  return papi_err;
	}
#endif
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp,&reason);
	PetscInt my_its = 0;
	if (reason<0) {
	  printf("Divergence: this should not happen.\n");
	} else {
	  KSPGetIterationNumber(ksp,&my_its);
	  //printf("\nConvergence in %d iterations.\n",(int)my_its);
	}

        ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
        if (ckrnorm){   /* Check residual for each rhs */
	  PetscScalar infinityNorm;
        //  if (trans) {
        //    ierr = MatMultTranspose(A,x,b2);CHKERRQ(ierr);
        //  } else {
            ierr = MatMult(A,x,b2);CHKERRQ(ierr);
        //  }
          ierr = VecAXPY(b2,-1.0,b);CHKERRQ(ierr);
	  //ierr = VecView(b2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          ierr = VecNorm(b2,NORM_2,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"  Number of iterations = %3D\n",its);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual norm %g\n",rnorm);CHKERRQ(ierr);
          //ierr = VecNorm(b2,NORM_INFINITY,&infinityNorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD," Infinity norm %g\n",infinityNorm);CHKERRQ(ierr);
        } 
        if (ckerror && !trans){  /* Check error for each rhs */
          /* ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
          ierr = VecAXPY(u,-1.0,x);CHKERRQ(ierr);
          ierr = VecNorm(u,NORM_2,&enorm);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"  Error norm %f\n",enorm);CHKERRQ(ierr);
          ierr = VecNorm(u,NORM_INFINITY,&enorm);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"  Error infinity norm %f\n",enorm);CHKERRQ(ierr);
        }
      
	//} /* while ( num_rhs-- ) */
	ierr = PetscGetTime(&tsolve2);CHKERRQ(ierr);
	tsolve = tsolve2 - tsolve1;

	/*
	  Write output (optinally using table for solver details).
	  - PetscPrintf() handles output for multiprocessor jobs 
          by printing from only one processor in the communicator.
	  - KSPView() prints information about the linear solver.
	*/
	if (table && ckrnorm) {
	  char        *matrixname,kspinfo[120];
	  PetscViewer viewer;

	  /*
	    Open a string viewer; then write info to it.
	  */
	  //ierr = PetscViewerStringOpen(PETSC_COMM_WORLD,kspinfo,120,&viewer);//CHKERRQ(ierr);
	  //ierr = KSPView(ksp,viewer);//CHKERRQ(ierr);
	  //kspinfo[0]=0;
	  ierr = PetscStrrchr(file[PetscPreLoadIt],'/',&matrixname);//CHKERRQ(ierr);
	  //ierr = PetscPrintf(PETSC_COMM_WORLD,"%-8.8s %3D %f %f %f %f %s \n",
	  //		   matrixname,its,rnorm,tsetup+tsolve,tsetup,tsolve,kspinfo);//CHKERRQ(ierr);
	  
	  //{PAPI_TOT_INS, PAPI_TOT_CYC, FP_COMP_OPS_EXE_SSE_FP, FP_COMP_OPS_EXE_X87, PAPI_FP_OPS,
	  //FP_COMP_OPS_EXE_SSE_FP_PACKED, FP_COMP_OPS_EXE_SSE_FP_SCALAR, FP_COMP_OPS_EXE_SSE_DOUBLE_PRECISION, FP_COMP_OPS_EXE_SSE_SINGLE_PRECISION, FP_COMP_OPS_EXE_SSE2_INTEGER};
#if USE_PAPI
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-8.8s %3D %f %f %ld %ld %ld %ld %ld %0.2lf %0.2lf %0.3lf\n",
			     matrixname,its,rnorm,tsolve, papi_values[0], papi_values[1], papi_values[2], papi_values[3], papi_values[4],
			     // papi_values[4], papi_values[5], papi_values[6], papi_values[7], papi_values[8], papi_values[9],
			     (papi_values[1]/tsolve)/1000000000, (papi_values[2]/tsolve)/1000000000, tsolve/its);//CHKERRQ(ierr);
#else
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s sggpu lex,%3D,%f\n",
			     file[2],its,tsolve);CHKERRQ(ierr);
#endif
	  
	  /*
	    Destroy the viewer
	  */
	  //ierr = PetscViewerDestroy(&viewer);//CHKERRQ(ierr);
	  //} 
	  
	  //ierr = PetscOptionsGetString(PETSC_NULL,"-solution",file[3],PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
	  //if (flg) {
	  //PetscViewer viewer;
	  //Vec         xstar;

	  //ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[3],FILE_MODE_READ,&viewer);CHKERRQ(ierr);
	  //ierr = VecCreate(PETSC_COMM_WORLD,&xstar);CHKERRQ(ierr);
	  //ierr = VecLoad(xstar,viewer);CHKERRQ(ierr);
	  //ierr = VecAXPY(xstar, -1.0, x);CHKERRQ(ierr);
	  //ierr = VecNorm(xstar, NORM_2, &enorm);CHKERRQ(ierr);
	  //ierr = PetscPrintf(PETSC_COMM_WORLD, "Error norm %A\n", enorm);CHKERRQ(ierr);
	  //ierr = VecDestroy(&xstar);CHKERRQ(ierr);
	  //ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	  //}
	  if (outputSoln) {
	    PetscViewer viewer;

	    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"solution.petsc",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
	    ierr = VecView(x, viewer);CHKERRQ(ierr);
	    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	  }
	
	  flg  = PETSC_FALSE;
	  ierr = PetscOptionsGetBool(PETSC_NULL, "-ksp_reason", &flg,PETSC_NULL);CHKERRQ(ierr);
	  if (flg){
	    KSPConvergedReason reason;
	    ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D\n", reason); 
	  }
	} /* while ( num_numfac-- ) */
      }


      if ( diagcomp ) {
	//  Solve linear system; we also explicitly time this stage.
	startTime = getTimeInMicroSeconds();
	ierr = PetscGetTime(&tsolve1);CHKERRQ(ierr);

	if ( useMyRtol )
	  ierr = KSPSetTolerances(kspdiag,myRtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

#if USE_PAPI
	if ( (papi_err = PAPI_read_counters(papi_events, 3)) != PAPI_OK){
	  printf("Error reading counters! (%d)\n", papi_err); fflush(NULL);
	  return papi_err;
	}
#endif
	
	ierr = KSPSolve(kspdiag,bdiag,xdiag);CHKERRQ(ierr);

#if USE_PAPI
	if ((papi_err=PAPI_stop_counters(papi_values, NUM_EVENTS)) != PAPI_OK){
	  printf("Error stoping counters! (%d)\n", papi_err); fflush(NULL);
	  return papi_err;
	}
#endif

	KSPConvergedReason reason;
	KSPGetConvergedReason(kspdiag,&reason);
	PetscInt my_its = 0;
	if (reason<0) {
	  printf("Divergence: this should not happen.\n");
	} else {
	  KSPGetIterationNumber(kspdiag,&my_its);
	  //printf("\nConvergence in %d iterations.\n",(int)my_its);
	}

        ierr = KSPGetIterationNumber(kspdiag,&its);CHKERRQ(ierr);
        if (ckrnorm){   /* Check residual for each rhs */
	  PetscScalar infinityNorm;
	  Vec xCopy, b2lex;
	  //  if (trans) {
	  //    ierr = MatMultTranspose(A,x,b2);CHKERRQ(ierr);
	  //  } else {
	  ierr = MatMult(Adiag,xdiag,b2diag);CHKERRQ(ierr);
	  //  }
          ierr = VecAXPY(b2diag,-1.0,bdiag);CHKERRQ(ierr);
          ierr = VecNorm(b2diag,NORM_2,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"  For Diagonal Numbering:\n");CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"  Number of iterations = %3D\n",its);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual norm %g\n",rnorm);CHKERRQ(ierr);
          //ierr = VecNorm(b2diag,NORM_INFINITY,&infinityNorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD," Infinity norm %g\n",infinityNorm);CHKERRQ(ierr);

	  //ierr = PetscPrintf(PETSC_COMM_WORLD,"Compare solution vector with lex numbering:\n");CHKERRQ(ierr);
	  ierr = VecDuplicate(x,&xCopy);CHKERRQ(ierr);
	  ierr = VecCopy(x,xCopy);CHKERRQ(ierr);

	  ierr = VecDuplicate(xdiag,&xlex);CHKERRQ(ierr);
	  ierr = MatConvertVecLexDiag_SeqSGGPU(A,xdiag,1,&xlex);

	  ierr = VecDuplicate(bdiag,&b2lex);CHKERRQ(ierr);
	  ierr = MatMult(A,xlex,b2lex);CHKERRQ(ierr);
          ierr = VecAXPY(b2lex,-1.0,b);CHKERRQ(ierr);
          ierr = VecNorm(b2lex,NORM_2,&rnorm);CHKERRQ(ierr);
	  //ierr = PetscPrintf(PETSC_COMM_WORLD,"xlex is xdiag converted to lex\n");CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"L2 norm of b - A*xlex %g\n",rnorm);CHKERRQ(ierr);
          //ierr = VecNorm(b2lex,NORM_INFINITY,&infinityNorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD," Infinity norm %g\n",infinityNorm);CHKERRQ(ierr);


	  //ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	  //ierr = VecView(xlex,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          //ierr = VecNorm(xCopy,NORM_2,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"L2 norm of x %f\n",rnorm);CHKERRQ(ierr);

          //ierr = VecNorm(xlex,NORM_2,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"L2 norm of xdiag %f\n",rnorm);CHKERRQ(ierr);

          //ierr = VecAXPY(xCopy,-1.0,xlex);CHKERRQ(ierr);
	  //ierr = VecView(xlex,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          //ierr = VecNorm(xCopy,NORM_2,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"L2 norm of difference %f\n",rnorm);CHKERRQ(ierr);
          //ierr = VecNorm(xCopy,NORM_INFINITY,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"Infinity norm of difference %f\n",rnorm);CHKERRQ(ierr);
        } 
        if (ckerror && !trans){  /* Check error for each rhs */
          /* ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
          ierr = VecAXPY(u,-1.0,xdiag);CHKERRQ(ierr);
          ierr = VecNorm(u,NORM_2,&enorm);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"  Error norm %f\n",enorm);CHKERRQ(ierr);
        }
      
	ierr = PetscGetTime(&tsolve2);CHKERRQ(ierr);
	tsolve = tsolve2 - tsolve1;

	if (table && ckrnorm) {
	  char        *matrixname,kspinfo[120];
	  PetscViewer viewer;

	  ierr = PetscStrrchr(file[PetscPreLoadIt],'/',&matrixname);//CHKERRQ(ierr);
#if USE_PAPI
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-8.8s %3D %f %f %ld %ld %ld %ld %ld %0.2lf %0.2lf %0.3lf\n",
			     matrixname,its,rnorm,tsolve, papi_values[0], papi_values[1], papi_values[2], papi_values[3], papi_values[4],
			     // papi_values[4], papi_values[5], papi_values[6], papi_values[7], papi_values[8], papi_values[9],
			     (papi_values[1]/tsolve)/1000000000, (papi_values[2]/tsolve)/1000000000, tsolve/its);//CHKERRQ(ierr);
#else
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s sggpu diag,%3D,%f\n",
			     file[2],its,tsolve);//CHKERRQ(ierr);
#endif
	  if (outputSoln) {
	    PetscViewer viewer;

	    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"solution.petsc",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
	    ierr = VecView(x, viewer);CHKERRQ(ierr);
	    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	  }
	
	  flg  = PETSC_FALSE;
	  ierr = PetscOptionsGetBool(PETSC_NULL, "-ksp_reason", &flg,PETSC_NULL);CHKERRQ(ierr);
	  if (flg){
	    KSPConvergedReason reason;
	    ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D\n", reason); 
	  }
       
	}

      }

      if ( aijcomp || aijonly ) {
	//  Solve linear system; we also explicitly time this stage.
	startTime = getTimeInMicroSeconds();
	ierr = PetscGetTime(&tsolve1);CHKERRQ(ierr);

	if ( useMyRtol )
	  ierr = KSPSetTolerances(kspaij,myRtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

#if USE_PAPI
	if ( (papi_err = PAPI_read_counters(papi_events, 3)) != PAPI_OK){
	  printf("Error reading counters! (%d)\n", papi_err); fflush(NULL);
	  return papi_err;
	}
#endif
	
	ierr = KSPSolve(kspaij,bNatural,xaij);CHKERRQ(ierr);

#if USE_PAPI
	if ((papi_err=PAPI_stop_counters(papi_values, NUM_EVENTS)) != PAPI_OK){
	  printf("Error stoping counters! (%d)\n", papi_err); fflush(NULL);
	  return papi_err;
	}
#endif

	KSPConvergedReason reason;
	KSPGetConvergedReason(kspaij,&reason);
	PetscInt my_its = 0;
	if (reason<0) {
	  printf("Divergence: this should not happen.\n");
	} else {
	  KSPGetIterationNumber(kspaij,&my_its);
	  //printf("\nConvergence in %d iterations.\n",(int)my_its);
	}

        ierr = KSPGetIterationNumber(kspaij,&its);CHKERRQ(ierr);
        if (ckrnorm){   /* Check residual for each rhs */
	  PetscScalar infinityNorm;
	  Vec xaijCopy;
	  //  if (trans) {
	  //    ierr = MatMultTranspose(A,x,b2);CHKERRQ(ierr);
	  //  } else {
	  ierr = MatMult(C,xaij,b2aij);CHKERRQ(ierr);
	  //  }
          ierr = VecAXPY(b2aij,-1.0,bNatural);CHKERRQ(ierr);
          ierr = VecNorm(b2aij,NORM_2,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"  Number of iterations = %3D\n",its);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual norm %g\n",rnorm);CHKERRQ(ierr);
          ierr = VecNorm(b2aij,NORM_INFINITY,&infinityNorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD," Infinity norm %g\n",infinityNorm);CHKERRQ(ierr);

	  ierr = VecDuplicate(xaij,&xaijCopy);CHKERRQ(ierr);
	  ierr = VecCopy(xaij,xaijCopy);CHKERRQ(ierr);
          ierr = VecAXPY(xaijCopy,-1.0,x);CHKERRQ(ierr);
          ierr = VecNorm(xaijCopy,NORM_2,&rnorm);CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_WORLD,"L2 norm of difference between x and xaij  %g\n",rnorm);CHKERRQ(ierr);
        } 
        if (ckerror && !trans){  /* Check error for each rhs */
          /* ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
          ierr = VecAXPY(u,-1.0,xaij);CHKERRQ(ierr);
          ierr = VecNorm(u,NORM_2,&enorm);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"  Error norm %f\n",enorm);CHKERRQ(ierr);
        }
      
	ierr = PetscGetTime(&tsolve2);CHKERRQ(ierr);
	tsolve = tsolve2 - tsolve1;

	if (table && ckrnorm) {
	  char        *matrixname,kspinfo[120];
	  PetscViewer viewer;

	  ierr = PetscStrrchr(file[PetscPreLoadIt],'/',&matrixname);//CHKERRQ(ierr);
#if USE_PAPI
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"%-8.8s %3D %f %f %ld %ld %ld %ld %ld %0.2lf %0.2lf %0.3lf\n",
			     matrixname,its,rnorm,tsolve, papi_values[0], papi_values[1], papi_values[2], papi_values[3], papi_values[4],
			     // papi_values[4], papi_values[5], papi_values[6], papi_values[7], papi_values[8], papi_values[9],
			     (papi_values[1]/tsolve)/1000000000, (papi_values[2]/tsolve)/1000000000, tsolve/its);//CHKERRQ(ierr);
#else
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s aij lex,%3D,%f\n",
			     file[2],its,tsolve);//CHKERRQ(ierr);
#endif
	  if (outputSoln) {
	    PetscViewer viewer;

	    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"solution.petsc",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
	    ierr = VecView(x, viewer);CHKERRQ(ierr);
	    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	  }
	
	  flg  = PETSC_FALSE;
	  ierr = PetscOptionsGetBool(PETSC_NULL, "-ksp_reason", &flg,PETSC_NULL);CHKERRQ(ierr);
	  if (flg){
	    KSPConvergedReason reason;
	    ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D\n", reason); 
	  }
       
	}

      }
    /* 
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
    */
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr); 
    ierr = MatDestroy(&A);CHKERRQ(ierr); ierr = VecDestroy(&b);CHKERRQ(ierr);
    ierr = VecDestroy(&u);CHKERRQ(ierr); ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b2);CHKERRQ(ierr);
    if ( aijcomp || aijonly ) { ierr = VecDestroy(&b2aij);CHKERRQ(ierr); }
    if ( aijcomp || aijonly ) { ierr = MatDestroy(&C);CHKERRQ(ierr); }
    if (flgB) { ierr = MatDestroy(&B);CHKERRQ(ierr); }
  PetscPreLoadEnd();
  /* -----------------------------------------------------------
                      End of linear solver loop
     ----------------------------------------------------------- */
#if USE_PAPI
  PAPI_shutdown();
#endif
  ierr = PetscFinalize();
  return 0;
}
