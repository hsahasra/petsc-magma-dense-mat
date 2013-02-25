static char help[] = "Reads a PETSc matrix from a file and verifies ilu(k) decomposition.\n\
  Input parameters include\n\
  -f <input_file> : matrix file to load\n\n";

#include <petscksp.h>
#include <sys/time.h>

#define FACTOR_TOL 0.01

long getTimeInMicroSeconds(){
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_sec*1000000L + (now.tv_usec);
}

long startTime = 0;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Mat 	         mataij, matsggpu, convaij;
  char           file[4][PETSC_MAX_PATH_LEN];     /* input file name */
  PetscViewer    fd;              /* viewer */
  MatFactorInfo  facInfo;
  PetscInt       level;
  PetscBool      flg, levelSet;
  PetscMPIInt    rank;

  printf("Initialize PETSc\n");
  PetscInitialize(&argc,&args,(char *)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-f",file[0],PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-level",&level,&levelSet);CHKERRQ(ierr);
  if ( !levelSet )
    level = PETSC_NULL;

  if ( flg ) {
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
    /* Load the matrix; then destroy the viewer. */
    ierr = MatCreate(PETSC_COMM_WORLD,&matsggpu);CHKERRQ(ierr);
    ierr = MatSetType(matsggpu,MATSEQSGGPU); CHKERRQ(ierr);
    ierr = MatLoad(matsggpu,fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  }

  ierr = MatConvert(matsggpu,MATSEQAIJ,MAT_INITIAL_MATRIX,&mataij);CHKERRQ(ierr);

  // ilu factor aij
  Mat factaij, aij_l, aij_u;
  IS rowIS = 0, colIS = 0;
  ierr = MatFactorInfoInitialize(&facInfo);CHKERRQ(ierr);
  facInfo.dt = 0.0;
  facInfo.fill = 1.0;
  facInfo.levels = level;
  ierr = MatGetOrdering(mataij,MATORDERINGNATURAL, &rowIS, &colIS);
  ierr = MatGetFactor(mataij,MATSOLVERPETSC,MAT_FACTOR_ILU,&factaij);
  ierr = MatILUFactorSymbolic(factaij,mataij,rowIS,colIS,&facInfo);CHKERRQ(ierr);
  ierr = MatLUFactorNumeric(factaij,mataij,&facInfo);CHKERRQ(ierr);
  ierr = MatConvertLU_SeqAIJ_SeqAIJ(factaij, &aij_l, &aij_u );CHKERRQ(ierr);

  // ilu factor sggpu
  Mat factsggpu, sggpu_l, sggpu_u, sggpu_product;
  ierr = MatGetOrdering(matsggpu,MATORDERINGNATURAL, &rowIS, &colIS);
  ierr = MatGetFactor(matsggpu,MATSOLVERPETSC,MAT_FACTOR_ILU,&factsggpu);
  ierr = MatILUFactorSymbolic(factsggpu,matsggpu,rowIS,colIS,&facInfo);CHKERRQ(ierr);
  ierr = MatLUFactorNumeric(factsggpu,matsggpu,&facInfo);CHKERRQ(ierr);
  ierr = MatConvertLU_SeqSGGPU_SeqAIJ(factsggpu, &sggpu_l, &sggpu_u );CHKERRQ(ierr);
  ierr = MatMatMult(sggpu_l,sggpu_u,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&sggpu_product);CHKERRQ(ierr);

  PetscBool matsEqual = PETSC_TRUE;;
  PetscInt i, j, jsg, aij_l_nrows, aij_u_nrows, sggpu_l_nrows, sggpu_u_nrows;
  PetscInt aij_l_ncols, aij_u_ncols, sggpu_l_ncols, sggpu_u_ncols;
  PetscInt nColsSg, nColsAij;
  const PetscInt *colsAij, *colsSg;
  const PetscScalar *valsAij, *valsSg;

  ierr = MatGetSize( aij_l, &aij_l_nrows, &aij_l_ncols );CHKERRQ(ierr);
  ierr = MatGetSize( sggpu_l, &sggpu_l_nrows, &sggpu_l_ncols );CHKERRQ(ierr);
  if ( aij_l_nrows != sggpu_l_nrows )
    printf("L matrices do not have the same number of rows\n");
  else if ( aij_l_ncols != sggpu_l_ncols ) {
    printf("L matrices do not have the same number of columns\n");
  }
  else {
    //    for ( i = 0; i < aij_l_nrows && matsEqual; i++ ) {
    for ( i = 0; i < aij_l_nrows; i++ ) {
      ierr = MatGetRow( aij_l, i, &nColsAij, &colsAij, &valsAij );CHKERRQ(ierr);
      ierr = MatGetRow( sggpu_l, i, &nColsSg, &colsSg, &valsSg );CHKERRQ(ierr);
      //if ( nColsAij != nColsSg ) {
      //printf("L matrices do not have the same number of columns in row %d\n", i);
      // matsEqual = PETSC_FALSE;
      // continue;
      // }
      //      for ( j = 0; j < nColsAij && matsEqual; j++ ) {
      jsg = 0;
      for ( j = 0; j < nColsAij; j++ ) {
	while ( (colsAij[j] > colsSg[jsg]) && (jsg < nColsSg) )
	  jsg++;
	if ( (jsg >= nColsSg) || (colsAij[j] < colsSg[jsg]) ) {
	  printf("L matrices row %d messed up at col %d\n", i, j);
	  matsEqual = PETSC_FALSE;
	  break;
	}
	else if ( ( fabs(valsAij[j] - valsSg[jsg]) > FACTOR_TOL ) &&
		  ( ( fabs(valsAij[j] - valsSg[jsg]) /
		      ( fabs(valsAij[j]) + fabs(valsSg[jsg]) ) ) > FACTOR_TOL ) ) {
	  printf("L matrices row %d Val %d not equal\n", i, j);
	  printf(" Values are %f and %f\n", valsAij[j], valsSg[jsg]);
	  matsEqual = PETSC_FALSE;
	}
	jsg++;
      }
      ierr = MatRestoreRow( aij_l, i, &nColsAij, &colsAij, &valsAij );CHKERRQ(ierr);
      ierr = MatRestoreRow( sggpu_l, i, &nColsSg, &colsSg, &valsSg );CHKERRQ(ierr);
    }
  }
  
  ierr = MatGetSize( aij_u, &aij_u_nrows, &aij_u_ncols );CHKERRQ(ierr);
  ierr = MatGetSize( sggpu_u, &sggpu_u_nrows, &sggpu_u_ncols );CHKERRQ(ierr);
  if ( aij_u_nrows != sggpu_u_nrows )
    printf("U matrices do not have the same number of rows\n");
  else if ( aij_u_ncols != sggpu_u_ncols ) {
    printf("U matrices do not have the same number of columns\n");
  }
  else {
    //for ( i = 0; i < aij_u_nrows && matsEqual; i++ ) {
    for ( i = 0; i < aij_u_nrows; i++ ) {
      ierr = MatGetRow( aij_u, i, &nColsAij, &colsAij, &valsAij );
      ierr = MatGetRow( sggpu_u, i, &nColsSg, &colsSg, &valsSg );
      //if ( nColsAij != nColsSg ) {
      //printf("U matrices do not have the same number of columns in row %d\n", i);
      //matsEqual = PETSC_FALSE;
      // continue;
      //}
      // for ( j = 0; j < nColsAij && matsEqual; j++ ) {
      jsg = 0;
      for ( j = 0; j < nColsAij; j++ ) {
	while ( (colsAij[j] > colsSg[jsg]) && (jsg < nColsSg) )
	  jsg++;
	if ( (jsg >= nColsSg) || (colsAij[j] < colsSg[jsg]) ) {
	  printf("U matrices row %d messed up at col %d\n", i, j);
	  matsEqual = PETSC_FALSE;
	  break;
	}
	else if ( ( fabs(valsAij[j] - valsSg[jsg]) > FACTOR_TOL ) &&
		  ( ( fabs(valsAij[j] - valsSg[jsg]) /
		      ( fabs(valsAij[j]) + fabs(valsSg[jsg]) ) ) > FACTOR_TOL ) ) {
	  printf("U matrices row %d Val %d not equal\n", i, j);
	  printf(" Values are %f and %f\n", valsAij[j], valsSg[jsg]);
	  matsEqual = PETSC_FALSE;
	}
	jsg++;
      }
      ierr = MatRestoreRow( aij_u, i, &nColsAij, &colsAij, &valsAij );
      ierr = MatRestoreRow( sggpu_u, i, &nColsSg, &colsSg, &valsSg );
    }
  }

  if ( matsEqual )
    printf("ILU(0) factorizations are the same\n");
  else
    printf("ILU(0) factorizations are different\n");

  ierr = PetscFinalize();
  return 0;
}

