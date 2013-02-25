/* Program usage:  exSG4 [-help] for all PETSc options
*/
static char help[] = "Simple program which does an ilu(0) factorization of a matrix in structgridgpu format. The original matrix is converted to a matrix of type aij. Both the original and the converted matrix are factored. The resulting factored matrices are compared for consistency. Performance is not an issue for factoring. Run time options: [-n] [-m] [-p] [-dim] [-REF] [-info 1 for more info]. All of these flags are required except for info. Note: It is preferable to use exSG2 while checking performance (especially for big inputs) as it has less memory footprint.\n\n";

#include <sys/time.h>
#include <petscconf.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include "../../impls/aij/seq/aij.h"
#include "../../impls/structgrid/matstructgrid.h"
#include <petscksp.h> // this includes all the below headers
//#include <petscsys.h >//      	- base PETSc routines   petscvec.h - vectors
//#include<petscmat.h>// 	- matrices
// #include<petscis.h>//     	- index sets            petscksp.h - Krylov subspace methods
//#include<petscviewer.h>// 	- viewers               petscpc.h  - preconditioners

//#define OMP// enable to check OPENMP version
#define NUM_THREADS 2
#define GPU //enable to check GPU versions

//#define PAPI// enable this to use PAPI directly without HPCToolkit and specify the required counters
#ifdef PAPI
#include"papi.h"
#define NUM_EVENTS 4
#endif

PetscInt m=1,n=1,p=1,dim=3,dof=1;
PetscInt nos;
PetscInt info=0;
PetscBool sgdump=PETSC_FALSE, comp_fact = PETSC_FALSE;
PetscReal normdiff = 1.0e-6;
long REP=1;
extern int OPENMP;
PetscRandom exRand;

double simple_rand() {
  PetscScalar result;
  PetscErrorCode ierr;
  ierr = PetscRandomGetValue( exRand, &result );CHKERRQ(ierr);
  return result;
}

double rtclock() {
  struct timezone tzp;
  struct timeval tp;
  gettimeofday (&tp, &tzp);
  return (1.0*tp.tv_sec + tp.tv_usec*1.0e-6);
}

#define GET_GRID_INDEX( m, n, im, in, ip ) ((im) + (in)*m + (ip)*m*n)

/*int correctness(double * y, double * yc, int size)
{
	int retval = 1;
	int i;
	for(i = 0; i < size; i++)
	{
		if((y[i]-yc[i])> normdiff)
		{
			printf("\n Values different at index: %d, y=%f, yc=%f",i,y[i],yc[i]);
			retval = 0;	
		}
	}
	return retval;
}
*/

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
#ifdef PAPI
int Events[NUM_EVENTS] = {PAPI_L1_DCM,PAPI_L2_TCM,PAPI_L3_TCM,PAPI_TLB_DM};
long_long values[NUM_EVENTS];
int e;
#endif
	dof=1;

	Mat 	mataij, matsggpu, convaij;
  	PetscErrorCode ierr;
  	PetscInt       i, nz=1, *dims, *starts, *rows, *cols;
	PetscScalar    *vals, one=1.0;
  	PetscMPIInt    size;
        char           sggpu_fname[PETSC_MAX_PATH_LEN];
        MatFactorInfo  facInfo;
        PetscBool      flg, gotTestFile = PETSC_FALSE, noWrapStencil = PETSC_FALSE;

        if (0) {  // runtime test for gpu?
          printf("No gpu available -- aborting\n");
          return 0;
        }

        printf("gpu available\n");

  	PetscInitialize(&argc,&args,(char *)0,help);
  	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

        ierr = PetscRandomCreate( PETSC_COMM_SELF, &exRand );
        ierr = PetscRandomSetType( exRand, PETSCRAND );
  	
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	To do: Can take command line arguments for m,n,p,dim and dof 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        printf("get options\n");
  	ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-p",&p,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-dim",&dim,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-dof",&dof,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-info",&info,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-sgdump",&sgdump,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-comp_fact",&comp_fact,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(PETSC_NULL,"-sgfile",sggpu_fname,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(PETSC_NULL,"-nowrap",&noWrapStencil,PETSC_NULL);CHKERRQ(ierr);
        if ( flg )
          gotTestFile = PETSC_TRUE;
	int rep=1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-REP",&rep,PETSC_NULL);CHKERRQ(ierr);
  	REP = (long)rep;

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Set dims[] and nz using n,dim and dof.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 	printf("dim=%d\n",dim);
	//dims = malloc(sizeof(PetscInt)*dim);
	dims = malloc(sizeof(PetscInt)*3);
	dims[2]=1;
	nz=1;
	if(dim>0){
	dims[0]=m;
	nz=nz*m;
	}
	if(dim>1){
	dims[1]=n;
	nz=nz*n;
	}
	if(dim>2){
	dims[2]=p;
	nz=nz*p;
	}
	nz=nz*dof;
	nos = (dim*2+1)*(2*dof-1);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create matrices.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	printf("Creating matrices of size %dx%d\n", nz, nz);
  	ierr = MatCreate(PETSC_COMM_WORLD,&mataij);CHKERRQ(ierr);
  	ierr = MatSetSizes(mataij,PETSC_DECIDE,PETSC_DECIDE,nz,nz);CHKERRQ(ierr);
  	MatSetType(mataij,MATSEQAIJ);
	MatSeqAIJSetPreallocation(mataij,nos,PETSC_NULL);

	ierr = MatCreate(PETSC_COMM_WORLD,&matsggpu);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsggpu,nz,nz,nz,nz);CHKERRQ(ierr);
	MatSetType(matsggpu,MATSEQSGGPU);
  	starts = malloc(sizeof(PetscInt)*dim);
	ierr = MatSetStencil(matsggpu,dim,dims,starts,dof);CHKERRQ(ierr);
	ierr = MatSeqSGGPUSetPreallocation(matsggpu,0,dof);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set values into matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscInt j,k,l,st;
	PetscInt nost = (dim*2+1);
	PetscInt inz, nnz = 0;
	PetscInt rowval = 0;
	cols = malloc(sizeof(PetscInt)*dof);
	vals = malloc(sizeof(PetscScalar)*dof);

	if ( noWrapStencil ) {
	  PetscInt im, in, ip, blockRow;
	  PetscInt *blockCols = malloc(sizeof(PetscInt)*nost);
	  for ( ip = 0; ip < p; ip++ )
	    for ( in = 0; in < n; in++ )
	      for ( im = 0; im < m; im++ ) {
		blockRow = GET_GRID_INDEX(m,n,im,in,ip);
		blockCols[0] = blockRow;
		nnz = 1;
		if ( ip > 0 )
		  blockCols[nnz++] = GET_GRID_INDEX(m,n,im,in,ip-1);
		if ( ip < p-1 )
		  blockCols[nnz++] = GET_GRID_INDEX(m,n,im,in,ip+1);
		if ( in > 0 )
		  blockCols[nnz++] = GET_GRID_INDEX(m,n,im,in-1,ip);
		if ( in < n-1 )
		  blockCols[nnz++] = GET_GRID_INDEX(m,n,im,in+1,ip);
		if ( im > 0 )
		  blockCols[nnz++] = GET_GRID_INDEX(m,n,im-1,in,ip);
		if ( im < m-1 )
		  blockCols[nnz++] = GET_GRID_INDEX(m,n,im+1,in,ip);
		for ( inz = 0; inz < nnz; inz++ ) {
		  for(l=0;l<dof;l++) {
		    rowval = blockRow*dof + l;
		    for(k=0;k<dof;k++) {
		      vals[k] = simple_rand();
		      cols[k] = blockCols[inz]*dof + k; 
		      // printf("val[%d,%d] = %f\n", rowval, cols[k], vals[k]);
		    }
		    ierr = MatSetValues(mataij,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		    ierr = MatSetValues(matsggpu,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
		  }
		}
	      }
	}
	else {
	  PetscInt lda2 = m*n;
	  PetscInt lda3 = m;

	  PetscInt *xval = malloc(sizeof(PetscInt)*nost);
	  PetscInt *idx = malloc(sizeof(PetscInt)*nost);
	  PetscInt *idy = malloc(sizeof(PetscInt)*nost);
	  PetscInt *idz = malloc(sizeof(PetscInt)*nost);
	  
	  PetscInt cnt=0;

	  idx[cnt] = 0; idy[cnt]=0; idz[cnt++]= 0;
	  if(dim>0)
	    {
	      idx[cnt] = -1; idy[cnt] = 0; idz[cnt++] = 0;
	      idx[cnt] = 1; idy[cnt] = 0; idz[cnt++] = 0;
	    }
	  if(dim>1)
	    {
	      idx[cnt] = 0; idy[cnt] = -1; idz[cnt++] = 0;
	      idx[cnt] = 0; idy[cnt] = 1; idz[cnt++] = 0;
	    }
	  if(dim>2)
	    {
	      idx[cnt] = 0; idy[cnt] = 0; idz[cnt++] = -1;
	      idx[cnt] = 0; idy[cnt] = 0; idz[cnt++] = 1;
	    }
	  for(st=0;st<nost;st++)
	    {
	      xval[st] = idx[st] + idy[st]*lda3 + idz[st]*lda2;
	      //printf("st=%d, idx=%d, idy=%d, idz=%d, xval=%d\n",st,idx[st],idy[st],idz[st],xval[st]);
	    }
	  printf("nost=%d, nos=%d, nz=%d\n",nost,nos,nz);
	  for(i=0;i<m*n*p;i++)
	    {
	      for(st=0;st<nost;st++)
        	{
		  j=xval[st]+i;
		  if(j>=0 && j<m*n*p)
		    {
		      //printf("i=%d,j=%d\n",i,j);
		      for(l=0;l<dof;l++)
			{
			  
			  for(k=0;k<dof;k++)	
			    {
			      vals[k] = simple_rand();
			      cols[k] = j*dof+k; 
			    }
			  rowval = i*dof+l;
			  ierr = MatSetValues(mataij,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
			  ierr = MatSetValues(matsggpu,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
			}
		    }
		}
	    }
	}

        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         AssemblyBegin/End as values can still remain in Cache
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatAssemblyBegin(mataij,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(mataij,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

        if ( sgdump ) {
          PetscViewer sgmatViewer, rhsViewer;
          Mat matFromFile;
          Vec RHS;
          PetscRandom rctx;
	  char sgmatFname[64], sgvecFname[64];

          ierr = VecCreateSeq(PETSC_COMM_WORLD,nz,&RHS);CHKERRQ(ierr);
          ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
          ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
          ierr = PetscRandomSetInterval(rctx,0.0,1.0);CHKERRQ(ierr);
          ierr = VecSetRandom(RHS,rctx);CHKERRQ(ierr);

	  // ierr = VecView(RHS,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	  // printf("Random RHS, nz = %d\n",nz);

/*           printf("\nOriginal sggpu matrix:\n"); */
/*           ierr = MatView(matsggpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */

          // dump matrix to a file
	  sprintf( sgmatFname, "sgmat-%d-%d-%d-%d.bin", m, n, p, dof );
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sgmatFname,FILE_MODE_WRITE, &sgmatViewer);CHKERRQ(ierr);
          ierr = MatView(matsggpu,sgmatViewer);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&sgmatViewer);CHKERRQ(ierr);

          // dump rhs vector to a file
	  sprintf( sgvecFname, "sgvec-%d-%d-%d-%d.bin", m, n, p, dof );
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sgvecFname,FILE_MODE_WRITE, &rhsViewer);CHKERRQ(ierr);
          ierr = VecView(RHS,rhsViewer);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&rhsViewer);CHKERRQ(ierr);

          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sgmatFname,FILE_MODE_READ, &sgmatViewer);CHKERRQ(ierr);
          ierr = MatCreate(PETSC_COMM_WORLD,&matFromFile);CHKERRQ(ierr);
          ierr = MatSetType(matFromFile,MATSEQSGGPU);CHKERRQ(ierr);
          ierr = MatLoad(matFromFile,sgmatViewer);CHKERRQ(ierr);
          printf("\nsggpu matrix from file:\n");
/*           ierr = MatView(matFromFile,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
          ierr = PetscViewerDestroy(&sgmatViewer);CHKERRQ(ierr);
          ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
          ierr = VecDestroy(&RHS);CHKERRQ(ierr);
        }

        ierr = MatConvert(matsggpu,MATSEQAIJ,MAT_INITIAL_MATRIX,&convaij);CHKERRQ(ierr);
        /* printf("\nOriginal aij matrix:\n"); */
        /* ierr = MatView(mataij,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
        /* printf("\nOriginal sggpu matrix:\n"); */
        /* ierr = MatView(matsggpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
        //printf("\nConverted sggpu matrix:\n");
        //ierr = MatView(convaij,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	for(i=0;i<dim;i++)
          printf("dims[%d] = %d\n", i,dims[i]);

	if ( comp_fact ) {
	  // ilu factor aij
	  Mat factaij, aij_l, aij_u;
	  IS rowIS = 0, colIS = 0;
	  ierr = MatFactorInfoInitialize(&facInfo);CHKERRQ(ierr);
	  facInfo.dt = 0.0;
	  facInfo.fill = 1.0;
	  facInfo.levels = PETSC_NULL;
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
	  PetscInt aij_l_nrows, aij_u_nrows, sggpu_l_nrows, sggpu_u_nrows;
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
	    for ( i = 0; i < aij_l_nrows && matsEqual; i++ ) {
	      ierr = MatGetRow( aij_l, i, &nColsAij, &colsAij, &valsAij );CHKERRQ(ierr);
	      ierr = MatGetRow( sggpu_l, i, &nColsSg, &colsSg, &valsSg );CHKERRQ(ierr);
	      if ( nColsAij != nColsSg ) {
		printf("L matrices do not have the same number of columns in row %d\n", i);
		matsEqual = PETSC_FALSE;
		continue;
	      }
	      for ( j = 0; j < nColsAij && matsEqual; j++ ) {
		if ( colsAij[j] != colsSg[j] ) {
		  printf("L matrices row %d Col %d not equal\n", i, j);
		  matsEqual = PETSC_FALSE;
		}
		else if ( ( fabs(valsAij[j] - valsSg[j]) > 0.001 ) &&
			  ( ( fabs(valsAij[j] - valsSg[j]) /
			      ( fabs(valsAij[j]) + fabs(valsSg[j]) ) ) > 0.001 ) ) {
		  printf("L matrices row %d Val %d not equal\n", i, j);
		  printf(" Values are %f and %f\n", valsAij[j], valsSg[j]);
		  matsEqual = PETSC_FALSE;
		}
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
	    for ( i = 0; i < aij_u_nrows && matsEqual; i++ ) {
	      ierr = MatGetRow( aij_u, i, &nColsAij, &colsAij, &valsAij );
	      ierr = MatGetRow( sggpu_u, i, &nColsSg, &colsSg, &valsSg );
	      if ( nColsAij != nColsSg ) {
		printf("U matrices do not have the same number of columns in row %d\n", i);
		matsEqual = PETSC_FALSE;
		continue;
	      }
	      for ( j = 0; j < nColsAij && matsEqual; j++ ) {
		if ( colsAij[j] != colsSg[j] ) {
		  printf("U matrices row %d Col %d not equal\n", i, j);
		  matsEqual = PETSC_FALSE;
		}
		else if ( ( fabs(valsAij[j] - valsSg[j]) > 0.001 ) &&
			  ( ( fabs(valsAij[j] - valsSg[j]) /
			      ( fabs(valsAij[j]) + fabs(valsSg[j]) ) ) > 0.001 ) ) {
		  printf("U matrices row %d Val %d not equal\n", i, j);
		  printf(" Values are %f and %f\n", valsAij[j], valsSg[j]);
		  matsEqual = PETSC_FALSE;
		}
	      }
	      ierr = MatRestoreRow( aij_u, i, &nColsAij, &colsAij, &valsAij );
	      ierr = MatRestoreRow( sggpu_u, i, &nColsSg, &colsSg, &valsSg );
	    }
	  }
	  
	  if ( matsEqual )
	    printf("ILU(0) factorizations are the same\n");

	  ierr = MatDestroy(&aij_l);CHKERRQ(ierr);
	  ierr = MatDestroy(&aij_u);CHKERRQ(ierr);
	  ierr = MatDestroy(&sggpu_l);CHKERRQ(ierr);
	  ierr = MatDestroy(&sggpu_u);CHKERRQ(ierr);
	  ierr = MatDestroy(&sggpu_product);CHKERRQ(ierr);
	}

	//CSR
        double start,end;
#ifdef PAPI
if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
#endif
	start = rtclock();

#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
for(e=0;e<NUM_EVENTS;e++)
	printf("Events[%d]= %lld\n",e,values[e]);

	end = rtclock();

if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
#endif

	ierr = MatDestroy(&matsggpu);CHKERRQ(ierr);
        ierr = MatDestroy(&mataij);CHKERRQ(ierr);
        ierr = MatDestroy(&convaij);CHKERRQ(ierr);

 	free(dims);
  	free(starts);
	free(cols);
	free(vals);
  	ierr = PetscFinalize();
  return 0;
}
