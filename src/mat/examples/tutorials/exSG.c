
/* Program usage:  mpiexec exSG [-help] [all PETSc options] [-n] dimension of vector(by default n=10)*/

static char help[] = "Simple program which does matrix vector multiplication using the default format aij and other formats namely, structgrid and structgridgpu. The resulting vectors are compared for consistency.\n\n";

#include <petscksp.h> // this includes all the below headers
//#include<petscsys.h >//      	- base PETSc routines   petscvec.h - vectors
//#include<petscmat.h>// 	- matrices
// #include<petscis.h>//     	- index sets            petscksp.h - Krylov subspace methods
//#include<petscviewer.h>// 	- viewers               petscpc.h  - preconditioners

int n,dim,dof,nos;
int normdiff = 1.0e-6;

unsigned int seed;
double simple_rand() {
  seed = (1103515245*seed+12345)%4294967296;
  return (1.0*seed)/4294967296;
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

	n=2;dim=3;dof=1;

	Vec            x, y,ysg, ysggpu;      
  	Mat            mat, matsg, matsggpu;           
  	PetscErrorCode ierr;
  	PetscInt       i,nz=1,*dims,*starts,*rows,*cols;
	PetscScalar    *vals,one=1.0;
  	PetscMPIInt    size;

  	PetscInitialize(&argc,&args,(char *)0,help);
  	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
  	
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	To do: Can take command line arguments for n,dim and dof 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
  
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Set nos, dims[] and nz using n,dim and dof.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	nos = dim*2 + 1;
 	dims = malloc(sizeof(PetscInt)*dim);
  	for(i=0;i<dim;i++)
	{
		dims[i]=n;
		nz= nz*n*dof;
	}

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  	ierr = VecSetSizes(x,PETSC_DECIDE,nz);CHKERRQ(ierr);
  	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
  	ierr = VecDuplicate(x,&ysg);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&ysggpu);CHKERRQ(ierr);
	
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create matrices.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatCreate(PETSC_COMM_WORLD,&mat);CHKERRQ(ierr);
  	ierr = MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,nz,nz);CHKERRQ(ierr);
  	ierr = MatCreate(PETSC_COMM_WORLD,&matsg);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsg,nz,nz,nz,nz);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&matsggpu);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsggpu,nz,nz,nz,nz);CHKERRQ(ierr);
	
  	//ierr = MatSetFromOptions(matsg);CHKERRQ(ierr);
  	MatSetType(mat,MATSEQAIJ);
  	MatSetType(matsg,MATSTRUCTGRID);
	MatSetType(matsggpu,MATSTRUCTGRIDGPU);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set stencils for Structgrid -matsg
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	starts = malloc(sizeof(PetscInt)*dim);
  	ierr = MatSetStencil(matsg,dim,dims,starts,dof);CHKERRQ(ierr);
	ierr = MatSetStencil(matsggpu,dim,dims,starts,dof);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set values into input vector and matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = VecSet(x,one);CHKERRQ(ierr);//this can be modified such that x holds random values
	//ierr = VecSetRandom(x,PETSC_NULL);

	rows = malloc(sizeof(PetscInt)*nz);
	cols = malloc(sizeof(PetscInt)*nos);
	vals = malloc(sizeof(PetscScalar)*nos);
	for(i=0;i<nz;i++)
		rows[i]=i;
	
	for(i=0;i<nos;i++)//This can be modified
	  vals[i]=1; //simple_rand();
	//This part should be changed. Right now it is hardcoded for n=2; dof=1; dim=3
	cols[0]=0;cols[1]=1;cols[2]=2;cols[3]=4;
   	ierr = MatSetValues(matsggpu,1,&rows[0],4,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(matsg,1,&rows[0],4,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[0],4,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	cols[0]=0;cols[1]=1;cols[2]=2;cols[3]=3;cols[4]=5;
	ierr = MatSetValues(matsggpu,1,&rows[1],5,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(matsg,1,&rows[1],5,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[1],5,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	cols[0]=0;cols[1]=1;cols[2]=2;cols[3]=3;cols[4]=4;cols[5]=6;
   	ierr = MatSetValues(matsggpu,1,&rows[2],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(matsg,1,&rows[2],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[2],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	cols[0]=1;cols[1]=2;cols[2]=3;cols[3]=4;cols[4]=5;cols[5]=7;
   	ierr = MatSetValues(matsggpu,1,&rows[3],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(matsg,1,&rows[3],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[3],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	cols[0]=0;cols[1]=2;cols[2]=3;cols[3]=4;cols[4]=5;cols[5]=6;
	ierr = MatSetValues(matsggpu,1,&rows[4],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(matsg,1,&rows[4],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[4],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	cols[0]=1;cols[1]=3;cols[2]=4;cols[3]=5;cols[4]=6;cols[5]=7;
   	ierr = MatSetValues(matsggpu,1,&rows[5],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues(matsg,1,&rows[5],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[5],6,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	cols[0]=2;cols[1]=4;cols[2]=5;cols[3]=6;cols[4]=7;
   	ierr = MatSetValues(matsggpu,1,&rows[6],5,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(matsg,1,&rows[6],5,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[6],5,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
	cols[0]=3;cols[1]=5;cols[2]=6;cols[3]=7;
   	ierr = MatSetValues(matsggpu,1,&rows[7],4,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(matsg,1,&rows[7],4,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   	ierr = MatSetValues(mat,1,&rows[7],4,cols,vals,INSERT_VALUES);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     AssemblyBegin/End as values can still remain in Cache
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatAssemblyBegin(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	printf("\nInputs:\n");
  	ierr = MatView(mat,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Compute solution vectors.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatMult(mat,x,y);CHKERRQ(ierr);
  	ierr = MatMult(matsg,x,ysg);CHKERRQ(ierr);
	ierr = MatMult(matsggpu,x,ysggpu);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	printf("\nOutput:\n");
  	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  	ierr = VecView(ysg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        ierr = VecView(ysggpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	PetscReal norm, normsg, normsggpu;
 	ierr = VecNorm(y,NORM_2,&norm); 
 	ierr = VecNorm(ysg,NORM_2,&normsg);
	ierr = VecNorm(ysggpu,NORM_2,&normsggpu);
	printf("Norm=%.3f\n",norm);
	printf("SG Norm=%.3f\n",normsg);
	printf("SG Norm=%.3f\n",normsggpu);
	if(norm-normsg > normdiff)
		printf("SG Test Failed\n");
	else if(norm-normsggpu > normdiff)
		printf("SGGPU Test Failed\n");
	else printf("Passed\n");
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Cleaning
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = VecDestroy(&x);CHKERRQ(ierr); 
	ierr = VecDestroy(&y);CHKERRQ(ierr);
  	ierr = VecDestroy(&ysg);CHKERRQ(ierr);
	ierr = VecDestroy(&ysggpu);CHKERRQ(ierr);
	ierr = MatDestroy(&mat);CHKERRQ(ierr);
	ierr = MatDestroy(&matsg);CHKERRQ(ierr);
	ierr = MatDestroy(&matsggpu);CHKERRQ(ierr);
	
  	ierr = PetscFinalize();
  return 0;
}
