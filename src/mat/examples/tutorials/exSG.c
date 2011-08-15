/* Program usage:  mpiexec ex1 [-help] for all PETSc options
*/
static char help[] = "Simple program which does matrix vector multiplication using the default format aij and other formats namely, structgrid and structgridgpu. The resulting vectors are compared for consistency. Options: [-n] [-m] [-p] [-dim] [-info 1 for more info]\n\n";

#include "../../impls/structgrid/matstructgrid.h"
#include <petscksp.h> // this includes all the below headers
//#include<petscsys.h >//      	- base PETSc routines   petscvec.h - vectors
//#include<petscmat.h>// 	- matrices
// #include<petscis.h>//     	- index sets            petscksp.h - Krylov subspace methods
//#include<petscviewer.h>// 	- viewers               petscpc.h  - preconditioners

PetscInt m=2,n=2,p=2,dim=3,dof=1;
PetscInt nos;
PetscInt info=0;
PetscReal normdiff = 1.0e-6;

double simple_rand() {
	int seed;
  	seed = (1103515245*seed+12345)%4294967296;
  	return (1.0*seed)/4294967296;
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

	dof=1;nos = dim*2 + 1;

	Vec            x, y, ysg, ygpu, ysggpu;      
  	Mat            mat, matgpu, matsg, matsggpu;           
  	PetscErrorCode ierr;
  	PetscInt       i, nz=1, *dims, *starts, *rows, *cols;
	PetscScalar    *vals, one=1.0;
  	PetscMPIInt    size;

  	PetscInitialize(&argc,&args,(char *)0,help);
  	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
  	
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	To do: Can take command line arguments for m,n,p,dim and dof 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-p",&p,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-dim",&dim,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-dof",&dof,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-info",&info,PETSC_NULL);CHKERRQ(ierr);
  
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Set dims[] and nz using n,dim and dof.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 	dims = malloc(sizeof(PetscInt)*dim);
	dims[0]=m;dims[1]=n;dims[2]=p;
	nz=m*n*p;
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
	ierr = VecDuplicate(x,&ygpu);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create matrices.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatCreate(PETSC_COMM_WORLD,&mat);CHKERRQ(ierr);
  	ierr = MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,nz,nz);CHKERRQ(ierr);
  	ierr = MatCreate(PETSC_COMM_WORLD,&matsg);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsg,nz,nz,nz,nz);CHKERRQ(ierr);
  	ierr = MatCreate(PETSC_COMM_WORLD,&matsggpu);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsggpu,nz,nz,nz,nz);CHKERRQ(ierr);
  	ierr = MatCreate(PETSC_COMM_WORLD,&matgpu);CHKERRQ(ierr);
  	ierr = MatSetSizes(matgpu,nz,nz,nz,nz);CHKERRQ(ierr);
  	//ierr = MatSetFromOptions(matsg);CHKERRQ(ierr);
  	MatSetType(mat,MATSEQAIJ);
  	MatSetType(matsg,MATSTRUCTGRID);
  	MatSetType(matsggpu,MATSTRUCTGRIDGPU);
  	MatSetType(matgpu,MATSEQAIJCUSP);
  
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set stencils for Structgrid -matsg
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	starts = malloc(sizeof(PetscInt)*dim);
  	ierr = MatSetStencil(matsg,dim,dims,starts,dof);CHKERRQ(ierr);
	ierr = MatSetStencil(matsggpu,dim,dims,starts,dof);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set values into input vector and matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	//ierr = VecSet(x,one);CHKERRQ(ierr);//this can be modified such that x holds random values
	ierr = VecSetRandom(x,PETSC_NULL);

	rows = malloc(sizeof(PetscInt)*nz);
	cols = malloc(sizeof(PetscInt)*nos);
	vals = malloc(sizeof(PetscScalar)*nos);
	for(i=0;i<nz;i++)
		rows[i]=i;
	
	Mat_SeqSG * sg = (Mat_SeqSG*) matsg->data;
	
	PetscInt k,l;
        PetscInt lda1 = m*n*p*dof;
        PetscInt lda2 = m*n*dof;
        PetscInt lda3 = m*dof;

	PetscInt *offset = malloc(sizeof(PetscInt)*nos);
	PetscInt *xval = malloc(sizeof(PetscInt)*nos);
	for(l=0;l<nos;l++)
        {
                offset[l] = l*lda1;
                xval[l] = sg->idx[l] + sg->idy[l]*lda3 + sg->idz[l]*lda2;
        }
	PetscInt count;
	for(i=0;i<nz;i++)
	{
		count=0;
		for(l=0;l<nos;l++)
        	{
                        vals[count] = simple_rand();
			if(xval[l]+i<nz)
				cols[count++] =  (xval[l]+i);    
		
		}
   		ierr = MatSetValues(mat,1,&rows[i],count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   		ierr = MatSetValues(matsg,1,&rows[i],count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   		ierr = MatSetValues(matsggpu,1,&rows[i],count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   		ierr = MatSetValues(matgpu,1,&rows[i],count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
			 
        }
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     AssemblyBegin/End as values can still remain in Cache
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matgpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matgpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if(info){
  	printf("\nInputs:\n");
  	ierr = MatView(mat,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  	ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Compute solution vectors.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatMult(mat,x,y);CHKERRQ(ierr);
  	ierr = MatMult(matsg,x,ysg);CHKERRQ(ierr);
	ierr = MatMult(matsggpu,x,ysggpu);CHKERRQ(ierr);
	ierr = MatMult(matgpu,x,ygpu);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	if(info){
	printf("\nOutput:\n");
	printf("Y - AIJ:\n");
  	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	printf("Y - Structgrid:\n");
  	ierr = VecView(ysg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	printf("Y - AIJ CUSP(GPU):\n");
        ierr = VecView(ygpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	printf("Y - Structgrid GPU:\n");
        ierr = VecView(ysggpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	}
	PetscReal norm, normsg,normsggpu,normgpu;
 	ierr = VecNorm(y,NORM_2,&norm); 
 	ierr = VecNorm(ysg,NORM_2,&normsg); 
 	ierr = VecNorm(ysggpu,NORM_2,&normsggpu); 
 	ierr = VecNorm(ygpu,NORM_2,&normgpu); 
	printf("AIJ Norm=%.3f\n",norm);
	printf("SG Norm=%.3f\n",normsg);
	printf("AIJ GPU Norm=%.3f\n",normgpu);
	printf("SGGPU Norm=%.3f\n",normsggpu);
	printf("Normdiff =%.3f\n",normdiff);
	if(abs(norm-normsg) > normdiff)
		printf("SG AVX/Openmp Test Failed\n");
	else 
		printf("SG AVX/Openmp Test Passed\n");
	if(abs(normgpu-normsggpu) > normdiff)
		printf("SG GPU Test Failed\n");
	else 
		printf("SG GPU Test Passed\n");
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
	ierr = MatDestroy(&matgpu);CHKERRQ(ierr);
	ierr = VecDestroy(&ygpu);CHKERRQ(ierr);

  	ierr = PetscFinalize();
  return 0;
}
