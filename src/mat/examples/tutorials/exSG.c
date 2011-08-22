/* Program usage:  mpiexec ex1 [-help] for all PETSc options
*/
static char help[] = "Simple program which does matrix vector multiplication using the default format aij and other formats namely, structgrid(avx, avx+openmp), aijcusp and structgridgpu. The resulting vectors are compared for consistency and also tested for performance. Options: [-n] [-m] [-p] [-dim] [-info 1 for more info]\n\n";

#include<sys/time.h>
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
PetscInt REP=1;
extern int OPENMP;

 double simple_rand() {
         int seed;
         seed = (1103515245*seed+12345)%4294967296;
         return (1000.0*seed)/4294967296;
 }

double rtclock() {
  struct timezone tzp;
  struct timeval tp;
  gettimeofday (&tp, &tzp);
  return (1.0*tp.tv_sec + tp.tv_usec*1.0e-6);
}

int correctness(double * y, double * yc, int size)
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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

	dof=1;nos = (dim*2+1)*(2*dof-1);//*dof?

	Vec  	x, y, ysg,ysgomp;
	Vec  	ygpu, ysggpu;      
  	Mat  	mat, matsg;
	Mat 	matgpu, matsggpu;           
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
	nz=m*n*p*dof;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  	ierr = VecSetSizes(x,PETSC_DECIDE,nz);CHKERRQ(ierr);
  	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
  	ierr = VecDuplicate(x,&ysg);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&ysgomp);CHKERRQ(ierr);
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
     Compute solution vectors and Performance test.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	double start,end;
	printf("Reps = %d, m=%d,n=%d,p=%d\n",REP,m,n,p);

	//CSR
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(mat,x,y);CHKERRQ(ierr);
	end = rtclock();
	printf("\nCSR :\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,REP*2*sg->stpoints*sg->nz/((end-start)*1024*1024*1024)); 
	fflush(stdout);	
	//SG (AVX)
	OPENMP = 0;
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysg);CHKERRQ(ierr);
	end = rtclock();
	printf("\nSG -AVX :\n");
	printf("Time =%.3f \nGFLOPS= %.3f\n",end-start,REP*2*sg->stpoints*sg->nz/((end-start)*1024*1024*1024)); 
 	
	//SG (AVX+OPENMP)
	OPENMP=1;
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysgomp);CHKERRQ(ierr);
	end = rtclock();
	printf("\nSG - AVX + OPENMP :\n");
	printf("Time =%.3f \nGFLOPS= %.3f\n",end-start,REP*2*sg->stpoints*sg->nz/((end-start)*1024*1024*1024)); 

	//CSR GPU
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matgpu,x,ygpu);CHKERRQ(ierr);
	end = rtclock();
	printf("\nCSR - GPU:\n");
	printf("Time =%.3f \nGFLOPS= %.3f",end-start,REP*2*sg->stpoints*sg->nz/((end-start)*1024*1024*1024)); 

	//SG (GPU)
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsggpu,x,ysggpu);CHKERRQ(ierr);
	end = rtclock();
	printf("\nSG - GPU:\n");
	printf("Time =%.3f \nGFLOPS= %.3f",end-start,REP*2*sg->stpoints*sg->nz/((end-start)*1024*1024*1024)); 


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	if(info){
	printf("\nOutput:\n");
	printf("Y - CSR:\n");
  	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	printf("Y - Structgrid AVX:\n");
  	ierr = VecView(ysg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	printf("Y - Structgrid AVX + OPENMP:\n");
  	ierr = VecView(ysgomp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	printf("Y - CSR CUSP(GPU):\n");
        ierr = VecView(ygpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	printf("Y - Structgrid GPU:\n");
        ierr = VecView(ysggpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	}

        //PetscScalar *xvec;
        //PetscInt xsize;
        //ierr = VecGetSize(ysggpu,&xsize);CHKERRQ(ierr);
        //ierr = VecGetArray(ysggpu,&xvec);CHKERRQ(ierr);
        //printf("xsize: %d\n",xsize);
        //if(xsize>10) xsize=10;
        //for(i=0;i<xsize;i++)printf("Xgpu[%d]: %f\n",i,xvec[i]);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Correctness test
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//	if(REP==1){call correctness instead?
		printf("\n\nCorrectness test : \n");
		PetscReal norm, normsg,normsgomp;
		PetscReal normsggpu;
		
		ierr = VecAXPY(ysg,-1,y);CHKERRQ(ierr);
//		ierr = VecView(ysg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	 	ierr = VecNorm(ysg,NORM_2,&normsg);CHKERRQ(ierr);
		printf("SG(AVX) Norm        = %.6f\n",normsg);
		if(normsg > normdiff)
			printf("SG AVX Test Failed\n");
		else 
			printf("SG AVX Test Passed\n");
  	
//		ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		ierr = VecAXPY(ysgomp,-1,y);CHKERRQ(ierr);
//		ierr = VecView(ysgomp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	 	ierr = VecNorm(ysgomp,NORM_2,&normsgomp); 
		printf("SG(AVX+OPENMP) Norm = %.6f\n",normsgomp);
		if(normsgomp > normdiff)
			printf("SG AVX+Openmp Test Failed\n");
		else 
			printf("SG AVX+Openmp Test Passed\n");

		ierr = VecAXPY(ysggpu,-1,ygpu);CHKERRQ(ierr);
	 	ierr = VecNorm(ysggpu,NORM_2,&normsggpu); 
		
		printf("SG-GPU Norm         = %.6f\n",normsggpu);

		if(normsggpu > normdiff)
			printf("SG GPU Test Failed\n");
		else 
			printf("SG GPU Test Passed\n");

//	}
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Cleaning
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = VecDestroy(&x);CHKERRQ(ierr); 
	ierr = VecDestroy(&y);CHKERRQ(ierr);
  	ierr = VecDestroy(&ysg);CHKERRQ(ierr);
  	ierr = VecDestroy(&ysgomp);CHKERRQ(ierr);
	ierr = MatDestroy(&mat);CHKERRQ(ierr);
	ierr = MatDestroy(&matsg);CHKERRQ(ierr);
	ierr = VecDestroy(&ygpu);CHKERRQ(ierr);
	ierr = VecDestroy(&ysggpu);CHKERRQ(ierr);
	ierr = MatDestroy(&matgpu);CHKERRQ(ierr);
	ierr = MatDestroy(&matsggpu);CHKERRQ(ierr);

  	ierr = PetscFinalize();
  return 0;
}
