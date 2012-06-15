/* Program usage:  mpiexec ex1 [-help] for all PETSc options
*/
static char help[] = "Simple program which does matrix vector multiplication using the default format aij and other formats namely, structgrid(avx, avx+openmp), aijcusp and structgridgpu. The resulting vectors are compared for consistency (when REP=1)and also tested for performance. Enable appropriate flags to check an implementation.( CSR,SG,OMP,GPU). Run time options: [-n] [-m] [-p] [-dim] [-REF] [-info 1 for more info]. All of these flags are required except for info. Note: It is preferable to use exSG2 while checking performance (especially for big inputs) as it has less memory footprint.\n\n";

#include <petscconf.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include "../../impls/aij/seq/aij.h"
#include<sys/time.h>
#include "../../impls/structgrid/matstructgrid.h"
#include <petscksp.h> // this includes all the below headers
//#include<petscsys.h >//      	- base PETSc routines   petscvec.h - vectors
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
PetscReal normdiff = 1.0e-6;
long REP=1;
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
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  	ierr = VecSetSizes(x,PETSC_DECIDE,nz);CHKERRQ(ierr);
  	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
  	ierr = VecDuplicate(x,&ysg);CHKERRQ(ierr);
#ifdef OMP
	ierr = VecDuplicate(x,&ysgomp);CHKERRQ(ierr);
#endif
#ifdef GPU
	ierr = VecDuplicate(x,&ysggpu);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&ygpu);CHKERRQ(ierr);
#endif
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create matrices.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatCreate(PETSC_COMM_WORLD,&mat);CHKERRQ(ierr);
  	ierr = MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,nz,nz);CHKERRQ(ierr);
  	MatSetType(mat,MATSEQAIJ);
	MatSeqAIJSetPreallocation(mat,nos,PETSC_NULL);
  	
	ierr   = MatCreate(PETSC_COMM_WORLD,&matsg);CHKERRQ(ierr);
  	ierr   = MatSetSizes(matsg,nz,nz,nz,nz);CHKERRQ(ierr);
	MatSetType(matsg,MATSTRUCTGRID);
#ifdef GPU 
	ierr   = MatCreate(PETSC_COMM_WORLD,&matsggpu);CHKERRQ(ierr);
  	ierr   = MatSetSizes(matsggpu,nz,nz,nz,nz);CHKERRQ(ierr);
  	ierr   = MatCreate(PETSC_COMM_WORLD,&matgpu);CHKERRQ(ierr);
  	ierr   = MatSetSizes(matgpu,nz,nz,nz,nz);CHKERRQ(ierr);
	MatSetType(matsggpu,MATSEQSGGPU);
	MatSetType(matgpu,MATSEQAIJCUSP);
        MatSetUp(matgpu);
        MatSetUp(matsggpu);
#endif
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set stencils for Structgrid -matsg
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	starts = malloc(sizeof(PetscInt)*dim);
  	ierr   = MatSetStencil(matsg,dim,dims,starts,dof);CHKERRQ(ierr);
	MatSetUpPreallocation_SeqSG(matsg);
#ifdef GPU
	ierr   = MatSetStencil(matsggpu,dim,dims,starts,dof);CHKERRQ(ierr);
	//MatSetUpPreallocation_SeqSG(matsggpu);
#endif
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set values into input vector and matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//  	ierr   = VecSet(x,1.0);CHKERRQ(ierr);//this can be modified such that x holds random values
	ierr   = VecSetRandom(x,PETSC_NULL);

	cols = malloc(sizeof(PetscInt)*dof);
	vals = malloc(sizeof(PetscScalar)*dof);
	
	PetscInt j,k,l,st;
        PetscInt lda2 = m*n;
        PetscInt lda3 = m;
	PetscInt rowval = 0;

	PetscInt nost;
	nost =  (dim*2+1);
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
   					ierr = MatSetValues(mat,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   					ierr = MatSetValues(matsg,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);

#ifdef GPU
					ierr = MatSetValues(matsggpu,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
                                        ierr = MatSetValues(matgpu,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
#endif
				}
			}
		}
        }

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     AssemblyBegin/End as values can still remain in Cache
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
#ifdef GPU
  	ierr = MatAssemblyBegin(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matgpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matgpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
#endif
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
	for(i=0;i<dim;i++)
	printf("dims[%d] = %d\n", i,dims[i]);
	printf("Rep=%d\n",REP);

	//CSR
#ifdef PAPI
if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
#endif
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(mat,x,y);CHKERRQ(ierr);
	end = rtclock();
#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
for(e=0;e<NUM_EVENTS;e++)
	printf("Events[%d]= %lld\n",e,values[e]);
#endif
	printf("\nCSR :\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 
	fflush(stdout);	
	//SG (AVX)
	OPENMP = 0;
#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
#endif
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysg);CHKERRQ(ierr);
	end = rtclock();
#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
for(e=0;e<NUM_EVENTS;e++)
	printf("Events[%d]= %lld\n",e,values[e]);
#endif
	printf("\nSG -AVX with Padding(original):\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 

#ifdef OMP
	//SG (AVX+OPENMP)
	OPENMP=1;
	
	omp_set_num_threads(NUM_THREADS);
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysgomp);CHKERRQ(ierr);
	end = rtclock();
	printf("\nSG - AVX + OPENMP :\n");
	printf("Threads= %d, Time =%.3f\n GFLOPS= %.3f\n",k,end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 
#endif

#ifdef GPU
	//CSR GPU
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matgpu,x,ygpu);CHKERRQ(ierr);
	end = rtclock();
	printf("\nCSR - GPU:\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 
         

	//SG (GPU)
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsggpu,x,ysggpu);CHKERRQ(ierr);
	end = rtclock();
	printf("\nSG - GPU:\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 
#endif

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	if(info){
	printf("\nOutput:\n");
	printf("Y - CSR:\n");
  	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	printf("Y - Structgrid AVX:\n");
  	ierr = VecView(ysg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#ifdef OMP	
	printf("Y - Structgrid AVX + OPENMP:\n");
  	ierr = VecView(ysgomp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
#ifdef GPU
	printf("Y - CSR CUSP(GPU):\n");
        ierr = VecView(ygpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	printf("Y - Structgrid GPU:\n");
        ierr = VecView(ysggpu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif	
	}


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Correctness test
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	if(REP==1){//correctness test only when REP =1, can be disabled if required
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
#ifdef OMP
		ierr = VecAXPY(ysgomp,-1,y);CHKERRQ(ierr);
//		ierr = VecView(ysgomp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	 	ierr = VecNorm(ysgomp,NORM_2,&normsgomp); 
		printf("SG(AVX+OPENMP) Norm = %.6f\n",normsgomp);
		if(normsgomp > normdiff)
			printf("SG AVX+Openmp Test Failed\n");
		else 
			printf("SG AVX+Openmp Test Passed\n");
#endif
#ifdef GPU
		ierr = VecAXPY(ysggpu,-1,ygpu);CHKERRQ(ierr);
	 	ierr = VecNorm(ysggpu,NORM_2,&normsggpu); 
		
		printf("SG-GPU Norm         = %.6f\n",normsggpu);

		if(normsggpu > normdiff)
			printf("SG GPU Test Failed\n");
		else 
			printf("SG GPU Test Passed\n");
#endif
	}
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Cleaning
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = VecDestroy(&x);CHKERRQ(ierr); 
	ierr = VecDestroy(&y);CHKERRQ(ierr);
  	ierr = VecDestroy(&ysg);CHKERRQ(ierr);
#ifdef OMP  
	ierr = VecDestroy(&ysgomp);CHKERRQ(ierr);
#endif
	ierr = MatDestroy(&mat);CHKERRQ(ierr);
	ierr = MatDestroy(&matsg);CHKERRQ(ierr);
#ifdef GPU
        ierr = VecDestroy(&ygpu);CHKERRQ(ierr);
	ierr = VecDestroy(&ysggpu);CHKERRQ(ierr);
	ierr = MatDestroy(&matgpu);CHKERRQ(ierr);
	ierr = MatDestroy(&matsggpu);CHKERRQ(ierr);
#endif
 	free(dims);
  	free(starts);
	free(cols);
	free(vals);
  	ierr = PetscFinalize();
  return 0;
}
