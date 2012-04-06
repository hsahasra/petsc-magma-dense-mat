/* Program usage:  mpiexec ex1 [-help] for all PETSc options
*/
static char help[] = "Simple program which does matrix vector multiplication using the default format aij and block structgrid. The resulting vectors are compared for consistency (when REP=1)and also tested for performance. Enable appropriate flags to check an implementation. Run time options: [-n] [-m] [-p] [-dim] [-REF] [-info 1 for more info]. All of these flags are required except for info. Note: It is preferable to use exSG2 while checking performance (especially for big inputs) as it has less memory footprint.\n\n";

#include<sys/time.h>
//#include "../../impls/blockstructgrid/matblockstructgrid.h"
#include <petscksp.h> // this includes all the below headers

PetscInt m=1,n=1,p=1,dim=3,dof=1;
PetscInt nos;
PetscInt info=0;
PetscReal normdiff = 1.0e-6;
long REP=1;

 double simple_rand() {
         int seed;
         seed = (1103515245*seed+12345)%4294967296;
         return 1.0;//(1000.0*seed)/4294967296;
 }

double rtclock() {
  struct timezone tzp;
  struct timeval tp;
  gettimeofday (&tp, &tzp);
  return (1.0*tp.tv_sec + tp.tv_usec*1.0e-6);
}

//#define OMP

#ifdef OMP
#include<omp.h>
extern int OPENMPB; 
#endif

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
	dof=1;

	Vec  	x, y, ybsg;
  	Mat  	mat, matbsg;
  	PetscErrorCode ierr;
  	PetscInt       i, nz=1, *dims, *starts, *cols;
	PetscScalar    *vals, one=1.0, *bvals;
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
  	ierr = VecDuplicate(x,&ybsg);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create matrices.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatCreate(PETSC_COMM_WORLD,&mat);CHKERRQ(ierr);
  	ierr = MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,nz,nz);CHKERRQ(ierr);
  	MatSetType(mat,MATSEQAIJ);
	//MatSeqAIJSetPreallocation(mat,nos,PETSC_NULL);
  	
	ierr = MatCreate(PETSC_COMM_WORLD,&matbsg);CHKERRQ(ierr);
  	ierr = MatSetSizes(matbsg,nz,nz,nz,nz);CHKERRQ(ierr);
	MatSetType(matbsg,MATBLOCKSTRUCTGRID);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set stencils for Structgrid -matsg
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	starts = malloc(sizeof(PetscInt)*dim);
  	ierr = MatSetStencil(matbsg,dim,dims,starts,dof);CHKERRQ(ierr);
	MatSetUpPreallocation_SeqBSG(matbsg);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set values into input vector and matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	//ierr = VecSet(x,1.0);CHKERRQ(ierr);//this can be modified such that x holds random values
	ierr = VecSetRandom(x,PETSC_NULL);

	cols = malloc(sizeof(PetscInt)*dof);
	vals = malloc(sizeof(PetscScalar)*dof);
	bvals = malloc(sizeof(PetscScalar)*dof*dof);
	
	PetscInt j,k,l,st;
        PetscInt lda2 = m*n;
        PetscInt lda3 = m;
	PetscInt rowval = 0;
	PetscInt bcols = 0;

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
	cnt = 0;
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
						bvals[l*dof+k] = vals[k];
						cols[k] = j*dof+k; 
					}
					rowval = i*dof+l;
   					ierr = MatSetValues(mat,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
				}
				bcols = j;
				rowval = i;
   				ierr = MatSetValues(matbsg,1,&rowval,1,&bcols,bvals,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
        }

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     AssemblyBegin/End as values can still remain in Cache
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matbsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matbsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if(info){
  	printf("\nInputs:\n");
  	//ierr = MatView(mat,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  	ierr = MatView(matbsg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
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
	//BSG (AVX)

#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
#endif
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matbsg,x,ybsg);CHKERRQ(ierr);
	end = rtclock();
#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("error\n");
for(e=0;e<NUM_EVENTS;e++)
	printf("Events[%d]= %lld\n",e,values[e]);
#endif
	printf("\nBlock SG -AVX:\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 

#ifdef OMP
	Vec ysgomp;
	ierr = VecDuplicate(x,&ysgomp);CHKERRQ(ierr);
	//SG (AVX+OPENMP)
	OPENMPB=1;
	for(k=4;k<5;k++){
	omp_set_num_threads(k);
#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("Sg start error\n");
#endif
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matbsg,x,ysgomp);CHKERRQ(ierr);
	end = rtclock();
#ifdef PAPI
if (PAPI_read_counters(sgvalues, NUM_EVENTS) != PAPI_OK)
	printf("sg stop error\n");
for(e=0;e<NUM_EVENTS;e++)
	printf("SG Events[%d]= %lld\n",e,sgvalues[e]);
#endif
	
	printf("\nSG - AVX + OPENMP :\n");
	printf("Threads= %d, Time =%.3f\n GFLOPS= %.3f\n",k,end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 
	
	}
#endif
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Print the input vector and matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	if(info){
	printf("\nOutput:\n");
  	ierr = MatView(matbsg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	printf("Y - CSR:\n");
//  	ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	printf("Y - Block Structgrid AVX:\n");
  //	ierr = VecView(ybsg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
#ifdef OMP
	printf("Y - Block Structgrid AVX + OMP:\n");
  //	ierr = VecView(ysgomp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
	}


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Correctness test
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	if(REP==1){//correctness test only when REP =1, can be disabled if required
		printf("\n\nCorrectness test : \n");
		PetscReal norm, normbsg;
		
		ierr = VecAXPY(ybsg,-1,y);CHKERRQ(ierr);
//		ierr = VecView(ysg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	 	ierr = VecNorm(ybsg,NORM_2,&normbsg);CHKERRQ(ierr);
		printf("BSG(AVX) Norm        = %.6f\n",normbsg);
		if(normbsg > normdiff)
			printf("BSG AVX Test Failed\n");
		else 
			printf("BSG AVX Test Passed\n");
  	
#ifdef OMP
		ierr = VecAXPY(ysgomp,-1,y);CHKERRQ(ierr);
	 	ierr = VecNorm(ysgomp,NORM_2,&normbsg);CHKERRQ(ierr);
		printf("BSG(AVX) + OMP Norm        = %.6f\n",normbsg);
		if(normbsg > normdiff)
			printf("BSG AVX + OMP Test Failed\n");
		else 
			printf("BSG AVX + OMP  Test Passed\n");
//		ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
	}
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Cleaning
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = VecDestroy(&x);CHKERRQ(ierr); 
	ierr = VecDestroy(&y);CHKERRQ(ierr);
  	ierr = VecDestroy(&ybsg);CHKERRQ(ierr);
#ifdef OMP
  	ierr = VecDestroy(&ysgomp);CHKERRQ(ierr);
#endif
 	free(dims);
  	free(starts);
	free(cols);
	free(vals);
	free(bvals);
  	ierr = PetscFinalize();
  return 0;
}
