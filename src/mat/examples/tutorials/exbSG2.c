/* Program usage:  mpiexec ex1 [-help] for all PETSc options
*/
static char help[] = "Simple program to test the performance of matmult for various matrix implementations. Enable appropriate flags to check an implementation.( CSR,SG,OMP,GPU). Run time options: [-n] [-m] [-p] [-dim] [-REF] [-info 1 for more info] Note: All of these flags are required except for info\n\n";

#include<sys/time.h>
#include "../../impls/blockstructgrid/seq/matblockstructgrid.h"
#include <petscksp.h> // this includes all the below headers
//#include<petscsys.h >//      	- base PETSc routines   petscvec.h - vectors
//#include<petscmat.h>// 	- matrices
// #include<petscis.h>//     	- index sets            petscksp.h - Krylov subspace methods
//#include<petscviewer.h>// 	- viewers               petscpc.h  - preconditioners
#define PAPI

#ifdef PAPI
#include"papi.h"
#define NUM_EVENTS 4
#endif

PetscReal normdiff = 1.0e-6;

unsigned int seed = 1;
 double simple_rand() {
         seed = (1103515245*seed+12345)%4294967296;
         return (1000.0*seed)/4294967296;
 }

double rtclock() {
  struct timezone tzp;
  struct timeval tp;
  gettimeofday (&tp, &tzp);
  return (1.0*tp.tv_sec + tp.tv_usec*1.0e-6);
}

#define CSR
#define SG
#define OMP
//#define GPU
//

#ifdef OMP
extern int OPENMPB; 
#endif

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

#ifdef PAPI
unsigned int Events[NUM_EVENTS] = {PAPI_LD_INS,PAPI_L1_DCM,PAPI_L2_TCM,PAPI_L3_TCM};
long_long values[NUM_EVENTS], sgvalues[NUM_EVENTS];
int e;
#endif
  	PetscMPIInt    size;
  	PetscErrorCode ierr;

	PetscInitialize(&argc,&args,(char *)0,help);
  	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

    /*Default values*/
	PetscInt m=1,n=1,p=1,dim=3,dof=1;
	PetscInt nos;
	PetscInt info=0;
	long REP=1;

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Command line arguments for m,n,p,dim, dof, REP 
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
	printf("Reps = %ld, m=%d,n=%d,p=%d, dim=%d, dof=%d\n",REP,m,n,p,dim, dof);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Set dims[], nz and nos.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscScalar    *vals;
	PetscInt       i, nz, *dims, *starts,  *cols;
 	dims = malloc(sizeof(PetscInt)*dim);
  	starts = malloc(sizeof(PetscInt)*dim);
	nos = (dim*2+1)*(dof*2-1);
	dims[0]=m;dims[1]=n;dims[2]=p;
	nz=m*n*p*dof;

	Vec x;
	double start,end;
  	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  	ierr = VecSetSizes(x,PETSC_DECIDE,nz);CHKERRQ(ierr);
  	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecSetRandom(x,PETSC_NULL);

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
		printf("st=%d, idx=%d, idy=%d, idz=%d, xval=%d\n",st,idx[st],idy[st],idz[st],xval[st]);
        }
	printf("nost=%d, nos=%d, nz=%d\n",nost,nos,nz);
/*CSR*/
#ifdef CSR
	Vec  	y;
  	Mat  	mat;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Create matrices.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatCreate(PETSC_COMM_WORLD,&mat);CHKERRQ(ierr);
  	ierr = MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,nz,nz);CHKERRQ(ierr);
  	MatSetType(mat,MATSEQAIJ);
	MatSeqAIJSetPreallocation(mat,nos,PETSC_NULL);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Set values into input vector and matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	//ierr = VecSet(x,one);CHKERRQ(ierr);//this can be modified such that x holds random values
	for(i=0;i<m*n*p;i++)
	{
		for(st=0;st<nost;st++)
        	{
			j=xval[st]+i;
			if(j>=0 && j<m*n*p)
			{
				for(l=0;l<dof;l++)
				{
					for(k=0;k<dof;k++)	
					{
						vals[k] = simple_rand();
						cols[k] = j*dof+k; 
					}
					rowval = i*dof+l;
   					ierr = MatSetValues(mat,1,&rowval,dof,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
				}
			}
		}
        }

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     AssemblyBegin/End as values can still remain in Cache
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Compute solution vectors and Performance test.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef PAPI
if (PAPI_start_counters((int *)Events, NUM_EVENTS) != PAPI_OK)
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
	printf("CSR Events[%d]= %lld\n",e,values[e]);
#endif

	printf("\nCSR :\n");
	PetscInt nz1 = ((nost*m*n*p) - (2*(1+m+m*n))) *dof*dof;
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nz1)/((end-start)*1024*1024*1024)); 
	fflush(stdout);	
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Cleaning
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = VecDestroy(&y);CHKERRQ(ierr);
	ierr = MatDestroy(&mat);CHKERRQ(ierr);

#endif

/*SG(AVX)*/
#ifdef SG
	Vec     ysg;
  	ierr = VecDuplicate(x,&ysg);CHKERRQ(ierr);
	Mat     matsg;
	PetscInt bcols;
	PetscScalar * bvals;
	PetscInt nop = ((nost*m*n*p) - (2*(1+m+m*n))) *dof*dof;
  	ierr = MatCreate(PETSC_COMM_WORLD,&matsg);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsg,nz,nz,nz,nz);CHKERRQ(ierr);
	MatSetType(matsg,MATBLOCKSTRUCTGRID);
  	ierr = MatSetStencil(matsg,dim,dims,starts,dof);CHKERRQ(ierr);
	bvals = malloc (sizeof(PetscScalar)*dof*dof);
	MatSetUpPreallocation_SeqBSG(matsg);
	for(i=0;i<m*n*p;i++)
	{
		for(st=0;st<nost;st++)
        	{
			j=xval[st]+i;
			if(j>=0 && j<m*n*p)
			{
				for(l=0;l<dof;l++)
				{
					for(k=0;k<dof;k++)	
					{
						bvals[l*dof+k] = simple_rand();
					}
				}
				rowval = i;
				bcols = j;
   				ierr = MatSetValues(matsg,1,&rowval,1,&bcols,bvals,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
        }
  	ierr = MatAssemblyBegin(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("Sg start error\n");
#endif
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysg);CHKERRQ(ierr);
	end = rtclock();
#ifdef PAPI
if (PAPI_read_counters(sgvalues, NUM_EVENTS) != PAPI_OK)
	printf("sg stop error\n");
for(e=0;e<NUM_EVENTS;e++)
	printf("SG Events[%d]= %lld\n",e,sgvalues[e]);
#endif
	
	printf("\nSG -AVX with Padding(original):\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nop)/((end-start)*1024*1024*1024)); 
	
	ierr = VecDestroy(&ysg);CHKERRQ(ierr);

#ifdef OMP
	Vec ysgomp;
	ierr = VecDuplicate(x,&ysgomp);CHKERRQ(ierr);
	//SG (AVX+OPENMP)
	OPENMPB=1;
	for(k=2;k<3;k++){
#ifdef PAPI
if (PAPI_read_counters(values, NUM_EVENTS) != PAPI_OK)
	printf("Sg start error\n");
#endif
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysgomp);CHKERRQ(ierr);
	end = rtclock();
#ifdef PAPI
if (PAPI_read_counters(sgvalues, NUM_EVENTS) != PAPI_OK)
	printf("sg stop error\n");
for(e=0;e<NUM_EVENTS;e++)
	printf("SG Events[%d]= %lld\n",e,sgvalues[e]);
#endif
	
	printf("\nSG - AVX + OPENMP :\n");
	printf("Threads= %d , Time =%.3f\n GFLOPS= %.3f\n",k,end-start,((long)REP*2*nop)/((end-start)*1024*1024*1024)); 
	
	}
  	ierr = VecDestroy(&ysgomp);CHKERRQ(ierr);
#endif

	ierr = MatDestroy(&matsg);CHKERRQ(ierr);
	free(bvals);
#endif

  	ierr = VecDestroy(&x);CHKERRQ(ierr); 
	//ierr = MatDestroy(&matip);CHKERRQ(ierr);

 	free(dims);
  	free(starts);
	free(xval);
	free(cols);
	free(vals);
  	
	ierr = PetscFinalize();
	return 0;
}
