/* Program usage:  mpiexec ex1 [-help] for all PETSc options
*/
static char help[] = "Simple program to test the performance of matmult for various matrix implementations. Enable appropriate flags to check an implementation.( CSR,SG,OMP,GPU). Run time options: [-n] [-m] [-p] [-dim] [-REF] [-info 1 for more info]\n\n";

#include<sys/time.h>
#include "../../impls/structgrid/matstructgrid.h"
#include <petscksp.h> // this includes all the below headers
//#include<petscsys.h >//      	- base PETSc routines   petscvec.h - vectors
//#include<petscmat.h>// 	- matrices
// #include<petscis.h>//     	- index sets            petscksp.h - Krylov subspace methods
//#include<petscviewer.h>// 	- viewers               petscpc.h  - preconditioners

PetscReal normdiff = 1.0e-6;
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

#define CSR
#define SG
//#define OMP
//#define GPU

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

OPENMP=0;

  	PetscMPIInt    size;
  	PetscErrorCode ierr;

	PetscInitialize(&argc,&args,(char *)0,help);
  	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

    /*Default values*/
	PetscInt m=2,n=2,p=2,dim=3,dof=1;
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
	printf("Reps = %d, m=%d,n=%d,p=%d, dim=%d, dof=%d\n",REP,m,n,p,dim, dof);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Set dims[], nz and nos.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscScalar    *vals;
	PetscInt       i, nz, *dims, *starts, *rows, *cols;
 	dims = malloc(sizeof(PetscInt)*dim);
  	starts = malloc(sizeof(PetscInt)*dim);
	nos = (dim*2+1)*(2*dof-1);
	dims[0]=m;dims[1]=n;dims[2]=p;
	nz=m*n*p*dof;

	/*Mat     matip;
	ierr = MatCreate(PETSC_COMM_WORLD,&matip);CHKERRQ(ierr);
  	ierr = MatSetSizes(matip,nz,nz,nz,nz);CHKERRQ(ierr);
	ierr = MatSetType(matip,MATSTRUCTGRID);//CHKERRQ(ierr);
  	ierr = MatSetStencil(matip,dim,dims,starts,dof);CHKERRQ(ierr);
	Mat_SeqSG * sg = (Mat_SeqSG*) matip->data;
	*/
	
	PetscInt *idx =  malloc (sizeof(PetscInt)*nos);
	PetscInt *idy =  malloc (sizeof(PetscInt)*nos);
	PetscInt *idz =  malloc (sizeof(PetscInt)*nos);
	int cnt=0;	
	idx[cnt] = 0; idy[cnt]=0; idz[cnt++]= 0;
	for(i=1;i<dof;i++)
	{
		idx[cnt] = i; idy[cnt] = 0; idz[cnt++] = 0;
		idx[cnt] = -i; idy[cnt] = 0; idz[cnt++] = 0;
	}
	if(dim>0)
	{	
		idx[cnt] = -dof; idy[cnt] = 0; idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			idx[cnt] = -dof+i; idy[cnt] = 0; idz[cnt++] = 0;
			idx[cnt] = -dof-i; idy[cnt] = 0; idz[cnt++] = 0;
		}
			
		idx[cnt] = dof; idy[cnt] = 0; idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			idx[cnt] = dof+i; idy[cnt] = 0; idz[cnt++] = 0;
			idx[cnt] = dof-i; idy[cnt] = 0; idz[cnt++] = 0;
		}	
	}
	if(dim>1)
	{
		idx[cnt] = 0; idy[cnt] = -1; idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			idx[cnt] = i; idy[cnt] = -1; idz[cnt++] = 0;
			idx[cnt] = -i; idy[cnt] = -1; idz[cnt++] = 0;
		}
			
		idx[cnt] = 0; idy[cnt] = 1; idz[cnt++] = 0;
		for(i=1;i<dof;i++)
		{
			idx[cnt] = i; idy[cnt] = 1; idz[cnt++] = 0;
			idx[cnt] = -i; idy[cnt] = 1; idz[cnt++] = 0;
		}	
	}
	if(dim>2)
	{
		idx[cnt] = 0; idy[cnt] = 0; idz[cnt++] = -1;
		for(i=1;i<dof;i++)
		{
			idx[cnt] = i; idy[cnt] = 0; idz[cnt++] = -1;
			idx[cnt] = -i; idy[cnt] = 0; idz[cnt++] = -1;
		}
		idx[cnt] = 0; idy[cnt] = 0; idz[cnt++] = 1;
	}
	PetscInt k,l;
        PetscInt lda1 = m*n*p*dof;
        PetscInt lda2 = m*n*dof;
        PetscInt lda3 = m*dof;
	PetscInt *offset = malloc(sizeof(PetscInt)*nos);
	PetscInt *xval = malloc(sizeof(PetscInt)*nos);
	cols = malloc(sizeof(PetscInt)*nos);
	vals = malloc(sizeof(PetscScalar)*nos);

	for(l=0;l<nos;l++)
        {
                offset[l] = l*lda1;
		printf("idx[%d]=%d\n",l,idx[l]);
		fflush(stdout);
                xval[l] = idx[l] + idy[l]*lda3 + idz[l]*lda2;
        }
/*	
	for(l=0;l<nos;l++)
        {
                offset[l] = l*lda1;
		printf("sg->idx[%d]=%d\n",l,sg->idx[l]);
		fflush(stdout);
                xval[l] = sg->idx[l] + sg->idy[l]*lda3 + sg->idz[l]*lda2;
        }
*/	
	PetscInt count;
	Vec x;
	double start,end;
  	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  	ierr = VecSetSizes(x,PETSC_DECIDE,nz);CHKERRQ(ierr);
  	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecSetRandom(x,PETSC_NULL);

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

	for(i=0;i<nz;i++)
	{
		count=0;
		for(l=0;l<nos;l++)
        	{
                        vals[count] = simple_rand();
			if(xval[l]+i<nz)
				cols[count++] =  (xval[l]+i);    
		
		}
   		ierr = MatSetValues(mat,1,&i,count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        }
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     AssemblyBegin/End as values can still remain in Cache
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  	ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Compute solution vectors and Performance test.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(mat,x,y);CHKERRQ(ierr);
	end = rtclock();

	printf("\nCSR :\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 
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
  	ierr = MatCreate(PETSC_COMM_WORLD,&matsg);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsg,nz,nz,nz,nz);CHKERRQ(ierr);
	MatSetType(matsg,MATSTRUCTGRID);
  	ierr = MatSetStencil(matsg,dim,dims,starts,dof);CHKERRQ(ierr);
	MatSetUpPreallocation_SeqSG(matsg);
	for(i=0;i<nz;i++)
	{
		count=0;
		for(l=0;l<nos;l++)
        	{
                        vals[count] = simple_rand();
			if(xval[l]+i<nz)
				cols[count++] =  (xval[l]+i);    
		
		}
   		ierr = MatSetValues(matsg,1,&i,count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        }
  	ierr = MatAssemblyBegin(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsg,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		

	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysg);CHKERRQ(ierr);
	end = rtclock();
	
	printf("\nSG -AVX with Padding(original):\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*nos*nz)/((end-start)*1024*1024*1024)); 
	
	ierr = VecDestroy(&ysg);CHKERRQ(ierr);

/*SG(AVX+OMP)*/
#ifdef OMP
	Vec ysgomp;
	ierr = VecDuplicate(x,&ysgomp);CHKERRQ(ierr);
	//SG (AVX+OPENMP)
	OPENMP=1;
	
	for(k=2;k<3;k++){
	omp_set_num_threads(k);
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsg,x,ysgomp);CHKERRQ(ierr);
	end = rtclock();
	printf("\nSG - AVX + OPENMP :\n");
	printf("Threads= %d, Time =%.3f\n GFLOPS= %.3f\n",k,end-start,((long)REP*2*sg->stpoints*sg->nz)/((end-start)*1024*1024*1024)); 
	}
  	ierr = VecDestroy(&ysgomp);CHKERRQ(ierr);
#endif
	ierr = MatDestroy(&matsg);CHKERRQ(ierr);
#endif

#ifdef GPU
/* GPU*/
	Vec  	ygpu, ysggpu;      
	Mat 	matgpu, matsggpu;           
	ierr = VecDuplicate(x,&ysggpu);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&ygpu);CHKERRQ(ierr);
  	ierr = MatCreate(PETSC_COMM_WORLD,&matsggpu);CHKERRQ(ierr);
  	ierr = MatSetSizes(matsggpu,nz,nz,nz,nz);CHKERRQ(ierr);
  	ierr = MatCreate(PETSC_COMM_WORLD,&matgpu);CHKERRQ(ierr);
  	ierr = MatSetSizes(matgpu,nz,nz,nz,nz);CHKERRQ(ierr);
	MatSetType(matsggpu,MATSTRUCTGRIDGPU);
	MatSetType(matgpu,MATSEQAIJCUSP);
	ierr = MatSetStencil(matsggpu,dim,dims,starts,dof);CHKERRQ(ierr);
	MatSetUpPreallocation_SeqSG(matsggpu);
	for(i=0;i<nz;i++)
	{
		count=0;
		for(l=0;l<nos;l++)
        	{
                        vals[count] = simple_rand();
			if(xval[l]+i<nz)
				cols[count++] =  (xval[l]+i);    
		
		}
   		ierr = MatSetValues(matgpu,1,&i,count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
   		ierr = MatSetValues(matsggpu,1,&i,count,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
        }
  
  	ierr = MatAssemblyBegin(matgpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matgpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(matsggpu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	//CSR GPU
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matgpu,x,ygpu);CHKERRQ(ierr);
	end = rtclock();
	printf("\nCSR - GPU:\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*sg->stpoints*sg->nz)/((end-start)*1024*1024*1024)); 

	//SG (GPU)
	start = rtclock();	
	for(i=0;i<REP;i++)
  		ierr = MatMult(matsggpu,x,ysggpu);CHKERRQ(ierr);
	end = rtclock();
	printf("\nSG - GPU:\n");
	printf("Time =%.3f\n GFLOPS= %.3f\n",end-start,((long)REP*2*sg->stpoints*sg->nz)/((end-start)*1024*1024*1024)); 

	ierr = VecDestroy(&ygpu);CHKERRQ(ierr);
	ierr = VecDestroy(&ysggpu);CHKERRQ(ierr);
	ierr = MatDestroy(&matgpu);CHKERRQ(ierr);
	ierr = MatDestroy(&matsggpu);CHKERRQ(ierr);
#endif
  	ierr = VecDestroy(&x);CHKERRQ(ierr); 
	//ierr = MatDestroy(&matip);CHKERRQ(ierr);

 	free(dims);
  	free(starts);
	free(offset);
	free(xval);
	free(cols);
	free(vals);
  	
	ierr = PetscFinalize();
	return 0;
}
