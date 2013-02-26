#include <../src/mat/impls/sggpu/mpi/mpisggpu.h>
#include <petscblaslapack.h>


#undef __FUNCT__
#define __FUNCT__ "MatSetUpMultiply_MPISGGPU"
PetscErrorCode MatSetUpMultiply_MPISGGPU(Mat A)
{
  Mat_MPISGGPU *mat = (Mat_MPISGGPU*)A->data;
  Mat_SeqSGGPU *mat_seq = mat->mat_seq;
	
  PetscErrorCode ierr;
  IS           from,to;
  Vec          gvec;
  PetscInt numprocs;

  PetscFunctionBegin;
  
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);
  /* Create local vector that is used to scatter into */
  ierr = VecCreateSeq(PETSC_COMM_SELF,A->cmap->N,&mat->lvec);CHKERRQ(ierr);

  /* Create temporary index set for building scatter gather */
  ierr = ISCreateStride(((PetscObject)A)->comm,A->cmap->N,0,1,&from);CHKERRQ(ierr);
  ierr = ISCreateStride(PETSC_COMM_SELF,A->cmap->N,0,1,&to);CHKERRQ(ierr);

  /* Create temporary global vector to generate scatter context */
  /* n    = mdn->cowners[mdn->rank+1] - mdn->cowners[mdn->rank]; */

   ierr = VecCreateMPIWithArray(((PetscObject)A)->comm,1,((A->cmap->N)/mat->size),A->cmap->N,PETSC_NULL,&gvec);CHKERRQ(ierr);

  /* Generate the scatter context */
  ierr = VecScatterCreate(gvec,from,mat->lvec,to,&mat->Mvctx);CHKERRQ(ierr);
  ierr = PetscLogObjectParent(A,mat->Mvctx);CHKERRQ(ierr);
  ierr = PetscLogObjectParent(A,mat->lvec);CHKERRQ(ierr);
  ierr = PetscLogObjectParent(A,from);CHKERRQ(ierr);
  ierr = PetscLogObjectParent(A,to);CHKERRQ(ierr);
  ierr = PetscLogObjectParent(A,gvec);CHKERRQ(ierr);

  ierr = ISDestroy(&to);CHKERRQ(ierr);
  ierr = ISDestroy(&from);CHKERRQ(ierr);
  ierr = VecDestroy(&gvec);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}
