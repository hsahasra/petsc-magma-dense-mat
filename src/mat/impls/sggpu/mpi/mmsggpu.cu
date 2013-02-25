#include <../src/mat/impls/sggpu/mpi/mpisggpu.h>
#include <petscblaslapack.h>


//This function is a copy of MatSetUpMultiply_Dense(). Need to edit this function to make it work for MPISGGPU. - Chekuri

#undef __FUNCT__
#define __FUNCT__ "MatSetUpMultiply_MPISGGPU"
PetscErrorCode MatSetUpMultiply_MPISGGPU(Mat mat)
{
  Mat_MPISGGPU *sggpu = (Mat_MPISGGPU*)mat->data;
  PetscErrorCode ierr;
  IS           from,to;
  Vec          gvec;
  PetscInt numprocs;

  PetscFunctionBegin;
  
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);
  /* Create local vector that is used to scatter into */
//  ierr = VecCreateSeq(PETSC_COMM_SELF,mat->cmap->N,&mat->lvec);CHKERRQ(ierr);

  /* Create temporary index set for building scatter gather */
//  ierr = ISCreateStride(((PetscObject)mat)->comm,mat->cmap->N,0,1,&from);CHKERRQ(ierr);
//  ierr = ISCreateStride(PETSC_COMM_SELF,mat->cmap->N,0,1,&to);CHKERRQ(ierr);

  /* Create temporary global vector to generate scatter context */
  /* n    = mdn->cowners[mdn->rank+1] - mdn->cowners[mdn->rank]; */

//  ierr = VecCreateMPIWithArray(((PetscObject)mat)->comm,1,((mat->cmap->N)/numprocs),mat->cmap->N,PETSC_NULL,&gvec);CHKERRQ(ierr);

  /* Generate the scatter context */
//  ierr = VecScatterCreate(gvec,from,sggpu->lvec,to,&sggpu->Mvctx);CHKERRQ(ierr);
//  ierr = PetscLogObjectParent(mat,sggpu->Mvctx);CHKERRQ(ierr);
//  ierr = PetscLogObjectParent(mat,sggpu->lvec);CHKERRQ(ierr);
//  ierr = PetscLogObjectParent(mat,from);CHKERRQ(ierr);
//  ierr = PetscLogObjectParent(mat,to);CHKERRQ(ierr);
//  ierr = PetscLogObjectParent(mat,gvec);CHKERRQ(ierr);

//  ierr = ISDestroy(&to);CHKERRQ(ierr);
//  ierr = ISDestroy(&from);CHKERRQ(ierr);
//  ierr = VecDestroy(&gvec);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}
