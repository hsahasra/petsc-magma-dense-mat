static char help[] = "Example for PetscFileUpload()";
#include <petscsys.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode  ierr;

  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = PetscFileUpload(PETSC_COMM_WORLD,"ex11.c","login.mcs.anl.gov","/home/bsmith/public_html/",NULL);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}

