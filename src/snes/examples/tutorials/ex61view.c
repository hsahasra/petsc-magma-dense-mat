static const char help[] = "Loads data generated by ex61 and VTK file suitable for Paraview or Visit.\n\n";

/*

      ./ex61gen  -random_seed <integer>                            :   creates ex61.random.<integer>
      ./ex61 [-random_seed <integer>] -graphicsfile other options  :   creates ex61.data[.<integer>]
      ./ex61view [-random_seed <integer>]                          :   creates ex61.data[.<integer>]_cv_%d.vtk  ex61.data[.<integer>]_eta_%d.vtk

    where [ ] indicates the optional seed value;  <integer> the seed value provided and %d the timestep

    The resulting .vtk files may be loaded into ParaView or ViSIT

*/

#include "petscdmda.h"


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  PetscErrorCode      ierr;
  Vec                 U,cv,eta;
  DM                  da,da2;
  PetscViewer         viewer,view_vtk_cv,view_vtk_eta;
  char                filename[PETSC_MAX_PATH_LEN],cv_filename[PETSC_MAX_PATH_LEN],eta_filename[PETSC_MAX_PATH_LEN];
  PetscBool           flg,sflg = PETSC_FALSE;
  PetscInt            i,n=10000;
  PetscInt            seed;

  PetscInitialize(&argc,&argv, (char*)0, help);
  ierr = PetscOptionsSetValue("-viewer_binary_skip_info","true");CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-f",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr = PetscOptionsGetInt(PETSC_NULL,"-random_seed",&seed,&sflg);CHKERRQ(ierr);
    if (!sflg) {
      ierr = PetscStrcpy(filename,"ex61.data");CHKERRQ(ierr);
    } else {
      sprintf(filename,"ex61.data.%d",seed);
    }
  }

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);

  /* Get physics and time parameters */
  ierr = DMCreate(PETSC_COMM_WORLD,&da);CHKERRQ(ierr);
  ierr = DMLoad(da,viewer);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&U);CHKERRQ(ierr);
  ierr = DMDAGetReducedDA(da,1,&da2);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da2,&cv);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da2,&eta);CHKERRQ(ierr);

  for (i=0; i<n; i++) {
    /* when this fails it simply means the file is finished */
    ierr = VecLoad(U,viewer);CHKERRQ(ierr);
    ierr = VecStrideGather(U,1,cv,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecStrideGather(U,4,eta,INSERT_VALUES);CHKERRQ(ierr);
    sprintf(cv_filename,"%s_cv_%d.vtk",filename,i);
    sprintf(eta_filename,"%s_eta_%d.vtk",filename,i);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,cv_filename,&view_vtk_cv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,eta_filename,&view_vtk_eta);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(view_vtk_cv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(view_vtk_eta, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
    ierr = DMView(da2,view_vtk_cv);CHKERRQ(ierr);
    ierr = DMView(da2,view_vtk_eta);CHKERRQ(ierr);
    ierr = VecView(cv,view_vtk_cv);CHKERRQ(ierr);
    ierr = VecView(eta,view_vtk_eta);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&view_vtk_cv);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&view_vtk_eta);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&cv);CHKERRQ(ierr);
  ierr = VecDestroy(&eta);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = DMDestroy(&da2);CHKERRQ(ierr);
  PetscFinalize();
  return 0;
}
