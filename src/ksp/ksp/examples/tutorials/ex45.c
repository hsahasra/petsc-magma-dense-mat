
/*
Laplacian in 3D. Modeled by the partial differential equation

   - Laplacian u = 1,0 < x,y,z < 1,

with boundary conditions

   u = 1 for x = 0, x = 1, y = 0, y = 1, z = 0, z = 1.

   This uses multigrid to solve the linear system

   The same as ex22.c except it does not use DMMG, it uses its replacement.
   See src/snes/examples/tutorials/ex50.c

   Can also be run with -pc_type exotic -ksp_type fgmres

*/

static char help[] = "Solves 3D Laplacian using multigrid.\n\n";

#include <petscksp.h>
#include <petscdmda.h>

extern PetscErrorCode ComputeMatrix(DM,Vec,Mat,Mat,MatStructure*);
extern PetscErrorCode ComputeRHS(DM,Vec,Vec);
extern PetscErrorCode ComputeInitialGuess(DM,Vec);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  KSP            ksp;
  PetscReal      norm;
  DM             da;
  Vec            x,b,r;
  Mat            A;

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,-7,-7,-7,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,0,&da);CHKERRQ(ierr);  
  ierr = DMSetInitialGuess(da,ComputeInitialGuess);CHKERRQ(ierr);
  ierr = DMSetFunction(da,ComputeRHS);CHKERRQ(ierr);
  ierr = DMSetJacobian(da,ComputeMatrix);CHKERRQ(ierr);
  ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);
  /*  ierr = KSPSetDMActive(ksp,PETSC_FALSE);CHKERRQ(ierr);*/
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = KSPGetSolution(ksp,&x);CHKERRQ(ierr);
  ierr = KSPGetRhs(ksp,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&r);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp,&A,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  ierr = MatMult(A,x,r);CHKERRQ(ierr);
  ierr = VecAXPY(r,-1.0,b);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Residual norm %G\n",norm);CHKERRQ(ierr); 

  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscFinalize();

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(DM dm,Vec x,Vec b)
{
  PetscErrorCode ierr;
  PetscInt       mx,my,mz;
  PetscScalar    h;
  PetscScalar    ***xx,***bb;
  PetscInt       i,j,k,xm,ym,zm,xs,ys,zs;
  PetscScalar    Hx,Hy,Hz,HxHydHz,HyHzdHx,HxHzdHy;

  PetscFunctionBegin;
  ierr = DMDAGetInfo(dm,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  h    = 1.0/((mx-1)*(my-1)*(mz-1));
  ierr = VecSet(b,h);CHKERRQ(ierr);

  if (x) {
    ierr = DMDAVecGetArray(dm,x,&xx);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm,b,&bb);CHKERRQ(ierr);

  ierr = DMDAGetInfo(dm,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);  
  Hx = 1.0 / (PetscReal)(mx-1); Hy = 1.0 / (PetscReal)(my-1); Hz = 1.0 / (PetscReal)(mz-1);
  HxHydHz = Hx*Hy/Hz; HxHzdHy = Hx*Hz/Hy; HyHzdHx = Hy*Hz/Hx;
  ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  for (k=zs; k<zs+zm; k++){
    for (j=ys; j<ys+ym; j++){
      for(i=xs; i<xs+xm; i++){
	if (!(i==0 || j==0 || k==0 || i==mx-1 || j==my-1 || k==mz-1)){
          bb[k][j][i] -= +HxHydHz*xx[k-1][j][i] + HxHzdHy*xx[k][j-1][i] + HyHzdHx*xx[k][j][i-1] + HyHzdHx*xx[k][j][i+1] + HxHzdHy*xx[k][j+1][i] + HxHydHz*xx[k+1][j][i];
        }
        bb[k][j][i] += 2.0*(HxHydHz + HxHzdHy + HyHzdHx)*xx[k][j][i];
      }
    }
  }
    ierr = DMDAVecRestoreArray(dm,x,&xx);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(dm,b,&bb);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
    
#undef __FUNCT__
#define __FUNCT__ "ComputeInitialGuess"
PetscErrorCode ComputeInitialGuess(DM dm,Vec b)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(b,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(DM dm,Vec x,Mat jac,Mat B,MatStructure *stflg)
{
  DM             da = dm;
  PetscErrorCode ierr;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  PetscScalar    v[7],Hx,Hy,Hz,HxHydHz,HyHzdHx,HxHzdHy;
  MatStencil     row,col[7];

  ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);  
  Hx = 1.0 / (PetscReal)(mx-1); Hy = 1.0 / (PetscReal)(my-1); Hz = 1.0 / (PetscReal)(mz-1);
  HxHydHz = Hx*Hy/Hz; HxHzdHy = Hx*Hz/Hy; HyHzdHx = Hy*Hz/Hx;
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  for (k=zs; k<zs+zm; k++){
    for (j=ys; j<ys+ym; j++){
      for(i=xs; i<xs+xm; i++){
        row.i = i; row.j = j; row.k = k;
	if (i==0 || j==0 || k==0 || i==mx-1 || j==my-1 || k==mz-1){
          v[0] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx);
	  ierr = MatSetValuesStencil(B,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
	} else {
	  v[0] = -HxHydHz;col[0].i = i; col[0].j = j; col[0].k = k-1;
	  v[1] = -HxHzdHy;col[1].i = i; col[1].j = j-1; col[1].k = k;
	  v[2] = -HyHzdHx;col[2].i = i-1; col[2].j = j; col[2].k = k;
	  v[3] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx);col[3].i = row.i; col[3].j = row.j; col[3].k = row.k;
	  v[4] = -HyHzdHx;col[4].i = i+1; col[4].j = j; col[4].k = k;
	  v[5] = -HxHzdHy;col[5].i = i; col[5].j = j+1; col[5].k = k;
	  v[6] = -HxHydHz;col[6].i = i; col[6].j = j; col[6].k = k+1;
	  ierr = MatSetValuesStencil(B,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  *stflg = SAME_NONZERO_PATTERN;
  return 0;
}


