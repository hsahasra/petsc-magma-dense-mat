
/*
   This file contains routines for Parallel vector operations.
 */
#include <petscconf.h>
PETSC_CUDA_EXTERN_C_BEGIN
#include <../src/vec/vec/impls/mpi/pvecimpl.h>   /*I  "petscvec.h"   I*/
PETSC_CUDA_EXTERN_C_END
#include <../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h>
static cudaError_t ccs[16];
static cudaError_t cms[16];

extern MPI_Op VecMax_Local_Op;
extern MPI_Op VecMin_Local_Op;



#undef __FUNCT__
#define __FUNCT__ "VecDestroy_MPIGPU"
PetscErrorCode VecDestroy_MPIGPU(Vec v)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Destroy the stashes: note the order - so that the tags are freed properly */
  ierr = VecStashDestroy_Private(&v->bstash);CHKERRQ(ierr);
  ierr = VecStashDestroy_Private(&v->stash);CHKERRQ(ierr);
  ierr = VecDestroy_SeqGPU(v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecSetValues_MPIGPU"
PetscErrorCode VecSetValues_MPIGPU(Vec xin,PetscInt ni,const PetscInt ix[],const PetscScalar y[],InsertMode addv)
{
  PetscErrorCode ierr;
  Vec_SeqGPU *x = (Vec_SeqGPU*)xin->data;
  PetscFunctionBegin;
  if (addv == ADD_VALUES && x->syncState == VEC_GPU) {
    ierr = VecCopyOverD2H(xin,x->cpuptr); CHKERRQ(ierr);
    cudaDeviceSynchronize();
  }
  ierr = VecSetValues_MPI(xin,ni,ix,y,addv);CHKERRQ(ierr);
  x->syncState = VEC_CPU;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecAssemblyEnd_MPIGPU"
PetscErrorCode VecAssemblyEnd_MPIGPU(Vec xin)
{
  PetscErrorCode ierr;
  Vec_SeqGPU *x = (Vec_SeqGPU*)xin->data;

  PetscFunctionBegin;
  ierr = VecAssemblyEnd_MPI(xin);CHKERRQ(ierr);
  x->syncState = VEC_CPU;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecMax_MPIGPU"
PetscErrorCode VecMax_MPIGPU(Vec xin, PetscInt *idx, PetscReal *z)
{
  PetscErrorCode ierr;
  PetscReal work;

  PetscFunctionBegin;
  ierr = VecMax_SeqGPU(xin,idx,&work);CHKERRQ(ierr); 
  if (!idx) {
    ierr = MPI_Allreduce(&work,z,1,MPIU_SCALAR,MPIU_MAX,((PetscObject)xin)->comm);CHKERRQ(ierr);
  } else {
    PetscReal work2[2],z2[2];
    PetscInt  rstart;
    rstart = xin->map->rstart;
    work2[0] = work;
    work2[1] = *idx + rstart;
    ierr = MPI_Allreduce(work2,z2,2,MPIU_REAL,VecMax_Local_Op,((PetscObject)xin)->comm);CHKERRQ(ierr);
    *z   = z2[0];
    *idx = (PetscInt)z2[1];
  }
  PetscFunctionReturn(0);

}


#undef __FUNCT__
#define __FUNCT__ "VecMin_MPIGPU"
PetscErrorCode VecMin_MPIGPU(Vec xin,PetscInt *idx,PetscReal *z)
{
  PetscErrorCode ierr;
  PetscReal      work;

  PetscFunctionBegin;
  /* Find the local Min */
  ierr = VecMin_SeqGPU(xin,idx,&work);CHKERRQ(ierr);

  /* Find the global Min */
  if (!idx) {
    ierr = MPI_Allreduce(&work,z,1,MPIU_REAL,MPIU_MIN,((PetscObject)xin)->comm);CHKERRQ(ierr);
  } else {
    PetscReal work2[2],z2[2];
    PetscInt  rstart;

    ierr = VecGetOwnershipRange(xin,&rstart,PETSC_NULL);CHKERRQ(ierr);
    work2[0] = work;
    work2[1] = *idx + rstart;
    ierr = MPI_Allreduce(work2,z2,2,MPIU_REAL,VecMin_Local_Op,((PetscObject)xin)->comm);CHKERRQ(ierr);
    *z   = z2[0];
    *idx = (PetscInt)z2[1];
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "VecNorm_MPIGPU"
PetscErrorCode VecNorm_MPIGPU(Vec xin,NormType type,PetscReal *z)
{
  PetscReal      sum,work = 0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (type == NORM_2 || type == NORM_FROBENIUS) {
    ierr = VecNorm_SeqGPU(xin,NORM_2,&work);
    work *= work;
    ierr = MPI_Allreduce(&work,&sum,1,MPIU_REAL,MPIU_SUM,((PetscObject)xin)->comm);CHKERRQ(ierr);
    *z = PetscSqrtReal(sum);
    //printf("VecNorm_MPIGPU : z=%1.5g\n",*z);
  } else if (type == NORM_1) {
    /* Find the local part */
    ierr = VecNorm_SeqGPU(xin,NORM_1,&work);CHKERRQ(ierr);
    /* Find the global max */
    ierr = MPI_Allreduce(&work,z,1,MPIU_REAL,MPIU_SUM,((PetscObject)xin)->comm);CHKERRQ(ierr);
  } else if (type == NORM_INFINITY) {
    /* Find the local max */
    ierr = VecNorm_SeqGPU(xin,NORM_INFINITY,&work);CHKERRQ(ierr);
    /* Find the global max */
    ierr = MPI_Allreduce(&work,z,1,MPIU_REAL,MPIU_MAX,((PetscObject)xin)->comm);CHKERRQ(ierr);
  } else if (type == NORM_1_AND_2) {
    PetscReal temp[2];
    ierr = VecNorm_SeqGPU(xin,NORM_1,temp);CHKERRQ(ierr);
    ierr = VecNorm_SeqGPU(xin,NORM_2,temp+1);CHKERRQ(ierr);
    temp[1] = temp[1]*temp[1];
    ierr = MPI_Allreduce(temp,z,2,MPIU_REAL,MPIU_SUM,((PetscObject)xin)->comm);CHKERRQ(ierr);
    z[1] = PetscSqrtReal(z[1]);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecView_MPIGPU_ASCII"
PetscErrorCode VecView_MPIGPU_ASCII(Vec xin,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt          i,work = xin->map->n,cnt,len;
  PetscMPIInt       j,n = 0,size,rank,tag = ((PetscObject)viewer)->tag;
  MPI_Status        status;
  PetscScalar       *values;
  PetscScalar *xarray;
  const char        *name;
  Vec_SeqGPU  *x=(Vec_SeqGPU*)xin->data;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (x->syncState == VEC_GPU) {
    ierr = VecCopyOverD2H(xin,x->cpuptr); CHKERRQ(ierr);
    cudaDeviceSynchronize();
  }
  xarray = x->cpuptr;
  /* determine maximum message to arrive */
  ierr = MPI_Comm_rank(((PetscObject)xin)->comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Reduce(&work,&len,1,MPIU_INT,MPI_MAX,0,((PetscObject)xin)->comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(((PetscObject)xin)->comm,&size);CHKERRQ(ierr);

  if (!rank) {
    ierr = PetscMalloc(len*sizeof(PetscScalar),&values);CHKERRQ(ierr);
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    /*
        MATLAB format and ASCII format are very similar except
        MATLAB uses %18.16e format while ASCII uses %g
    */
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      ierr = PetscObjectGetName((PetscObject)xin,&name);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s = [\n",name);CHKERRQ(ierr);
      for (i=0; i<xin->map->n; i++) {
#if defined(PETSC_USE_COMPLEX)
        if (PetscImaginaryPart(xarray[i]) > 0.0) {
          ierr = PetscViewerASCIIPrintf(viewer,"%18.16e + %18.16ei\n",PetscRealPart(xarray[i]),PetscImaginaryPart(xarray[i]));CHKERRQ(ierr);
        } else if (PetscImaginaryPart(xarray[i]) < 0.0) {
          ierr = PetscViewerASCIIPrintf(viewer,"%18.16e - %18.16ei\n",PetscRealPart(xarray[i]),-PetscImaginaryPart(xarray[i]));CHKERRQ(ierr);
        } else {
          ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",PetscRealPart(xarray[i]));CHKERRQ(ierr);
        }
#else
        ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)xarray[i]);CHKERRQ(ierr);
#endif
      }
      /* receive and print messages */
      for (j=1; j<size; j++) {
        ierr = MPI_Recv(values,(PetscMPIInt)len,MPIU_SCALAR,j,tag,((PetscObject)xin)->comm,&status);CHKERRQ(ierr);
        ierr = MPI_Get_count(&status,MPIU_SCALAR,&n);CHKERRQ(ierr);
        for (i=0; i<n; i++) {
#if defined(PETSC_USE_COMPLEX)
          if (PetscImaginaryPart(values[i]) > 0.0) {
            ierr = PetscViewerASCIIPrintf(viewer,"%18.16e + %18.16e i\n",PetscRealPart(values[i]),PetscImaginaryPart(values[i]));CHKERRQ(ierr);
          } else if (PetscImaginaryPart(values[i]) < 0.0) {
            ierr = PetscViewerASCIIPrintf(viewer,"%18.16e - %18.16e i\n",PetscRealPart(values[i]),-PetscImaginaryPart(values[i]));CHKERRQ(ierr);
          } else {
            ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",PetscRealPart(values[i]));CHKERRQ(ierr);
          }
#else
          ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",values[i]);CHKERRQ(ierr);
#endif
        }
      }
      ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);

    } else {
      ierr = PetscObjectPrintClassNamePrefixType((PetscObject)xin,viewer,"Vector Object");CHKERRQ(ierr);
      if (format != PETSC_VIEWER_ASCII_COMMON) {ierr = PetscViewerASCIIPrintf(viewer,"Process [%d]\n",rank);CHKERRQ(ierr);}
      cnt = 0;
      for (i=0; i<xin->map->n; i++) {
        if (format == PETSC_VIEWER_ASCII_INDEX) {
          ierr = PetscViewerASCIIPrintf(viewer,"%D: ",cnt++);CHKERRQ(ierr);
        }
#if defined(PETSC_USE_COMPLEX)
        if (PetscImaginaryPart(xarray[i]) > 0.0) {
          ierr = PetscViewerASCIIPrintf(viewer,"%g + %g i\n",PetscRealPart(xarray[i]),PetscImaginaryPart(xarray[i]));CHKERRQ(ierr);
        } else if (PetscImaginaryPart(xarray[i]) < 0.0) {
          ierr = PetscViewerASCIIPrintf(viewer,"%g - %g i\n",PetscRealPart(xarray[i]),-PetscImaginaryPart(xarray[i]));CHKERRQ(ierr);
        } else {
          ierr = PetscViewerASCIIPrintf(viewer,"%g\n",PetscRealPart(xarray[i]));CHKERRQ(ierr);
        }
#else
        ierr = PetscViewerASCIIPrintf(viewer,"%g\n",(double)xarray[i]);CHKERRQ(ierr);
#endif
      }
      /* receive and print messages */
      for (j=1; j<size; j++) {
        ierr = MPI_Recv(values,(PetscMPIInt)len,MPIU_SCALAR,j,tag,((PetscObject)xin)->comm,&status);CHKERRQ(ierr);
        ierr = MPI_Get_count(&status,MPIU_SCALAR,&n);CHKERRQ(ierr);
        if (format != PETSC_VIEWER_ASCII_COMMON) {
          ierr = PetscViewerASCIIPrintf(viewer,"Process [%d]\n",j);CHKERRQ(ierr);
        }
        for (i=0; i<n; i++) {
          if (format == PETSC_VIEWER_ASCII_INDEX) {
            ierr = PetscViewerASCIIPrintf(viewer,"%D: ",cnt++);CHKERRQ(ierr);
          }
#if defined(PETSC_USE_COMPLEX)
          if (PetscImaginaryPart(values[i]) > 0.0) {
            ierr = PetscViewerASCIIPrintf(viewer,"%g + %g i\n",PetscRealPart(values[i]),PetscImaginaryPart(values[i]));CHKERRQ(ierr);
          } else if (PetscImaginaryPart(values[i]) < 0.0) {
            ierr = PetscViewerASCIIPrintf(viewer,"%g - %g i\n",PetscRealPart(values[i]),-PetscImaginaryPart(values[i]));CHKERRQ(ierr);
          } else {
            ierr = PetscViewerASCIIPrintf(viewer,"%g\n",PetscRealPart(values[i]));CHKERRQ(ierr);
          }
#else
          ierr = PetscViewerASCIIPrintf(viewer,"%g\n",(double)values[i]);CHKERRQ(ierr);
#endif
        }
      }
    }
    ierr = PetscFree(values);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      /* this may be a collective operation so make sure everyone calls it */
      ierr = PetscObjectGetName((PetscObject)xin,&name);CHKERRQ(ierr);
    }
    /* send values */
    ierr = MPI_Send((void*)xarray,xin->map->n,MPIU_SCALAR,0,tag,((PetscObject)xin)->comm);CHKERRQ(ierr);
  }
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecView_MPIGPU"
PetscErrorCode VecView_MPIGPU(Vec xin,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      iascii,isbinary,isdraw;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);

  if (iascii){
    ierr = VecView_MPIGPU_ASCII(xin,viewer);CHKERRQ(ierr);
    /*
  } else if (isbinary) {
    ierr = VecView_MPIGPU_Binary(xin,viewer);CHKERRQ(ierr);

  } else if (isdraw) {
    PetscViewerFormat format;

    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_DRAW_LG) {
      ierr = VecView_MPIGPU_Draw_LG(xin,viewer);CHKERRQ(ierr);
    } else {
      ierr = VecView_MPIGPU_Draw(xin,viewer);CHKERRQ(ierr);
    }
     */
  } else SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Viewer type %s not supported for this object",((PetscObject)viewer)->type_name);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecDot_MPIGPU"
PetscErrorCode VecDot_MPIGPU(Vec xin,Vec yin,PetscScalar *z)
{
  PetscScalar    sum,work;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecDot_SeqGPU(xin,yin,&work);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&work,&sum,1,MPIU_SCALAR,MPIU_SUM,((PetscObject)xin)->comm);CHKERRQ(ierr);
  *z = sum;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecTDot_MPIGPU"
PetscErrorCode VecTDot_MPIGPU(Vec xin,Vec yin,PetscScalar *z)
{
  PetscScalar    sum,work;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecTDot_SeqGPU(xin,yin,&work);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&work,&sum,1,MPIU_SCALAR,MPIU_SUM,((PetscObject)xin)->comm);CHKERRQ(ierr);
  *z   = sum;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecMDot_MPIGPU"
PetscErrorCode VecMDot_MPIGPU(Vec xin,PetscInt nv,const Vec y[],PetscScalar *z)
{
  PetscScalar    awork[128],*work = awork;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (nv > 128) {
    ierr = PetscMalloc(nv*sizeof(PetscScalar),&work);CHKERRQ(ierr);
  }
  ierr = VecMDot_SeqGPU(xin,nv,y,work);CHKERRQ(ierr);
  ierr = MPI_Allreduce(work,z,nv,MPIU_SCALAR,MPIU_SUM,((PetscObject)xin)->comm);CHKERRQ(ierr);
  if (nv > 128) {
    ierr = PetscFree(work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*MC
   VECMPIGPU - VECMPIGPU = "mpicusp" - The basic parallel vector, modified to use CUSP

   Options Database Keys:
. -vec_type mpicusp - sets the vector type to VECMPIGPU during a call to VecSetFromOptions()

  Level: beginner

.seealso: VecCreate(), VecSetType(), VecSetFromOptions(), VecCreateMpiWithArray(), VECMPI, VecType, VecCreateMPI(), VecCreateMpi()
M*/


#undef __FUNCT__  
#define __FUNCT__ "VecDuplicate_MPIGPU"
PetscErrorCode VecDuplicate_MPIGPU(Vec win,Vec *v)
{
  PetscErrorCode ierr;
  Vec_SeqGPU     *s;
  Vec V; 

  PetscFunctionBegin;
  ierr = VecCreate(((PetscObject)win)->comm,&V);CHKERRQ(ierr);
  ierr = PetscLayoutReference(win->map,&V->map);CHKERRQ(ierr);
  ierr = PetscNewLog(V,Vec_SeqGPU,&s);CHKERRQ(ierr);

  V->data = (void*)s;

  ierr = PetscMemcpy(V->ops,win->ops,sizeof(struct _VecOps));CHKERRQ(ierr);


  /* New vector should inherit stashing property of parent */
  V->stash.donotstash = win->stash.donotstash;
  V->stash.ignorenegidx = win->stash.ignorenegidx;

  ierr = PetscOListDuplicate(((PetscObject)win)->olist,&((PetscObject)V)->olist);CHKERRQ(ierr);
  ierr = PetscFListDuplicate(((PetscObject)win)->qlist,&((PetscObject)V)->qlist);CHKERRQ(ierr);
  V->map->bs    = win->map->bs;
  V->bstash.bs = win->bstash.bs;


  /* Set up local and device storage */
  s->syncState      = VEC_UNALLOC;
  s->unplacedarray=PETSC_NULL;
  s->array_allocated=PETSC_NULL;
  s->array=PETSC_NULL;
  /* create an associated stream */
  cms[0] = cudaStreamCreate(&(s->streamid));
  /* allocate the variable for vector size */
  cms[1]=cudaMalloc((void**)&(s->length),sizeof(int));
  /* send vec length size to device */
  ccs[0]=cudaMemcpyAsync((void*)s->length,
               (void*)&(V->map->n),sizeof(int),cudaMemcpyHostToDevice,s->streamid);
  /* allocate the vector on device */
  cms[2]=cudaMalloc((void**)&(s->devptr),V->map->n*sizeof(double));
  ccs[1]=cudaMemsetAsync((void*)s->devptr,0,V->map->n*sizeof(double),s->streamid);
  /* allocate the variable for vector offsets */
  cms[3]=cudaMalloc((void**)&(s->offset),sizeof(int));
  /* allocate the variable for vector segment length */
  cms[4]=cudaMalloc((void**)&(s->segment),sizeof(int));
  /* allocate the variable for vector single value result */
  cms[5]=cudaMalloc((void**)&(s->zval),sizeof(double));
  cms[6]=cudaMalloc((void**)&(s->scalar),sizeof(double));
  /* using pinned memory (could be a resource hog with very large arrays) */
  ierr = PinnedMalloc(&(s->cpuptr),V->map->n*sizeof(double));CHKERRQ(ierr);


  ierr = VecCheckCUDAStatus(ccs[0],"Copy H2D devlength in VecCreate_SeqGPU");CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(ccs[1],"on device cudaMemSet VecCreate_SeqGPU"); CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(cms[0],"on cudaStreamCreate VecCreate_SeqGPU");  CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(cms[1],"Alloc devlength in VecCreate_SeqGPU");   CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(cms[2],"Alloc of devptr in VecCreate_SeqGPU");   CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(cms[3],"Alloc devoffset in VecCreate_SeqGPU");   CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(cms[4],"Alloc dev segment in VecCreate_SeqGPU"); CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(cms[5],"Alloc dev zval in VecCreate_SeqGPU");    CHKERRQ(ierr);
  ierr = VecCheckCUDAStatus(cms[6],"Alloc dev scalar in VecCreate_SeqGPU");    CHKERRQ(ierr);


  /* ierr = PetscMalloc(V->map->n*sizeof(PetscScalar),&(s->cpuptr)); */
  ierr = PetscMemzero(s->cpuptr,V->map->n*sizeof(double));CHKERRQ(ierr);
  s->syncState=VEC_ALLOC;


  /* change type_name appropriately */
  ierr = PetscObjectChangeTypeName((PetscObject)V,VECMPIGPU);CHKERRQ(ierr);


  ierr = PetscOListDuplicate(((PetscObject)win)->olist,&((PetscObject)V)->olist);CHKERRQ(ierr);
  ierr = PetscFListDuplicate(((PetscObject)win)->qlist,&((PetscObject)V)->qlist);CHKERRQ(ierr);
  V->map->bs    = win->map->bs;
  V->bstash.bs = win->bstash.bs;

  *v = V;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecDotNorm2_MPIGPU"
PetscErrorCode VecDotNorm2_MPIGPU(Vec s,Vec t,PetscScalar *dp,PetscScalar *nm)
{
  PetscErrorCode  ierr;
  PetscScalar     work[2],sum[2];
  PetscFunctionBegin;
  ierr    = VecDotNorm2_SeqGPU(s,t,work,work+1);CHKERRQ(ierr);
  ierr    = MPI_Allreduce(&work,&sum,2,MPIU_SCALAR,MPIU_SUM,((PetscObject)s)->comm);CHKERRQ(ierr);
  *dp     = sum[0];
  *nm     = sum[1];
  //printf("VecDotNorm2_MPIGPU=%1.5g,%1.5g\n",PetscRealPart(*dp),PetscImaginaryPart(*dp));
  //printf("VecDotNorm2_MPIGPU=%1.5g,%1.5g\n",PetscRealPart(*nm),PetscImaginaryPart(*nm));
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "VecCreate_MPIGPU"
PetscErrorCode  VecCreate_MPIGPU(Vec V)
{
  PetscErrorCode ierr;
  Vec_SeqGPU* s = PETSC_NULL;
  PetscFunctionBegin;
  ierr = PetscNewLog(V,Vec_SeqGPU,&s);CHKERRQ(ierr);
  V->data = (void*)s;
  ierr = PetscLayoutSetUp(V->map);CHKERRQ(ierr);

  V->stash.insertmode  = NOT_SET_VALUES;
  /* create the stashes. The block-size for bstash is set later when
     VecSetValuesBlocked is called.
  */
  ierr = VecStashCreate_Private(((PetscObject)V)->comm,1,&V->stash);CHKERRQ(ierr);
  ierr = VecStashCreate_Private(((PetscObject)V)->comm,V->map->bs,&V->bstash);CHKERRQ(ierr);

  
  s->syncState      = VEC_UNALLOC;
  s->unplacedarray=PETSC_NULL;
  s->array_allocated=PETSC_NULL;
  s->array=PETSC_NULL;
  /* create an associated stream */
  cms[0] = cudaStreamCreate(&(s->streamid));
  /* allocate the variable for vector size */
  cms[1]=cudaMalloc((void**)&(s->length),sizeof(int));
  /* send vec length size to device */
  ccs[0]=cudaMemcpyAsync((void*)s->length,
               (void*)&(V->map->n),sizeof(int),cudaMemcpyHostToDevice,s->streamid);
  /* allocate the vector on device */
  cms[2]=cudaMalloc((void**)&(s->devptr),V->map->n*sizeof(double));
  ccs[1]=cudaMemsetAsync((void*)s->devptr,0,V->map->n*sizeof(double),s->streamid);
  /* allocate the variable for vector offsets */
  cms[3]=cudaMalloc((void**)&(s->offset),sizeof(int));
  /* allocate the variable for vector segment length */
  cms[4]=cudaMalloc((void**)&(s->segment),sizeof(int));
  /* allocate the variable for vector single value result */
  cms[5]=cudaMalloc((void**)&(s->zval),sizeof(double));
  cms[6]=cudaMalloc((void**)&(s->scalar),sizeof(double));
  /* using pinned memory (could be a resource hog with very large arrays) */
  ierr = PinnedMalloc(&(s->cpuptr),V->map->n*sizeof(double));CHKERRQ(ierr);

  /* ierr = PetscMalloc(V->map->n*sizeof(PetscScalar),&(s->cpuptr)); */
  ierr = PetscMemzero(s->cpuptr,V->map->n*sizeof(double));CHKERRQ(ierr);
  s->syncState=VEC_ALLOC;

  ierr = PetscObjectChangeTypeName((PetscObject)V,VECMPIGPU);CHKERRQ(ierr);
  V->ops->dotnorm2        = VecDotNorm2_MPIGPU;
  V->ops->waxpy           = VecWAXPY_SeqGPU;
  V->ops->duplicate       = VecDuplicate_MPIGPU;
  V->ops->dot             = VecDot_MPIGPU;
  V->ops->mdot            = VecMDot_MPIGPU;
  V->ops->tdot            = VecTDot_MPIGPU;
  V->ops->norm            = VecNorm_MPIGPU;
  V->ops->view            = VecView_MPIGPU;
  V->ops->max             = VecMax_MPIGPU;
  V->ops->min             = VecMin_MPIGPU;
  V->ops->destroy         = VecDestroy_MPIGPU;

  V->ops->setvalues       = VecSetValues_MPIGPU;
  V->ops->assemblybegin   = VecAssemblyBegin_MPI;
  V->ops->assemblyend     = VecAssemblyEnd_MPIGPU;

  V->ops->getarray        = VecGetArray_SeqGPU;
  V->ops->restorearray    = VecRestoreArray_SeqGPU;
  V->ops->getsize         = VecGetSize_SeqGPU;
  V->ops->duplicatevecs   = VecDuplicateVecs_SeqGPU;
  V->ops->destroyvecs     = VecDestroyVecs_SeqGPU;
  V->ops->scale           = VecScale_SeqGPU;
  V->ops->copy            = VecCopy_SeqGPU;
  V->ops->set             = VecSet_SeqGPU;
  V->ops->swap            = VecSwap_SeqGPU;
  V->ops->axpy            = VecAXPY_SeqGPU;
  V->ops->axpby           = VecAXPBY_SeqGPU;
  V->ops->maxpy           = VecMAXPY_SeqGPU;
  V->ops->aypx            = VecAYPX_SeqGPU;
  V->ops->axpbypcz        = VecAXPBYPCZ_SeqGPU;
  V->ops->pointwisemult   = VecPointwiseMult_SeqGPU;
  V->ops->setrandom       = VecSetRandom_SeqGPU;
  V->ops->replacearray    = VecReplaceArray_SeqGPU;
  V->ops->dot_local       = VecDot_SeqGPU;
  V->ops->tdot_local      = VecTDot_SeqGPU;
  V->ops->norm_local      = VecNorm_SeqGPU;
  V->ops->mdot_local      = VecMDot_SeqGPU;
  V->ops->pointwisedivide = VecPointwiseDivide_SeqGPU;
  /* place array?
     reset array?
     get values?
  */
  ierr = VecSet(V,0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "VecCreate_GPU"
PetscErrorCode  VecCreate_GPU(Vec v)
{
  PetscErrorCode ierr;
  PetscMPIInt    size;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(((PetscObject)v)->comm,&size);CHKERRQ(ierr);
  if (size == 1) {
    ierr = VecSetType(v,VECSEQGPU);CHKERRQ(ierr);
  } else {
    ierr = VecSetType(v,VECMPIGPU);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
EXTERN_C_END





