


#include <../src/mat/impls/aij/seq/aij.h>          /*I "petscmat.h" I*/
PETSC_CUDA_EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "MatMult_SeqAIJ"
PetscErrorCode MatMult_SeqAIJ(Mat A,Vec xx,Vec yy){
  PetscFunctionBegin;
  cudaDeviceSynchronize();
  Mat_SeqAIJ         *a = (Mat_SeqAIJ*)A->data;
  PetscErrorCode     ierr;
  //PetscInt           matsize=A->rmap->n;
  PetscInt           i,n,m,nnz,*rowoffsets,*cindices;
  PetscScalar        *aa;
  cusparseHandle_t   handle=0;
  cusparseMatDescr_t descrip=0;
  cudaError_t        cs;
  cusparseStatus_t   csrs=CUSPARSE_STATUS_SUCCESS;
  
  /* set up cusparse library handle and environment */
  if(cusparseCreate(&handle)!=CUSPARSE_STATUS_SUCCESS){
     printf("cusparse handle creation error.\nExiting...\n");
     PetscFunctionReturn(PETSC_ERR_LIB);
  }

  if(cusparseCreateMatDescr(&descrip)!=CUSPARSE_STATUS_SUCCESS){
     printf("cusparse matrix descriptor creation error.\nExiting...\n");
     PetscFunctionReturn(PETSC_ERR_LIB);
  }
  cusparseSetMatType(descrip,CUSPARSE_MATRIX_TYPE_GENERAL);/* default anyways... */
  cusparseSetMatIndexBase(descrip,CUSPARSE_INDEX_BASE_ZERO);

  aa  = a->a;                     /* nonzero elements */
  rowoffsets  = a->i;             /* pointer to beginning of each row */
  cindices = a->j;                /* column indices */
  nnz = a->nz;                    /* nonzeros */

  /* declare and allocate device csr memory and dense vectors x, y */
  ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
  int* dev_csrRowOffsets;
  int* dev_csrIndices;
  double* dev_dataA;

  /* Allocate CSR device memory */

  cs=cudaMalloc((void**)&dev_csrRowOffsets,(m+1)*sizeof(int));
  if(cs!=cudaSuccess)printf("Error1: %s\n",cudaGetErrorString(cs));

  cs=cudaMalloc((void**)&dev_csrIndices,nnz*sizeof(int));
  if(cs!=cudaSuccess)printf("Error2: %s\n",cudaGetErrorString(cs));

  cs=cudaMalloc((void**)&dev_dataA,nnz*sizeof(double));
  if(cs!=cudaSuccess)printf("Error3: %s\n",cudaGetErrorString(cs));

  /* Send off data to device */

  cs=cudaMemcpy(dev_csrRowOffsets,rowoffsets,(m+1)*sizeof(int),cudaMemcpyHostToDevice);
  if(cs!=cudaSuccess)printf("Error4: %s\n",cudaGetErrorString(cs));

  cs=cudaMemcpy(dev_csrIndices,cindices,nnz*sizeof(int),cudaMemcpyHostToDevice);
  if(cs!=cudaSuccess)printf("Error5: %s\n",cudaGetErrorString(cs));

  cs=cudaMemcpy(dev_dataA,aa,nnz*sizeof(double),cudaMemcpyHostToDevice);
  if(cs!=cudaSuccess)printf("Error6: %s\n",cudaGetErrorString(cs));

  Vec_SeqGPU *xd=(Vec_SeqGPU*)xx->data;
  Vec_SeqGPU *yd=(Vec_SeqGPU*)yy->data;

  /*if(yd->syncState == VEC_GPU || yd->syncState == VEC_SYNCHED){
    cs=cudaMemcpy(yd->cpuptr,yd->devptr,yy->map->n*sizeof(double),cudaMemcpyDeviceToHost);
    if(cs!=cudaSuccess)printf("Error7: %s\n",cudaGetErrorString(cs));
    }else */
  if(yd->syncState == VEC_CPU){
    cs=cudaMemcpy(yd->devptr,yd->cpuptr,yy->map->n*sizeof(double),cudaMemcpyHostToDevice);
    if(cs!=cudaSuccess)printf("Error8: %s\n",cudaGetErrorString(cs));
  }
  //for(i=0;i<yy->map->n;i++){
  //   if(yd->cpuptr[i]!=0.)printf("preMM Y[%d]: %e\n",i,yd->cpuptr[i]);
  //}
  //ierr = VecCheck_SeqGPU(yy);CHKERRQ(ierr);

  /*if(xd->syncState == VEC_GPU || xd->syncState == VEC_SYNCHED){
    cs=cudaMemcpy(xd->cpuptr,xd->devptr,xx->map->n*sizeof(double),cudaMemcpyDeviceToHost);
    if(cs!=cudaSuccess)printf("Error9: %s\n",cudaGetErrorString(cs));
    }else*/
   if(xd->syncState == VEC_CPU){
    cs=cudaMemcpy(xd->devptr,xd->cpuptr,xx->map->n*sizeof(double),cudaMemcpyHostToDevice);
    if(cs!=cudaSuccess)printf("Error10: %s\n",cudaGetErrorString(cs));
  }
  // for(i=0;i<xx->map->n;i++){
  //    if(xd->cpuptr[i]!=0.)printf("preMM X[%d]: %e\n",i,xd->cpuptr[i]);
  //}
  cudaDeviceSynchronize();
  csrs=cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                      m,n,1.0,descrip,
		      dev_dataA,dev_csrRowOffsets,dev_csrIndices,
		      xd->devptr,0.,yd->devptr);

  //printf("Error code %d returned from cusparseDcsrmv call: ",csrs);
  cudaDeviceSynchronize();
  if(csrs!=CUSPARSE_STATUS_SUCCESS) {
    printf("SpMV cuspare lib failed.\n");
    PetscFunctionReturn(PETSC_ERR_LIB);
  }
  yd->syncState = VEC_GPU;

  if(dev_csrRowOffsets)cudaFree(dev_csrRowOffsets);
  if(dev_csrIndices)cudaFree(dev_csrIndices);
  if(dev_dataA)cudaFree(dev_dataA);
  //cudaDeviceSynchronize();

  //cs=cudaMemcpy(yd->cpuptr,yd->devptr,yy->map->n*sizeof(double),cudaMemcpyDeviceToHost);
  //if(cs!=cudaSuccess)printf("Error11: %s\n",cudaGetErrorString(cs));

  //for(i=0;i<yy->map->n;i++){
  //   if(yd->cpuptr[i]!=0.)printf("postMM Y[%d]: %e\n",i,yd->cpuptr[i]);
  //}

  PetscFunctionReturn(0);
}
PETSC_CUDA_EXTERN_C_END


