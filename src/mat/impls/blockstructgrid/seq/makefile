
ALL: lib

CFLAGS   = -xavx -openmp
FFLAGS   =
SOURCEC  = matblockstructgrid.c matmultbintrinsics_sse.c matmultbintrinsics_avx.c 
SOURCEF  = 
SOURCEH  = matblockstructgrid.h
OBJSC    = matblockstructgrid.o matmultbintrinsics_sse.o matmultbintrinsics_avx.o
OBJSF    = 
LIBBASE  = libpetscmat
DIRS     = 
MANSEC   = Mat
DIRS     =
LOCDIR   = src/mat/impls/blockstructgrid/seq

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test
