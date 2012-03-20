
ALL: lib

CFLAGS   = -mavx -fopenmp -ggdb
FFLAGS   =
SOURCEC  = matblockstructgrid.c matmultbintrinsics_ver1.c 
SOURCEF  = 
SOURCEH  = matblockstructgrid.h
OBJSC    = matblockstructgrid.o matmultbintrinsics_ver1.o
OBJSF    = 
LIBBASE  = libpetscmat
DIRS     = 
MANSEC   = Mat
DIRS     =
LOCDIR   = src/mat/impls/blockstructgrid/seq

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test
