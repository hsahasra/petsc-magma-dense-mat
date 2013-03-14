#requirespackage 'PETSC_HAVE_CUDA'
ALL: lib

CFLAGS   = 
FFLAGS   =
SOURCECU = mpigpu.cu
SOURCEF  =
SOURCEH  = #mpigpu.h
LIBBASE  = libpetscvec
MANSEC   = Vec
LOCDIR   = src/vec/vec/impls/mpi/mpigpu/

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test
