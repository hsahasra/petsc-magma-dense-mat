#requirespackage 'PETSC_HAVE_CUDA'

ALL: lib

CPPFLAGS   =
FFLAGS   =
SOURCEC  =
SOURCECU = vecgpu.cu
SOURCEF  =
SOURCEH  = gpuvecimpl.h
OBJSCU   = vecgpu.o
OBJSF    =
LIBBASE  = libpetscvec
MANSEC   = Vec
LOCDIR   = src/vec/vec/impls/seq/seqgpu/
DIRS     =

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test
