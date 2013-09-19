#requirespackage   'PETSC_HAVE_SUPERLU_MT'

ALL: lib

CFLAGS   = ${SUPERLU_MT_INCLUDE}
FFLAGS   =
SOURCEC  = superlu_mt.c
SOURCEF  =
SOURCEH  =
LIBBASE  = libpetscmat
DIRS     =
MANSEC   = Mat
LOCDIR   = src/mat/impls/aij/seq/superlu_mt/

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test
