#requirespackage  'PETSC_HAVE_SSL'

ALL: lib

CFLAGS    =
FFLAGS    =
CPPFLAGS  =
SOURCEC	  = client.c
SOURCEF	  =
SOURCEH	  =
MANSEC	  = Sys
LIBBASE	  = libpetscsys
DIRS	  =
LOCDIR	  = src/sys/webclient/

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test

petscwebclient2: petscwebclient2.o  chkopts
	${CLINKER} -o petscwebclient2 petscwebclient2.o ${PETSC_SYS_LIB} -lssl -lcrypto
	${RM} -f petscwebclient2.o

googledriveauthorize: .o  chkopts
	${CLINKER} -o googledriveauthorize googledriveauthorize.o ${PETSC_SYS_LIB} -lssl -lcrypto
	${RM} -f googledriveauthorize.o


