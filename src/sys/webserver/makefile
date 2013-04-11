#requirespackage  'PETSC_HAVE_SERVER'
#
SPECIALLIB     = yes
SPECIALFASTLIB = yes
speciallib:
	-@echo "BEGINNING TO COMPILE PetscWebServer "
	-@if [ ! -d "${PETSC_DIR}/${PETSC_ARCH}/bin" ] ; then ${MKDIR}  ${PETSC_DIR}/${PETSC_ARCH}/bin; fi
	-@${OMAKE}  petscwebserver PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR}
	-@${OMAKE}  javascript PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR}
specialfastlib: speciallib

ALL: lib

CFLAGS    =
FFLAGS    =
CPPFLAGS  =
SOURCEC	  =
SOURCEF	  =
SOURCEH	  =
MANSEC	  = Sys
LIBBASE	  = libpetscsys
DIRS	  =
LOCDIR	  = src/sys/objects/webserver

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test

petscwebserver: petscwebserver.o webserve.o  chkopts
	@-lib=`echo  ${AMS_LIB} | sed "s?amspub?amsacc?g"`; \
        ${CLINKER} -o petscwebserver petscwebserver.o webserve.o $$lib ${PETSC_SYS_LIB}
	${RM} -f petscwebserver.o webserve.o
	${MV} -f petscwebserver ${PETSC_DIR}/${PETSC_ARCH}/bin

petscwebclient: petscwebclient.o  chkopts
	${CLINKER} -o petscwebclient petscwebclient.o ${PETSC_SYS_LIB}
	${RM} -f petscwebclient.o

javascript:
	cd ${PETSC_DIR}/include/web; \
           ${PYTHON} ${PYJSBUILD} --frame=AMSSnoopObjects.html --enable-debug --enable-wrap-calls --output . AMSSnoopObjects.py ;\
           ${PYTHON} ${PYJSBUILD} --frame=AMSSnoopStack.html --enable-debug --enable-wrap-calls --output . AMSSnoopStack.py
