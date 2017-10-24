# This makefile relies on the variables
# matlabroot and mexext being set. These may be set by running
# sgwt_compile_mex.m within matlab

MEX = ${matlabroot}/bin/mex
INC = -I ${matlabroot}/extern/include
LIB = -L ${matlabroot}/extern/lib
CC = g++
USEOPENMP=yes
#USEOPENMP=no
#DEBUG=yes
DEBUG=no

## using openmp
ifeq ($(USEOPENMP),yes)
	OPENMPFLAGS=-fopenmp
	CXXLIBS=-lgomp
endif

ifeq ($(DEBUG),yes)
	DBFLAGS=-g
	OPTFLAGS=
else
	DBFLAGS=
	OPTFLAGS=-O3
endif

CFLAGS = -Wall $(DBFLAGS) $(OPTFLAGS) $(OPENMPFLAGS) -fPIC
MEXOPTS = CXXLIBS="\$$CXXLIBS $(CXXLIBS)" CXX="g++" CC="gcc" LD="g++" -largeArrayDims

OBJ = cheby_op.o

all : cheby_op_mex cheby_op_adjoint_mex 

cheby_op_mex : ${OBJ} cheby_op_mex.cpp
	${MEX} ${MEXOPTS} cheby_op_mex.cpp ${OBJ}

cheby_op_adjoint_mex : ${OBJ} cheby_op_adjoint_mex.cpp
	${MEX} ${MEXOPTS} cheby_op_adjoint_mex.cpp ${OBJ}

%.o : %.cpp %.h
	${CC} ${INC} -c ${CFLAGS} $<	

clean : 
	rm -f *.o *.${mexext}
