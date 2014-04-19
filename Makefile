MKLROOT = /opt/intel/composerxe-2011.0.084/mkl
CC = gcc -m64
EXE_DIR = /common1/feltus/co-expression_networks/software/bin

CCFLAGS = -g
LDFLAGS = -lm  -lgsl -lgslcblas -g -Wall -fopenmp -lpthread -fopenmp
#LDFLAGS = -lm -lgsl -lgslcblas -g -Wall -fopenmp -lpthread -fopenmp  \
#	-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core 

OBJS = ccm.o
#OBJS = ccm.o RandomMatrixModeling.v3.1.o RandomMatrix.v2-3.o
EXE = ccm rmm

all: ${OBJS}
	${CC} ccm.o  ${LDFLAGS} -o ccm 
#	${CC} RandomMatrixModeling.v3.1.o RandomMatrix.v2-3.o  ${LDFLAGS} -o rmm

ccm.o: ccm.c
	${CC} -c ${CCFLAGS} ccm.c

#RandomMatrixModeling.v3.1.o: RandomMatrixModeling.v3.1.c RandomMatrix.h 
#	${CC} -c ${CCFLAGS} RandomMatrixModeling.v3.1.c

#RandomMatrix.v2-3.o: RandomMatrix.v2-3.c
#	${CC} -c ${CCFLAGS} RandomMatrix.v2-3.c

clean:
	rm -f ${OBJS} ${EXE}

install: all
	install -m 0755 ccm ${EXE_DIR}
	install -m 0755 rmm ${EXE_DIR}
	install -m 0755 parse_pearson_bin.pl ${EXE_DIR}
