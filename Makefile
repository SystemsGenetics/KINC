MKLROOT = /opt/intel/composerxe-2011.0.084/mkl
CC = gcc -m64
EXE_DIR = /common1/feltus/co-expression_networks/software/bin

CCFLAGS = -g

LDFLAGS =  -g -Wall -lm -lgsl -lgslcblas -fopenmp -llapack -lblas -lpthread -lm

OBJS = ccm.o rmm.o
EXE = ccm rmm

all: ${OBJS}
	${CC} ccm.o  ${LDFLAGS} -o ccm 
	${CC} rmm.o  ${LDFLAGS} -o rmm

ccm.o: ccm.c
	${CC} -c ${CCFLAGS} ccm.c

rmm.o: rmm.c rmm.h 
	${CC} -c ${CCFLAGS} rmm.c


clean:
	rm -f ${OBJS} ${EXE}

install: all
	install -m 0755 ccm ${EXE_DIR}
	install -m 0755 rmm ${EXE_DIR}
	install -m 0755 parse_pearson_bin.pl ${EXE_DIR}
