MKLROOT = 
CC = gcc -m64
EXE_DIR = 

CCFLAGS =
INCLUDES = 
LDFLAGS = -Wall -lm -lgsl -lgslcblas -llapack -lblas -lpthread

OBJS = similarity/bspline_mi.o similarity/pearson.o similarity.o threshold.o RMTGeneNet.o
EXE = RMTGeneNet

all: ${OBJS}
	${CC} ${OBJS} ${LDFLAGS} -o RMTGeneNet

similarity/pearson.o: similarity/pearson.c similarity/pearson.h
	${CC} -c ${CCFLAGS} ${INCLUDES} similarity/pearson.c -o similarity/pearson.o

similarity/bspline_mi.o: similarity/bspline_mi.c similarity/bspline_mi.h
	${CC} -c ${CCFLAGS} ${INCLUDES} similarity/bspline_mi.c -o similarity/bspline_mi.o

similarity.o: similarity.c similarity.h
	${CC} -c ${CCFLAGS} ${INCLUDES} similarity.c
	
threshold.o: threshold.c threshold.h
	${CC} -c ${CCFLAGS} ${INCLUDES} threshold.c

RMTGeneNet.o: RMTGeneNet.c RMTGeneNet.h
	${CC} -c ${CCFLAGS} ${INCLUDES} RMTGeneNet.c

clean:
	rm -f ${OBJS} ${EXE}

install: all
	install -m 0755 RMTGeneNet ${EXE_DIR}
