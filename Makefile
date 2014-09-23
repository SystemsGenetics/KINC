MKLROOT = 
CC = gcc -m64
EXE_DIR = 

CCFLAGS =
INCLUDES = 
LDFLAGS = -Wall -lm -lgsl -lgslcblas -llapack -lblas -lpthread

OBJS = \
  stats/stats.o \
  stats/kurtosis.o \
  stats/sfrancia.o \
  stats/swilk.o \
  stats/royston.o \
  similarity/bspline_mi.o \
  similarity/pearson.o \
  similarity/spearman.o \
  preprocess.o \
  similarity.o \
  threshold.o \
  extract.o \
  kinc.o
EXE = KINC

all: ${OBJS}
	${CC} ${OBJS} ${LDFLAGS} -o kinc

stats/stats.o: stats/stats.c stats/stats.h
	${CC} -c ${CCFLAGS} ${INCLUDES} stats/stats.c -o stats/stats.o

stats/kurtosis.o: stats/kurtosis.c stats/kurtosis.h
	${CC} -c ${CCFLAGS} ${INCLUDES} stats/kurtosis.c -o stats/kurtosis.o

stats/sfrancia.o: stats/sfrancia.c stats/sfrancia.h
	${CC} -c ${CCFLAGS} ${INCLUDES} stats/sfrancia.c -o stats/sfrancia.o

stats/swilk.o: stats/swilk.c stats/swilk.h
	${CC} -c ${CCFLAGS} ${INCLUDES} stats/swilk.c -o stats/swilk.o

stats/royston.o: stats/royston.c stats/royston.h
	${CC} -c ${CCFLAGS} ${INCLUDES} stats/royston.c -o stats/royston.o

similarity/spearman.o: similarity/spearman.c similarity/spearman.h
	${CC} -c ${CCFLAGS} ${INCLUDES} similarity/spearman.c -o similarity/spearman.o

similarity/pearson.o: similarity/pearson.c similarity/pearson.h
	${CC} -c ${CCFLAGS} ${INCLUDES} similarity/pearson.c -o similarity/pearson.o

similarity/bspline_mi.o: similarity/bspline_mi.c similarity/bspline_mi.h
	${CC} -c ${CCFLAGS} ${INCLUDES} similarity/bspline_mi.c -o similarity/bspline_mi.o

preprocess.o: preprocess.c preprocess.h
	${CC} -c ${CCFLAGS} ${INCLUDES} preprocess.c
	
similarity.o: similarity.c similarity.h
	${CC} -c ${CCFLAGS} ${INCLUDES} similarity.c

threshold.o: threshold.c threshold.h
	${CC} -c ${CCFLAGS} ${INCLUDES} threshold.c

extract.o: extract.c extract.h
	${CC} -c ${CCFLAGS} ${INCLUDES} extract.c

kink.o: kink.c kink.h
	${CC} -c ${CCFLAGS} ${INCLUDES} kink.c

clean:
	rm -f ${OBJS} ${EXE}

install: all
	install -m 0755 kink ${EXE_DIR}
