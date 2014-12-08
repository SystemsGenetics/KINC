MKLROOT = 
CC = mpicc -m64
EXE_DIR = 

MPI_INCLUDES =  $(shell mpicc --showme:compile)
MPI_LDLINK = $(shell mpicc --showme:link)

CFLAGS = -g -Wall
INCLUDES = 
LDFLAGS = -Wall -lm -lgsl -lgslcblas -llapack -lblas -lpthread 

OBJS = \
  error.o \
  misc.o \
  vector.o \
  stats/stats.o \
  stats/kurtosis.o \
  stats/sfrancia.o \
  stats/swilk.o \
  stats/royston.o \
  stats/outlier.o \
  stats/meanshift.o \
  similarity/bspline_mi.o \
  similarity/pearson.o \
  similarity/spearman.o \
  dimreduce.o \
  similarity.o \
  threshold.o \
  extract.o \
  kinc.o
EXE = kinc

all: ${OBJS}
	${CC} ${OBJS} ${LDFLAGS} ${MPI_LDLINK} -o ${EXE}

misc.o: misc.c misc.h
	${CC} -c ${CFLAGS} ${INCLUDES} misc.c

vector.o: vector.c vector.h
	${CC} -c ${CFLAGS} ${INCLUDES} vector.c

error.o: error.c error.h
	${CC} -c ${CFLAGS} ${INCLUDES} error.c

stats/stats.o: stats/stats.c stats/stats.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/stats.c -o stats/stats.o

stats/kurtosis.o: stats/kurtosis.c stats/kurtosis.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/kurtosis.c -o stats/kurtosis.o

stats/sfrancia.o: stats/sfrancia.c stats/sfrancia.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/sfrancia.c -o stats/sfrancia.o

stats/swilk.o: stats/swilk.c stats/swilk.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/swilk.c -o stats/swilk.o

stats/royston.o: stats/royston.c stats/royston.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/royston.c -o stats/royston.o

stats/outlier.o: stats/outlier.c stats/outlier.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/outlier.c -o stats/outlier.o

stats/meanshift.o: stats/meanshift.c stats/meanshift.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/meanshift.c -o stats/meanshift.o

similarity/spearman.o: similarity/spearman.c similarity/spearman.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/spearman.c -o similarity/spearman.o

similarity/pearson.o: similarity/pearson.c similarity/pearson.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/pearson.c -o similarity/pearson.o

similarity/bspline_mi.o: similarity/bspline_mi.c similarity/bspline_mi.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/bspline_mi.c -o similarity/bspline_mi.o

dimreduce.o: dimreduce.c dimreduce.h
	${CC} -c ${CFLAGS} ${INCLUDES} ${MPI_INCLUDES} dimreduce.c

similarity.o: similarity.c similarity.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity.c

threshold.o: threshold.c threshold.h
	${CC} -c ${CFLAGS} ${INCLUDES} threshold.c

extract.o: extract.c extract.h
	${CC} -c ${CFLAGS} ${INCLUDES} extract.c

kinc.o: kinc.c kinc.h
	${CC} -c ${CFLAGS} ${INCLUDES} ${MPI_INCLUDES} kinc.c

clean:
	rm -f ${OBJS} ${EXE}

install: all
	install -m 0755 kink ${EXE_DIR}
