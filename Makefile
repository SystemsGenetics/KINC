MKLROOT = 
CC = mpic++ -m64
EXE_DIR = 

MPI_INCLUDES =  $(shell mpic++ --showme:compile)
MPI_LDLINK = $(shell mpic++ --showme:link)

# Debugging CFLAGS
CFLAGS = -g -Wall -fno-inline
# Non-debugging CFLAGS
#CFLAGS = -Wall 
INCLUDES = -I/usr/local/include
LDFLAGS = -Wall -O3 -lm -lgsl -lgslcblas -llapack -lblas -lpthread -lmixmod -lmixmod_newmat -g

OBJS = \
  error.o \
  misc.o \
  vector.o \
  ematrix.o \
  stats/stats.o \
  stats/outlier.o \
  stats/meanshift.o \
  similarity/bspline_mi.o \
  similarity/pearson.o \
  similarity/spearman.o \
  clustering/clusters.o \
  clustering/mixmod.o \
  dimreduce.o \
  similarity.o \
  threshold.o \
  extract.o \
  kinc.o
EXE = kinc

all: ${OBJS}
	${CC} ${OBJS} ${LDFLAGS} ${MPI_LDLINK} -o ${EXE}

misc.o: misc.cpp misc.h
	${CC} -c ${CFLAGS} ${INCLUDES} misc.cpp

vector.o: vector.cpp vector.h
	${CC} -c ${CFLAGS} ${INCLUDES} vector.cpp

error.o: error.cpp error.h
	${CC} -c ${CFLAGS} ${INCLUDES} error.cpp
	
ematrix.o: ematrix.cpp ematrix.h
	${CC} -c ${CFLAGS} ${INCLUDES} ematrix.cpp

stats/stats.o: stats/stats.cpp stats/stats.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/stats.cpp -o stats/stats.o

stats/kurtosis.o: stats/kurtosis.cpp stats/kurtosis.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/kurtosis.cpp -o stats/kurtosis.o

stats/sfrancia.o: stats/sfrancia.cpp stats/sfrancia.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/sfrancia.cpp -o stats/sfrancia.o

stats/swilk.o: stats/swilk.cpp stats/swilk.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/swilk.cpp -o stats/swilk.o

stats/royston.o: stats/royston.cpp stats/royston.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/royston.cpp -o stats/royston.o

stats/outlier.o: stats/outlier.cpp stats/outlier.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/outlier.cpp -o stats/outlier.o

stats/meanshift.o: stats/meanshift.cpp stats/meanshift.h
	${CC} -c ${CFLAGS} ${INCLUDES} stats/meanshift.cpp -o stats/meanshift.o

similarity/spearman.o: similarity/spearman.cpp similarity/spearman.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/spearman.cpp -o similarity/spearman.o

similarity/pearson.o: similarity/pearson.cpp similarity/pearson.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/pearson.cpp -o similarity/pearson.o

similarity/bspline_mi.o: similarity/bspline_mi.cpp similarity/bspline_mi.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/bspline_mi.cpp -o similarity/bspline_mi.o

clustering/clusters.o: clustering/clusters.cpp clustering/clusters.h
	${CC} -c ${CFLAGS} ${INCLUDES} clustering/clusters.cpp -o clustering/clusters.o

clustering/mixmod.o: clustering/mixmod.cpp clustering/mixmod.h
	${CC} -c ${CFLAGS} ${INCLUDES} clustering/mixmod.cpp -o clustering/mixmod.o

dimreduce.o: dimreduce.cpp dimreduce.h
	${CC} -c ${CFLAGS} ${INCLUDES} ${MPI_INCLUDES} dimreduce.cpp

similarity.o: similarity.cpp similarity.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity.cpp

threshold.o: threshold.cpp threshold.h
	${CC} -c ${CFLAGS} ${INCLUDES} threshold.cpp

extract.o: extract.cpp extract.h
	${CC} -c ${CFLAGS} ${INCLUDES} extract.cpp

kinc.o: kinc.cpp kinc.h
	${CC} -c ${CFLAGS} ${INCLUDES} ${MPI_INCLUDES} kinc.cpp

clean:
	rm -f ${OBJS} ${EXE}

install: all
	install -m 0755 kink ${EXE_DIR}
