MKLROOT = 
#CC = mpic++ -m64
CC = g++ -m64
EXE_DIR = 

#MPI_INCLUDES =  $(shell mpic++ --showme:compile)
#MPI_LDLINK = $(shell mpic++ --showme:link)
MPI_INCLUDES =  
MPI_LDLINK = 

# Debugging CFLAGS
CFLAGS = -g -Wall -fno-inline
# Non-debugging CFLAGS
#CFLAGS = -Wall 
INCLUDES = -I/usr/local/include
LDFLAGS = -Wall -O3 -lm -lgsl -lgslcblas -llapack -lblas -lpthread -lmixmod -lmixmod_newmat -g

OBJS = \
  general/error.o \
  general/misc.o \
  general/vector.o \
  ematrix/EMatrix.o \
  similarity/PairWiseSet.o \
  similarity/PairWiseSimilarity.o \
  similarity/MISimilarity.o \
  similarity/PearsonSimilarity.o \
  similarity/SpearmanSimilarity.o \
  similarity/clustering/PairWiseCluster.o \
  similarity/clustering/PairWiseClusterList.o \
  similarity/clustering/PairWiseClusterWriter.o \
  similarity/clustering/PairWiseClustering.o \
  similarity/clustering/MixtureModelPWClusters.o \
  similarity/clustering/MixtureModelClustering.o \
  similarity/Similarity.o \
  threshold/ThresholdMethod.o \
  threshold/RMTThreshold.o \
  extract/SimilarityMatrix.o \
  extract/SimMatrixBinary.o \
  extract/SimMatrixTabCluster.o \
  kinc.o
EXE = kinc

all: ${OBJS}
	${CC} ${OBJS} ${LDFLAGS} ${MPI_LDLINK} -o ${EXE}

general/misc.o: general/misc.cpp general/misc.h
	${CC} -c ${CFLAGS} ${INCLUDES} general/misc.cpp -o general/misc.o

general/vector.o: general/vector.cpp general/vector.h
	${CC} -c ${CFLAGS} ${INCLUDES} general/vector.cpp -o general/vector.o

general/error.o: general/error.cpp general/error.h
	${CC} -c ${CFLAGS} ${INCLUDES} general/error.cpp -o general/error.o
	
ematrix/EMatrix.o: ematrix/EMatrix.cpp ematrix/EMatrix.h
	${CC} -c ${CFLAGS} ${INCLUDES} ematrix/EMatrix.cpp -o ematrix/EMatrix.o

#stats/stats.o: stats/stats.cpp stats/stats.h
#	${CC} -c ${CFLAGS} ${INCLUDES} stats/stats.cpp -o stats/stats.o
#
#stats/kurtosis.o: stats/kurtosis.cpp stats/kurtosis.h
#	${CC} -c ${CFLAGS} ${INCLUDES} stats/kurtosis.cpp -o stats/kurtosis.o
#
#stats/sfrancia.o: stats/sfrancia.cpp stats/sfrancia.h
#	${CC} -c ${CFLAGS} ${INCLUDES} stats/sfrancia.cpp -o stats/sfrancia.o
#
#stats/swilk.o: stats/swilk.cpp stats/swilk.h
#	${CC} -c ${CFLAGS} ${INCLUDES} stats/swilk.cpp -o stats/swilk.o
#
#stats/royston.o: stats/royston.cpp stats/royston.h
#	${CC} -c ${CFLAGS} ${INCLUDES} stats/royston.cpp -o stats/royston.o
#
#stats/outlier.o: stats/outlier.cpp stats/outlier.h
#	${CC} -c ${CFLAGS} ${INCLUDES} stats/outlier.cpp -o stats/outlier.o

similarity/PairWiseSet.o: similarity/PairWiseSet.cpp similarity/PairWiseSet.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/PairWiseSet.cpp -o similarity/PairWiseSet.o

similarity/PairWiseSimilarity.o: similarity/PairWiseSimilarity.cpp similarity/PairWiseSimilarity.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/PairWiseSimilarity.cpp -o similarity/PairWiseSimilarity.o

similarity/SpearmanSimilarity.o: similarity/SpearmanSimilarity.cpp similarity/SpearmanSimilarity.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/SpearmanSimilarity.cpp -o similarity/SpearmanSimilarity.o

similarity/PearsonSimilarity.o: similarity/PearsonSimilarity.cpp similarity/PearsonSimilarity.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/PearsonSimilarity.cpp -o similarity/PearsonSimilarity.o

similarity/MISimilarity.o: similarity/MISimilarity.cpp similarity/MISimilarity.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/MISimilarity.cpp -o similarity/MISimilarity.o

#clustering/meanshift.o: clustering/meanshift.cpp clustering/meanshift.h
#	${CC} -c ${CFLAGS} ${INCLUDES} clustering/meanshift.cpp -o clustering/meanshift.o

similarity/clustering/PairWiseCluster.o: similarity/clustering/PairWiseCluster.cpp similarity/clustering/PairWiseCluster.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/clustering/PairWiseCluster.cpp -o similarity/clustering/PairWiseCluster.o

similarity/clustering/PairWiseClusterList.o: similarity/clustering/PairWiseCluster.cpp similarity/clustering/PairWiseClusterList.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/clustering/PairWiseClusterList.cpp -o similarity/clustering/PairWiseClusterList.o

similarity/clustering/PairWiseClusterWriter.o: similarity/clustering/PairWiseCluster.cpp similarity/clustering/PairWiseClusterWriter.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/clustering/PairWiseClusterWriter.cpp -o similarity/clustering/PairWiseClusterWriter.o

similarity/clustering/PairWiseClustering.o: similarity/clustering/PairWiseClustering.cpp similarity/clustering/PairWiseClustering.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/clustering/PairWiseClustering.cpp -o similarity/clustering/PairWiseClustering.o

similarity/clustering/MixtureModelPWClusters.o: similarity/clustering/MixtureModelPWClusters.cpp similarity/clustering/MixtureModelPWClusters.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/clustering/MixtureModelPWClusters.cpp -o similarity/clustering/MixtureModelPWClusters.o

similarity/clustering/MixtureModelClustering.o: similarity/clustering/MixtureModelClustering.cpp similarity/clustering/MixtureModelClustering.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/clustering/MixtureModelClustering.cpp -o similarity/clustering/MixtureModelClustering.o

similarity/Similarity.o: similarity/Similarity.cpp similarity/Similarity.h
	${CC} -c ${CFLAGS} ${INCLUDES} similarity/Similarity.cpp -o similarity/Similarity.o

threshold/ThresholdMethod.o: threshold/ThresholdMethod.cpp threshold/ThresholdMethod.h
	${CC} -c ${CFLAGS} ${INCLUDES} threshold/ThresholdMethod.cpp -o threshold/ThresholdMethod.o

threshold/RMTThreshold.o: threshold/RMTThreshold.cpp threshold/RMTThreshold.h
	${CC} -c ${CFLAGS} ${INCLUDES} threshold/RMTThreshold.cpp -o threshold/RMTThreshold.o

extract/SimilarityMatrix.o: extract/SimilarityMatrix.cpp extract/SimilarityMatrix.h
	${CC} -c ${CFLAGS} ${INCLUDES} extract/SimilarityMatrix.cpp -o extract/SimilarityMatrix.o

extract/SimMatrixBinary.o: extract/SimMatrixBinary.cpp extract/SimMatrixBinary.h
	${CC} -c ${CFLAGS} ${INCLUDES} extract/SimMatrixBinary.cpp -o extract/SimMatrixBinary.o

extract/SimMatrixTabCluster.o: extract/SimMatrixTabCluster.cpp extract/SimMatrixTabCluster.h
	${CC} -c ${CFLAGS} ${INCLUDES} extract/SimMatrixTabCluster.cpp -o extract/SimMatrixTabCluster.o

kinc.o: kinc.cpp kinc.h
	${CC} -c ${CFLAGS} ${INCLUDES} ${MPI_INCLUDES} kinc.cpp -o kinc.o

clean:
	rm -f ${OBJS} ${EXE}

install: all
	install -m 0755 kink ${EXE_DIR}
