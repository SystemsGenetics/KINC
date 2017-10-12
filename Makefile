# build parameters
DEBUG  ?= 0
CXX    ?= g++

# compiler flags
GSLDIR    ?= /software/gsl/2.3
MIXMODDIR ?= $(HOME)/software/libmixmod

CXXFLAGS = -I $(GSLDIR)/include \
           -I $(MIXMODDIR)/include \
           -std=c++11

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -Wall -fno-inline
else
CXXFLAGS += -O3 -Wno-unused-result
endif

# linker flags
LDFLAGS = -lm \
          -L $(GSLDIR) -lgsl -lgslcblas \
          -llapack -lblas -lpthread \
          -L $(MIXMODDIR)/lib -lmixmod -lmixmod_newmat

# object files, executables
OBJS = \
  ematrix/EMatrix.o \
  extract/SimilarityMatrix.o \
  extract/SimMatrixBinary.o \
  extract/SimMatrixTabCluster.o \
  extract/RunExtract.o \
  general/error.o \
  general/misc.o \
  general/vector.o \
  indexer/Indexer.o \
  indexer/IndexQuery.o \
  indexer/RunIndex.o \
  indexer/RunQuery.o \
  similarity/PairWiseSet.o \
  similarity/methods/PairWiseSimilarity.o \
  similarity/methods/MISimilarity.o \
  similarity/methods/PearsonSimilarity.o \
  similarity/methods/SpearmanSimilarity.o \
  similarity/clustering/PairWiseCluster.o \
  similarity/clustering/PairWiseClusterList.o \
  similarity/clustering/PairWiseClusterWriter.o \
  similarity/clustering/PairWiseClustering.o \
  similarity/clustering/MixtureModelPWClusters.o \
  similarity/clustering/MixtureModelClustering.o \
  similarity/RunSimilarity.o \
  stats/stats.o \
  stats/outlier.o \
  stats/swilk.o \
  stats/kurtosis.o \
  stats/sfrancia.o \
  stats/royston.o \
  threshold/methods/ThresholdMethod.o \
  threshold/methods/RMTThreshold.o \
  threshold/RunThreshold.o \
  kinc.o
BIN = kinc

all: $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $(BIN)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(BIN)
