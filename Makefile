# build parameters
DEBUG  ?= 0
SINGLE_PRECISION_MIXMOD ?= 1
CXX    ?= g++
LINK.o = $(LINK.cc)  # use the C++ linker

# compiler flags
GSLDIR    ?= /software/gsl/2.3
MIXMODDIR ?= $(HOME)/software/libmixmod

CPPFLAGS = -I $(GSLDIR)/include \
           -I $(MIXMODDIR)/include
CXXFLAGS = -std=c++11

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -Wall -fno-inline
else
CXXFLAGS += -O3 -Wno-unused-result
endif

ifeq ($(SINGLE_PRECISION_MIXMOD), 1)
CPPFLAGS += -DSINGLE_PRECISION_MIXMOD=1
endif

# linker flags and libraries
LDFLAGS = -L $(GSLDIR)/lib \
          -L $(MIXMODDIR)/lib
LDLIBS = -lm \
         -lgsl -lgslcblas \
         -llapack -lblas -lpthread \
         -lmixmod -lmixmod_newmat

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

$(BIN): $(OBJS)

# a C++ source file depends on its matching header file
%.o: %.cpp %.h
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(BIN)
