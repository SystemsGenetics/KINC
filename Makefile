# build parameters
DEBUG ?= 0

# compiler, compiler flags
CXX = g++
CXXFLAGS = -std=c++11
LDFLAGS = -lm -lgsl -lgslcblas -llapack -lblas -lpthread -lmixmod -lmixmod_newmat

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -pg -Wall
else
CXXFLAGS += -O3 -Wno-unused-result
endif

# output files
OBJS = \
	general/error.o \
	general/misc.o \
	general/vector.o \
	stats/stats.o \
	stats/outlier.o \
	stats/swilk.o \
	stats/kurtosis.o \
	stats/sfrancia.o \
	stats/royston.o \
	ematrix/EMatrix.o \
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
	indexer/Indexer.o \
	indexer/IndexQuery.o \
	indexer/RunIndex.o \
	indexer/RunQuery.o \
	threshold/methods/ThresholdMethod.o \
	threshold/methods/RMTThreshold.o \
	threshold/RunThreshold.o \
	extract/SimilarityMatrix.o \
	extract/SimMatrixBinary.o \
	extract/SimMatrixTabCluster.o \
	extract/RunExtract.o \
	kinc.o
BINS = kinc

all: echo $(BINS)

echo:
	$(info DEBUG     = $(DEBUG))
	$(info CXX       = $(CXX))
	$(info CXXFLAGS  = $(CXXFLAGS))
	$(info LDFLAGS   = $(LDFLAGS))

kinc: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.c %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(BINS)
