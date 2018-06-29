# Basic Settings
TARGET = kinccore
TEMPLATE = lib
CONFIG += staticlib

# Build settings
DESTDIR = $$PWD/../../build/libs/

# Qt libraries
QT += core

# Preprocessor defines
DEFINES += QT_DEPRECATED_WARNINGS

# Used to ignore useless warnings from OpenCL
QMAKE_CXXFLAGS += -Wno-ignored-attributes

# Source files
SOURCES += \
   analyticfactory.cpp \
   ccmatrix.cpp \
   correlationmatrix.cpp \
   datafactory.cpp \
   exportcorrelationmatrix_input.cpp \
   exportcorrelationmatrix.cpp \
   exportexpressionmatrix_input.cpp \
   exportexpressionmatrix.cpp \
   expressionmatrix.cpp \
   extract_input.cpp \
   extract.cpp \
   importcorrelationmatrix_input.cpp \
   importcorrelationmatrix.cpp \
   importexpressionmatrix_input.cpp \
   importexpressionmatrix.cpp \
   pairwise_clustering.cpp \
   pairwise_correlation.cpp \
   pairwise_gmm.cpp \
   pairwise_index.cpp \
   pairwise_kmeans.cpp \
   pairwise_linalg.cpp \
   pairwise_matrix.cpp \
   pairwise_pearson.cpp \
   pairwise_spearman.cpp \
   rmt_input.cpp \
   rmt.cpp \
   similarity_input.cpp \
   similarity_opencl_fetchpair.cpp \
   similarity_opencl_gmm.cpp \
   similarity_opencl_kmeans.cpp \
   similarity_opencl_pearson.cpp \
   similarity_opencl_spearman.cpp \
   similarity_opencl_worker.cpp \
   similarity_opencl.cpp \
   similarity_resultblock.cpp \
   similarity_serial.cpp \
   similarity_workblock.cpp \
   similarity.cpp

# Header files
HEADERS += \
   analyticfactory.h \
   ccmatrix.h \
   correlationmatrix.h \
   datafactory.h \
   exportcorrelationmatrix_input.h \
   exportcorrelationmatrix.h \
   exportexpressionmatrix_input.h \
   exportexpressionmatrix.h \
   expressionmatrix.h \
   extract_input.h \
   extract.h \
   importcorrelationmatrix_input.h \
   importcorrelationmatrix.h \
   importexpressionmatrix_input.h \
   importexpressionmatrix.h \
   pairwise_clustering.h \
   pairwise_correlation.h \
   pairwise_gmm.h \
   pairwise_index.h \
   pairwise_kmeans.h \
   pairwise_linalg.h \
   pairwise_matrix.h \
   pairwise_pearson.h \
   pairwise_spearman.h \
   rmt_input.h \
   rmt.h \
   similarity_input.h \
   similarity_opencl_fetchpair.h \
   similarity_opencl_gmm.h \
   similarity_opencl_kmeans.h \
   similarity_opencl_pearson.h \
   similarity_opencl_spearman.h \
   similarity_opencl_worker.h \
   similarity_opencl.h \
   similarity_resultblock.h \
   similarity_serial.h \
   similarity_workblock.h \
   similarity.h
