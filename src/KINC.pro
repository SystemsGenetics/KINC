# General build variables
TARGET = kinc
TEMPLATE = app
CONFIG += c++11

GUI = $$(GUI)

isEmpty(GUI) {
   GUI = 0
}

# Qt libraries
QT += core

# external libraries
LIBS += -L$$(HOME)/software/lib -lacecore -laceconsole -lgsl -lgslcblas -llapack -llapacke
INCLUDEPATH += $$(HOME)/software/include

# HACK
INCLUDEPATH += $$(HOME)/software/include/ace

# Preprocessor defines
DEFINES += QT_DEPRECATED_WARNINGS

equals(GUI, 1) {
   DEFINES += ACE_GUI
}

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
   main.cpp \
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
   # similarity.cpp

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
   # similarity.h

RESOURCES += \
   opencl.qrc
