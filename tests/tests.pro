# General build variables
TARGET = tests
TEMPLATE = app
CONFIG += c++11 debug

# Qt libraries
QT += core testlib

# ACE and other external libraries
unix {
   LIBS += -lOpenCL -L/usr/local/lib64/ -L$$(HOME)/software/lib -lacecore -lgsl -lgslcblas
   INCLUDEPATH += $$(HOME)'/software/include'
}

win32 {
   LIBS += -L$$(OCL_ROOT)'/lib/x86' -lopencl -L'C:/ACE/lib' -lacecore0
   INCLUDEPATH += 'C:/ACE/include' $$(OCL_ROOT)'/include'
   DEPENDPATH += $$(OCL_ROOT)'/include'
}

INCLUDEPATH += ../src

# Preprocessor defines
DEFINES += QT_DEPRECATED_WARNINGS

# Source files
SOURCES += \
   ../src/analyticfactory.cpp \
   ../src/datafactory.cpp \
   ../src/expressionmatrix.cpp \
   ../src/importexpressionmatrix.cpp \
   ../src/exportexpressionmatrix.cpp \
   ../src/correlationmatrix.cpp \
   ../src/importcorrelationmatrix.cpp \
   ../src/exportcorrelationmatrix.cpp \
   ../src/similarity.cpp \
   ../src/rmt.cpp \
   ../src/genepair_index.cpp \
   ../src/genepair_base.cpp \
   ../src/genepair_clustering.cpp \
   ../src/genepair_gmm.cpp \
   ../src/genepair_kmeans.cpp \
   ../src/genepair_linalg.cpp \
   ../src/genepair_correlation.cpp \
   ../src/genepair_pearson.cpp \
   ../src/genepair_spearman.cpp \
   ../src/ccmatrix.cpp \
   ../src/extract.cpp \
   testexpressionmatrix.cpp \
   main.cpp

HEADERS += \
   ../src/analyticfactory.h \
   ../src/datafactory.h \
   ../src/expressionmatrix.h \
   ../src/importexpressionmatrix.h \
   ../src/exportexpressionmatrix.h \
   ../src/correlationmatrix.h \
   ../src/importcorrelationmatrix.h \
   ../src/exportcorrelationmatrix.h \
   ../src/similarity.h \
   ../src/rmt.h \
   ../src/genepair_index.h \
   ../src/genepair_base.h \
   ../src/genepair_clustering.h \
   ../src/genepair_gmm.h \
   ../src/genepair_kmeans.h \
   ../src/genepair_linalg.h \
   ../src/genepair_correlation.h \
   ../src/genepair_pearson.h \
   ../src/genepair_spearman.h \
   ../src/ccmatrix.h \
   ../src/extract.h \
   testexpressionmatrix.h
