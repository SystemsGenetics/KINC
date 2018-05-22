# General build variables
TARGET = tests
TEMPLATE = app
CONFIG += c++11 debug

# Qt libraries
QT += core testlib

# ACE and other external libraries
unix {
	LIBS += -lOpenCL -L/usr/local/lib64/ -L$$(HOME)/software/lib -lacecore -lgsl -lgslcblas -llapack -llapacke
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
	../src/pairwise_index.cpp \
	../src/pairwise_matrix.cpp \
	../src/pairwise_clustering.cpp \
	../src/pairwise_gmm.cpp \
	../src/pairwise_kmeans.cpp \
	../src/pairwise_linalg.cpp \
	../src/pairwise_correlation.cpp \
	../src/pairwise_pearson.cpp \
	../src/pairwise_spearman.cpp \
	../src/ccmatrix.cpp \
	../src/extract.cpp \
	testclustermatrix.cpp \
	testcorrelationmatrix.cpp \
	testexportcorrelationmatrix.cpp \
	testexportexpressionmatrix.cpp \
	testexpressionmatrix.cpp \
	testimportcorrelationmatrix.cpp \
	testimportexpressionmatrix.cpp \
	testrmt.cpp \
	testsimilarity.cpp \
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
	../src/pairwise_index.h \
	../src/pairwise_matrix.h \
	../src/pairwise_clustering.h \
	../src/pairwise_gmm.h \
	../src/pairwise_kmeans.h \
	../src/pairwise_linalg.h \
	../src/pairwise_correlation.h \
	../src/pairwise_pearson.h \
	../src/pairwise_spearman.h \
	../src/ccmatrix.h \
	../src/extract.h \
	testclustermatrix.h \
	testcorrelationmatrix.h \
	testexportcorrelationmatrix.h \
	testexportexpressionmatrix.h \
	testexpressionmatrix.h \
	testimportcorrelationmatrix.h \
	testimportexpressionmatrix.h \
	testrmt.h \
	testsimilarity.h
