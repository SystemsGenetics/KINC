# General build variables
TARGET = tests
TEMPLATE = app
CONFIG += c++11 debug

# Qt libraries
QT += core testlib

# external libraries
LIBS += -lOpenCL -L/usr/local/lib64/ -L$$(HOME)/software/lib -lacecore -lgsl -lgslcblas -llapack -llapacke
INCLUDEPATH += $$(HOME)/software/include
INCLUDEPATH += ../src

# HACK
INCLUDEPATH += $$(HOME)/software/include/ace

# Preprocessor defines
DEFINES += QT_DEPRECATED_WARNINGS

# Source files
SOURCES += \
	../src/analyticfactory.cpp \
	../src/ccmatrix.cpp \
	../src/correlationmatrix.cpp \
	../src/datafactory.cpp \
	# ../src/exportcorrelationmatrix.cpp \
	../src/exportexpressionmatrix_input.cpp \
	../src/exportexpressionmatrix.cpp \
	../src/expressionmatrix.cpp \
	# ../src/extract.cpp \
	# ../src/importcorrelationmatrix.cpp \
	../src/importexpressionmatrix_input.cpp \
	../src/importexpressionmatrix.cpp \
	../src/pairwise_clustering.cpp \
	../src/pairwise_correlation.cpp \
	../src/pairwise_gmm.cpp \
	../src/pairwise_index.cpp \
	../src/pairwise_kmeans.cpp \
	../src/pairwise_linalg.cpp \
	../src/pairwise_matrix.cpp \
	../src/pairwise_pearson.cpp \
	../src/pairwise_spearman.cpp \
	# ../src/rmt.cpp \
	# ../src/similarity.cpp \
	testclustermatrix.cpp \
	testcorrelationmatrix.cpp \
	# testexportcorrelationmatrix.cpp \
	# testexportexpressionmatrix.cpp \
	testexpressionmatrix.cpp \
	# testimportcorrelationmatrix.cpp \
	# testimportexpressionmatrix.cpp \
	# testrmt.cpp \
	# testsimilarity.cpp \
	main.cpp

HEADERS += \
	../src/analyticfactory.h \
	../src/ccmatrix.h \
	../src/correlationmatrix.h \
	../src/datafactory.h \
	../src/expressionmatrix.h \
	# ../src/extract.h \
	# ../src/exportcorrelationmatrix.h \
	../src/exportexpressionmatrix_input.h \
	../src/exportexpressionmatrix.h \
	# ../src/importcorrelationmatrix.h \
	../src/importexpressionmatrix_input.h \
	../src/importexpressionmatrix.h \
	../src/pairwise_clustering.h \
	../src/pairwise_correlation.h \
	../src/pairwise_gmm.h \
	../src/pairwise_index.h \
	../src/pairwise_kmeans.h \
	../src/pairwise_linalg.h \
	../src/pairwise_matrix.h \
	../src/pairwise_pearson.h \
	../src/pairwise_spearman.h \
	# ../src/rmt.h \
	# ../src/similarity.h \
	testclustermatrix.h \
	testcorrelationmatrix.h \
	# testexportcorrelationmatrix.h \
	# testexportexpressionmatrix.h \
	testexpressionmatrix.h
	# testimportcorrelationmatrix.h \
	# testimportexpressionmatrix.h \
	# testrmt.h \
	# testsimilarity.h
