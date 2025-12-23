
# Include common settings
include (../KINC.pri)

# Basic settings
QT += testlib
TARGET = kinc-tests
TEMPLATE = app
CONFIG += debug

# Source files
SOURCES += \
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
	testclustermatrix.h \
	testcorrelationmatrix.h \
	testexportcorrelationmatrix.h \
	testexportexpressionmatrix.h \
	testexpressionmatrix.h \
	testimportcorrelationmatrix.h \
	testimportexpressionmatrix.h \
	testrmt.h \
	testsimilarity.h

# Installation instructions
isEmpty(PREFIX) { PREFIX = /usr/local }
program.path = $${PREFIX}/bin
program.files = $${OUT_PWD}/$${TARGET}
INSTALLS += program
