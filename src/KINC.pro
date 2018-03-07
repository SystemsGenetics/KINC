# General build variables
TARGET = kinc
TEMPLATE = app
CONFIG += c++11
GUI = 1

# Qt libraries
QT += core

equals(GUI, 1) {
   QT += gui

   greaterThan(QT_MAJOR_VERSION, 4) {
      QT += widgets
   }
}

# ACE and other external libraries
unix {
   LIBS += -lOpenCL -L/usr/local/lib64/ -L$$(HOME)/software/lib -lacecore -lgsl -lgslcblas
   INCLUDEPATH += $$(HOME)'/software/include'

   equals(GUI, 1) {
      LIBS += -lacegui
   }
   else {
      LIBS += -laceconsole
   }
}

win32 {
   LIBS += -L$$(OCL_ROOT)'/lib/x86' -lopencl -L'C:/ACE/lib' -lacecore0
   INCLUDEPATH += 'C:/ACE/include' $$(OCL_ROOT)'/include'
   DEPENDPATH += $$(OCL_ROOT)'/include'

   equals(GUI, 1) {
      LIBS += -lacegui0
   }
   else {
      LIBS += -laceconsole0
   }
}

# Preprocessor defines
DEFINES += QT_DEPRECATED_WARNINGS

equals(GUI, 1) {
   DEFINES += ACE_GUI
}

# Source files
SOURCES += \
    main.cpp \
    analyticfactory.cpp \
    datafactory.cpp \
    expressionmatrix.cpp \
    importexpressionmatrix.cpp \
    exportexpressionmatrix.cpp \
    correlationmatrix.cpp \
    importcorrelationmatrix.cpp \
    exportcorrelationmatrix.cpp \
    spearman.cpp \
    rmt.cpp \
    pearson.cpp \
    correlationbase.cpp \
    genepair_vector.cpp \
    genepair_base.cpp \
    genepair_gmm.cpp \
    genepair_kmeans.cpp \
    genepair_linalg.cpp \
    ccmatrix.cpp \
    kmeans.cpp \
    gmm.cpp \
    extract.cpp

HEADERS += \
    analyticfactory.h \
    datafactory.h \
    expressionmatrix.h \
    importexpressionmatrix.h \
    exportexpressionmatrix.h \
    correlationmatrix.h \
    importcorrelationmatrix.h \
    exportcorrelationmatrix.h \
    spearman.h \
    rmt.h \
    pearson.h \
    correlationbase.h \
    genepair_vector.h \
    genepair_base.h \
    genepair_gmm.h \
    genepair_kmeans.h \
    genepair_linalg.h \
    ccmatrix.h \
    kmeans.h \
    gmm.h \
    extract.h

RESOURCES += \
    opencl.qrc
