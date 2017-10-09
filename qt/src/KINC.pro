#-------------------------------------------------
#
# Project created by QtCreator 2017-08-24T19:05:56
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = KINC
TEMPLATE = app
CONFIG += c++11

unix: LIBS += -lOpenCL -L/usr/local/lib64/ -lacecore -lacegui
win32: LIBS += -L$$(OCL_ROOT)'/lib/x86' -lopencl -L'C:/ACE/lib' -lacecore0 -lacegui0
win32: INCLUDEPATH += 'C:/ACE/include' $$(OCL_ROOT)'/include'
win32: DEPENDPATH += $$(OCL_ROOT)'/include'

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
    main.cpp \
    analyticfactory.cpp \
    datafactory.cpp \
    expressionmatrix.cpp \
    importexpressionmatrix.cpp \
    correlationmatrix.cpp \
    spearman.cpp \
    rmt.cpp

HEADERS += \
    analyticfactory.h \
    datafactory.h \
    expressionmatrix.h \
    importexpressionmatrix.h \
    correlationmatrix.h \
    spearman.h \
    rmt.h

RESOURCES += \
    opencl.qrc

