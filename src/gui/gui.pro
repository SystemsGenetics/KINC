
# Include common settings
include (../KINC.pri)

# Basic settings
QT += gui widgets
TARGET = qkinc
TEMPLATE = app

# External libraries
LIBS += -lacegui

# Compiler defines
DEFINES += GUI=1

# Source files
SOURCES += \
    ../main.cpp

# Installation instructions
isEmpty(PREFIX) { PREFIX = /usr/local }
program.path = $${PREFIX}/bin
program.files = $${OUT_PWD}/$${TARGET}
INSTALLS += program
