
# Include common settings
include (../KINC.pri)

# Basic settings
QT += gui widgets
TARGET = qkinc

# External libraries
LIBS += -lacegui

# Compiler defines
DEFINES += GUI=1

# Installation instructions
isEmpty(PREFIX) { PREFIX = /usr/local }
program.path = $${PREFIX}/bin
program.files = $${PWD}/../../build/gui/$${TARGET}
INSTALLS += program
