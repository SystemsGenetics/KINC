
# Include common settings
include (../KINC.pri)

# Basic settings
TARGET = kinc

# External libraries
LIBS += -laceconsole

# Compiler defines
DEFINES += GUI=0

# Source files
SOURCES += \
    ../main.cpp

# Installation instructions
isEmpty(PREFIX) { PREFIX = /usr/local }
program.path = $${PREFIX}/bin
program.files = $${PWD}/../../build/cli/$${TARGET}
INSTALLS += program
