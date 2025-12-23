
# Include common settings
include (../KINC.pri)

# Basic settings
TARGET = kinc
TEMPLATE = app

# External libraries
LIBS += -lacecli

# Compiler defines
DEFINES += GUI=0

# Source files
SOURCES += \
    ../main.cpp

# Installation instructions
isEmpty(PREFIX) { PREFIX = /usr/local }
program.path = $${PREFIX}/bin
program.files = $${OUT_PWD}/$${TARGET}
INSTALLS += program
