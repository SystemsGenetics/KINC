
# Include common settings
include (../KINC.pri)

# Basic settings
TARGET = kinc

# External libraries
LIBS += -laceconsole

# Compiler defines
DEFINES += GUI=0

# Installation instructions
isEmpty(PREFIX) { PREFIX = /usr/local }
program.path = $${PREFIX}/bin
program.files = $${PWD}/../../build/cli/$${TARGET}
INSTALLS += program
