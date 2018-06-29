
# Default settings for MPI CXX include
isEmpty(MPICXX) { MPICXX = "yes" }

# Versions
MAJOR_VERSION = 0
MINOR_VERSION = 0
REVISION = 999
VERSION = $${MAJOR_VERSION}.$${MINOR_VERSION}.$${REVISION}

# Version compiler defines
DEFINES += \
    QT_DEPRECATED_WARNINGS \
    MAJOR_VERSION=$${MAJOR_VERSION} \
    MINOR_VERSION=$${MINOR_VERSION} \
    REVISION=$${REVISION}

# Basic settings
QT += core
TEMPLATE = app
QMAKE_CXX = mpic++
CONFIG += c++11

# Compiler defines
DEFINES += QT_DEPRECATED_WARNINGS

# External libraries
LIBS += -lmpi
equals(MPICXX,"yes") { LIBS += -lmpi_cxx }
LIBS += -lacecore -lOpenCL -lgsl -lgslcblas -L$${PWD}/../../build/libs -lkinccore

# Used to ignore useless warnings with OpenCL
QMAKE_CXXFLAGS += -Wno-ignored-attributes

# Source files
SOURCES += \
    ../main.cpp \

# Resource files
RESOURCES += \
    ../opencl.qrc
