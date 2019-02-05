
# Versions
MAJOR_VERSION = 3
MINOR_VERSION = 2
REVISION = 2

VERSION = $${MAJOR_VERSION}.$${MINOR_VERSION}.$${REVISION}

# Basic settings
QT += core
QMAKE_CXX = mpic++
CONFIG += c++11

# Compiler flags
QMAKE_CXXFLAGS += -Wno-ignored-attributes

# Preprocessor defines
DEFINES += \
    QT_DEPRECATED_WARNINGS \
    MAJOR_VERSION=$${MAJOR_VERSION} \
    MINOR_VERSION=$${MINOR_VERSION} \
    REVISION=$${REVISION}

# Include directories
INCLUDEPATH += /usr/include/openblas

# Default settings for optional linker flags
isEmpty(LINK_LAPACKE) { LINK_LAPACKE = 1 }
isEmpty(LINK_MPI_CXX) { LINK_MPI_CXX = 1 }

# External libraries
LIBS += \
    -L$${PWD}/../build/libs -lkinccore \
    -lacecore \
    -lgsl -lopenblas \
    -lOpenCL -lmpi

equals(LINK_LAPACKE, 1) { LIBS += -llapacke }
equals(LINK_MPI_CXX, 1) { LIBS += -lmpi_cxx }

# Resource files
RESOURCES += \
    ../opencl.qrc
