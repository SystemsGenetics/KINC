
# Versions
KINC_MAJOR_VERSION = 3
KINC_MINOR_VERSION = 3
KINC_REVISION = 0

VERSION = $${KINC_MAJOR_VERSION}.$${KINC_MINOR_VERSION}.$${KINC_REVISION}

# Version compiler defines
DEFINES += \
    QT_DEPRECATED_WARNINGS \
    KINC_MAJOR_VERSION=$${KINC_MAJOR_VERSION} \
    KINC_MINOR_VERSION=$${KINC_MINOR_VERSION} \
    KINC_REVISION=$${KINC_REVISION}

# Basic settings
QT += core
QMAKE_CXX = mpic++
CONFIG += c++11

# Compiler flags
QMAKE_CXXFLAGS += -Wno-ignored-attributes

# Preprocessor defines
DEFINES += \
    QT_DEPRECATED_WARNINGS \
    KINC_MAJOR_VERSION=$${KINC_MAJOR_VERSION} \
    KINC_MINOR_VERSION=$${KINC_MINOR_VERSION} \
    KINC_REVISION=$${KINC_REVISION}

# Default settings for optional linker flags
isEmpty(LINK_LAPACKE) { LINK_LAPACKE = 1 }
isEmpty(LINK_MPI_CXX) { LINK_MPI_CXX = 1 }

# Set CUDA directory if not already set, first from environment variable if it exists or a default
# if it does not
isEmpty(CUDADIR) {
    CUDADIR = $$(CUDADIR)
    isEmpty(CUDADIR) {
        CUDADIR = /usr/local/cuda
    }
}

INCLUDEPATH += $${CUDADIR}/include
DEPENDPATH += $${CUDADIR}/include

INCLUDEPATH += /usr/include/openblas

# External libraries
LIBS += \
    -L$${PWD}/../build/libs -lkinccore \
    -lacecore \
    -lgsl -lopenblas \
    -L$${CUDADIR}/lib64 -lcuda -lnvrtc -lcusolver -fopenmp \
    -lOpenCL -lmpi

equals(LINK_LAPACKE, 1) { LIBS += -llapacke }
equals(LINK_MPI_CXX, 1) { LIBS += -lmpi_cxx }

# Resource files
RESOURCES += \
    ../resources.qrc
