
DEBUG ?= 0
INSTALL_PREFIX ?= /usr/local
LINK_MPI_CXX ?= 

MAKE = make
BUILD = build
SRC = src

QMAKEFLAGS = PREFIX=$(INSTALL_PREFIX)  LINK_MPI_CXX=$(LINK_MPI_CXX)

ifeq ($(DEBUG), 1)
QMAKEFLAGS += CONFIG+=debug
endif

all: kinc

$(BUILD):
	mkdir -p $(BUILD)

kinc: $(SRC)/**/*.h $(SRC)/**/*.cpp $(SRC)/**/*.cl $(SRC)/**/*.cu | $(BUILD)
	rm -rf $(BUILD)/cli $(BUILD)/gui
	cd $(BUILD) && qmake ../$(SRC)/KINC.pro $(QMAKEFLAGS)
	+$(MAKE) -C $(BUILD)

install: kinc
	+$(MAKE) -C $(BUILD) qmake_all
	+$(MAKE) -C $(BUILD) install
	cp bin/* $(INSTALL_PREFIX)/bin

clean:
	rm -rf $(BUILD)/* $(BUILD)/.qmake.stash
