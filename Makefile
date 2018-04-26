
MAKE = make

BUILD_CLI = build-cli
BUILD_GUI = build-gui
BUILD_TESTS = build-tests
SRC = src
TESTS = tests
BINS = kinc-cli kinc-gui $(BUILD_TESTS)/tests

all: $(BINS)

kinc-cli: $(SRC)/*.h $(SRC)/*.cpp $(SRC)/opencl/*.cl
	cd $(BUILD_CLI) && GUI=0 qmake ../$(SRC)
	+$(MAKE) -C $(BUILD_CLI)
	cp $(BUILD_CLI)/kinc $@

kinc-gui: $(SRC)/*.h $(SRC)/*.cpp $(SRC)/opencl/*.cl
	cd $(BUILD_GUI) && GUI=1 qmake ../$(SRC)
	+$(MAKE) -C $(BUILD_GUI)
	cp $(BUILD_GUI)/kinc $@

$(BUILD_TESTS)/tests: $(TESTS)/*.h $(TESTS)/*.cpp
	cd $(BUILD_TESTS) && qmake ../$(TESTS)
	+$(MAKE) -C $(BUILD_TESTS)

clean:
	rm -f $(BUILD_CLI)/* $(BUILD_GUI)/* $(BUILD_TESTS)/* $(BINS)
