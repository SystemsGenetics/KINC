
MAKE = make

BUILD_CLI = build-cli
BUILD_GUI = build-gui
BUILD_TESTS = build-tests
SRC = src
TESTS = tests
BINS = kinc-cli kinc-tests

all: $(BINS)

kinc-cli: $(SRC)/*.h $(SRC)/*.cpp $(SRC)/opencl/*.cl
	cd $(BUILD_CLI) && GUI=0 qmake ../$(SRC)
	+$(MAKE) -C $(BUILD_CLI)
	cp $(BUILD_CLI)/kinc $@

kinc-gui: $(SRC)/*.h $(SRC)/*.cpp $(SRC)/opencl/*.cl
	cd $(BUILD_GUI) && GUI=1 qmake ../$(SRC)
	+$(MAKE) -C $(BUILD_GUI)
	cp $(BUILD_GUI)/kinc $@

kinc-tests: $(TESTS)/*.h $(TESTS)/*.cpp
	cd $(BUILD_TESTS) && qmake ../$(TESTS)
	+$(MAKE) -C $(BUILD_TESTS)
	cp $(BUILD_TESTS)/tests $@

clean:
	rm -f $(BUILD_CLI)/* $(BUILD_GUI)/* $(BUILD_TESTS)/* $(BINS)
