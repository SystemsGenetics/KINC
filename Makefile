
MAKE = make

BUILD_CLI = build-cli
BUILD_GUI = build-gui
SRC = src
BINS = kinc-cli kinc-gui

all: $(BINS)

kinc-cli: $(SRC)/*.h $(SRC)/*.cpp $(SRC)/opencl/*.cl
	cd $(BUILD_CLI) && GUI=0 qmake ../src
	+$(MAKE) -C $(BUILD_CLI)
	cp $(BUILD_CLI)/kinc $@

kinc-gui: $(SRC)/*.h $(SRC)/*.cpp $(SRC)/opencl/*.cl
	cd $(BUILD_GUI) && GUI=1 qmake ../src
	+$(MAKE) -C $(BUILD_GUI)
	cp $(BUILD_GUI)/kinc $@

clean:
	rm -f $(BUILD_CLI)/* $(BUILD_GUI)/* $(BINS)
