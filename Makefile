
MAKE = make

BUILD = build
SRC = src
BINS = kinc

all: $(BINS)

kinc: $(SRC)/*.h $(SRC)/*.cpp $(SRC)/opencl/*.cl
	cd $(BUILD) && qmake ../src
	+$(MAKE) -C $(BUILD)
	cp $(BUILD)/kinc .

clean:
	rm -f $(BUILD)/* $(BINS)
