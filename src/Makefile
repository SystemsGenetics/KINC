.RECIPEPREFIX=+
BUILD=../build/
RUN=../run/



CLPP=../../OpenCLPP
CXX=g++
CXXFLAGS=-g `sdl-config --libs` -I$(CLPP)/include
LDFLAGS=-L$(CLPP)/lib -lOpenCL -lclpp `sdl-config --libs`



NAME=kinc

SOURCES=\
console.cpp

UNITS=\
unit.cpp



PROGF=$(RUN)$(NAME)
PROGRAM=$(SOURCES) main.cpp
UTEST=$(SOURCES) $(UNITS)
ALLSRC=$(SOURCES) $(UNITS) main.cpp
ALLDPDS=$(patsubst %.cpp,$(BUILD)%.d,$(ALLSRC))
DPDS=$(patsubst %.cpp,$(BUILD)%.d,$(PROGRAM))
OBJS=$(patsubst %.cpp,$(BUILD)%.o,$(PROGRAM))
UDPDS=$(patsubst %.cpp,$(BUILD)%.d,$(UTEST))
UOBJS=$(patsubst %.cpp,$(BUILD)%.o,$(UTEST))



all: program test
program: $(PROGF)
test: $(RUN)unit

include $(ALLDPDS)

$(PROGF): $(OBJS) $(DPDS)
+@echo "Building program."
+@$(CXX) $(OBJS) $(LDFLAGS) -o $@

$(RUN)unit: $(UOBJS) $(UDPDS)
+@echo "Building unit test."
+@$(CXX) $(UOBJS) $(LDFLAGS) -o $@

depend: $(ALLDPDS)
+@echo Done.

$(BUILD)%.o : %.cpp
+@echo "Building object $@"
+@$(CXX) $(CXXFLAGS) -c $< -o $(BUILD)$@

$(BUILD)%.d: %.cpp
+@echo "Building depend $@"
+@echo -n "$@ $(BUILD)" > $@
+@$(CXX) -MM $< >> $@

.PHONY=clean
clean:
+@echo "Cleaning all."
+@rm -f $(BUILD)*.o $(PROGF)
