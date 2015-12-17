.RECIPEPREFIX := +
BUILD := ../build/
RUN := ../run/



CXX := g++
CXXFLAGS := -g `sdl-config --libs`
LDFLAGS := -lOpenCL `sdl-config --libs`



NAME := kinc

SOURCES := \
console.cpp \
terminal.cpp \
consolestream.cpp \
exception.cpp

UNITS := \
unit.cpp



progf=$(RUN)$(NAME)
program=$(SOURCES) main.cpp
utest=$(SOURCES) $(UNITS)
allsrc=$(SOURCES) $(UNITS) main.cpp
alldpds=$(patsubst %.cpp,$(BUILD)%.d,$(allsrc))
dpds=$(patsubst %.cpp,$(BUILD)%.d,$(program))
objs=$(patsubst %.cpp,$(BUILD)%.o,$(program))
udpds=$(patsubst %.cpp,$(BUILD)%.d,$(utest))
uobjs=$(patsubst %.cpp,$(BUILD)%.o,$(utest))



all: program test
program: $(progf)
test: $(RUN)unit

include $(alldpds)

$(progf): $(objs) $(dpds)
+@echo "Building program."
+@$(CXX) $(objs) $(LDFLAGS) -o $@

$(RUN)unit: $(uobjs) $(udpds)
+@echo "Building unit test."
+@$(CXX) $(uobjs) $(LDFLAGS) -o $@

depend: $(alldpds)
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
+@rm -f $(BUILD)*.o $(progf)
