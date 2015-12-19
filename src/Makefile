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



progf := $(RUN)$(NAME)
program := $(SOURCES) main.cpp
utest := $(SOURCES) $(UNITS)
allsrc := $(SOURCES) $(UNITS) main.cpp
alldpds := $(addprefix $(BUILD),$(allsrc:%.cpp=%.d))
dpds := $(addprefix $(BUILD),$(program:%.cpp=%.d))
objs := $(addprefix $(BUILD),$(program:%.cpp=%.o))
udpds := $(addprefix $(BUILD),$(utest:%.cpp=%.d))
uobjs := $(addprefix $(BUILD),$(utest:%.cpp=%.o))



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
