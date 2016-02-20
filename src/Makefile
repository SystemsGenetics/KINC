.RECIPEPREFIX := +
MAKEFLAGS := $(MAKEFLAGS) --no-print-directory
CXX := g++
NAME := kinc
FILES := KINC.files



build := ../build/
run := ../run/

acxxflags := $(CXXFLAGS) -g -std=c++11 $(shell sdl-config --cflags)
aldflags := $(LDFLAGS)
aldlibs := $(LDLIBS) -lOpenCL $(shell sdl-config --libs)

progf := $(run)$(NAME)

allsrc := $(filter %.cpp,$(shell cat $(FILES)))
utest := $(filter-out main.cpp,$(allsrc))
program := $(filter-out %.t.cpp,$(allsrc))

alldpds := $(addprefix $(build),$(allsrc:%.cpp=%.d))

dpds := $(addprefix $(build),$(program:%.cpp=%.d))
objs := $(addprefix $(build),$(program:%.cpp=%.o))

udpds := $(addprefix $(build),$(utest:%.cpp=%.d))
uobjs := $(addprefix $(build),$(utest:%.cpp=%.o))

plugin := plugins/plugins.a



.PHONY: clean all program test check plugins plugclean

all: program test
program: $(progf)
test: $(run)unit

include $(alldpds)

$(progf): $(objs) $(dpds) $(plugin)
+@echo "Building program."
+@$(CXX) $(objs) $(aldflags) $(aldlibs) -o $@

$(run)unit: $(uobjs) $(udpds) $(plugin)
+@echo "Building unit test."
+@$(CXX) $(uobjs) $(aldflags) $(aldlibs) -o $@

depend: $(alldpds)
+@echo Done.

plugins:
+@cd plugins && $(MAKE)

$(plugin):
+@cd plugins && $(MAKE)

$(build)%.o : %.cpp
+@echo "Building object $@"
+@$(CXX) $(acxxflags) -c $< -o $(build)$@

$(build)%.d: %.cpp
+@echo "Building depend $@"
+@echo -n "$@ $(build)" > $@
+@$(CXX) $(acxxflags) -MM $< >> $@

check: test
+@echo "Running unit tests..."
+@cd $(run) && ./unit

clean:
+@echo "Cleaning all."
+@rm -f $(build)*.o $(progf)

depclean:
+@echo "Cleaning all dependency files."
+@rm -f $(build)*.d

plugclean:
+@cd plugins && $(MAKE) plugclean
