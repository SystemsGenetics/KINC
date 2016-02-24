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
datahdr := plugins/data/plugin.h
analytichdr := plugins/analytic/plugin.h



.PHONY: clean all program test check plugins pclean

all: program test
program: $(progf)
test: $(run)unit

$(analytichdr):
+@echo "Generating analytic plugin header."
+@cd plugins/analytic && ./gen.sh

$(datahdr):
+@echo "Generating data plugin header."
+@cd plugins/data && ./gen.sh

include $(alldpds)

$(progf): $(objs) $(dpds) $(plugin)
+@echo "Building program."
+@$(CXX) $(objs) $(plugin) $(aldflags) $(aldlibs) -o $@

$(run)unit: $(uobjs) $(udpds) $(plugin)
+@echo "Building unit test."
+@$(CXX) $(uobjs) $(plugin) $(aldflags) $(aldlibs) -o $@

depend: $(alldpds)
+@echo Done.

plugins:
+@cd plugins && $(MAKE)

$(plugin):
+@cd plugins && $(MAKE)

$(build)%.o : %.cpp
+@echo "Building object $@"
+@$(CXX) $(acxxflags) -c $< -o $(build)$@

$(build)%.d: %.cpp $(datahdr) $(analytichdr)
+@echo "Building depend $@"
+@echo -n "$@ $(build)" > $@
+@$(CXX) $(acxxflags) -MM $< >> $@

check: test
+@cd $(run) && ./unit

clean:
+@echo "Cleaning all."
+@rm -f $(build)*.o $(progf)

depclean:
+@echo "Cleaning all dependency files."
+@rm -f $(build)*.d

pclean:
+@cd plugins && $(MAKE) clean
