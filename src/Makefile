.RECIPEPREFIX := +
MAKEFLAGS := $(MAKEFLAGS) --no-print-directory
CXX := g++
NAME := kinc
FILES := KINC.files



build := ../build/
run := ../run/

acxxflags := $(CXXFLAGS) -g -std=c++11
aldflags := $(LDFLAGS)
aldlibs := $(LDLIBS) -lOpenCL

progf := $(run)$(NAME)

allsrc := $(filter %.cpp,$(shell cat $(FILES)))
utest := $(filter-out main.cpp,$(allsrc))
program := $(filter-out %.t.cpp,$(allsrc))

alldpds := $(addprefix $(build),$(allsrc:%.cpp=%.d))

dpds := $(addprefix $(build),$(program:%.cpp=%.d))
objs := $(addprefix $(build),$(program:%.cpp=%.o))

udpds := $(addprefix $(build),$(utest:%.cpp=%.d))
uobjs := $(addprefix $(build),$(utest:%.cpp=%.o))

plugins := plugins/plugins.a



.PHONY: clean all program test check plugins pclean

all: program test
program: $(progf)
test: $(run)unit

include $(alldpds)

$(progf): $(objs) $(dpds) $(plugins)
+@echo "Building program."
+@$(CXX) $(objs) $(plugins) $(aldflags) $(aldlibs) -o $@

$(run)unit: $(uobjs) $(udpds)
+@echo "Building unit test."
+@$(CXX) $(uobjs) $(aldflags) $(aldlibs) -o $@

depend: $(alldpds)
+@echo Done.

plugins:
+@cd plugins && $(MAKE)

$(plugins):
+@cd plugins && $(MAKE)

$(build)%.o : %.cpp
+@echo "Building object $@"
+@$(CXX) $(acxxflags) -c $< -o $(build)$@

$(build)%.d: %.cpp
+@echo "Building depend $@"
+@echo -n "$@ $(build)" > $@
+@$(CXX) $(acxxflags) -MM $< >> $@

check: test
+@cd $(run)/tmp && ../unit

clean:
+@echo "Cleaning all."
+@rm -f $(build)*.o $(progf)

depclean:
+@echo "Cleaning all dependency files."
+@rm -f $(build)*.d

pclean:
+@cd plugins && $(MAKE) clean
