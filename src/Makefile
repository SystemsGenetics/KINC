.RECIPEPREFIX := +
MAKEFLAGS := $(MAKEFLAGS) --no-print-directory
CXX := g++
AR := ar
GENCL := ../../ACE/bin/mkcl.sh
NAME := kinc
FILES := KINC.files



build := ../build/
run := ../bin/

acxxflags := $(CXXFLAGS) -g -std=c++11 -I../../ACE
aldflags := $(LDFLAGS)
aldlibs := $(LDLIBS) -lOpenCL -L../../ACE/lib -lace -lgsl -lgslcblas
ucxxflags := $(acxxflags) -D UNIT_TEST

progf := $(run)$(NAME)
unitf := $(run)test

src := $(filter %.cpp,$(shell cat $(FILES)))
dpds := $(addprefix $(build),$(src:%.cpp=%.d))
objs := $(addprefix $(build),$(src:%.cpp=%.o))
dobjs := $(addprefix $(build),$(src:%.cpp=%.do))

.PHONY: clean depclean all program unit



all: program unit
program: $(progf)
unit: $(unitf)

ifndef DIGNORE
include $(dpds)
endif

$(progf): $(objs) $(dpds)
+@echo "Building program."
+@$(CXX) $(objs) $(aldflags) $(aldlibs) -o $@

$(unitf): $(dobjs) $(dpds)
+@echo "Building unit tests."
+@$(CXX) $(dobjs) $(aldflags) $(aldlibs) -o $@

depend: $(dpds)
+@echo Done.

$(build)%.o : %.cpp
+@echo "Building object $(*F)."
+@$(CXX) $(acxxflags) -c $< -o $@

$(build)%.do : %.cpp
+@echo "Building unit test object $(*F)."
+@$(CXX) $(ucxxflags) -c $< -o $@

%.cl.h: %.cl
+@echo "Building kernel $(*F)."
+@$(GENCL) $<

$(build)%.d: %.cpp
+@echo "Building depend $(*F)."
+@echo -n "$@ $(build)" > $@
+@$(CXX) $(acxxflags) -MM $< >> $@

clean:
+@echo "Cleaning all."
+@rm -f $(build)*.o $(build)*.do $(unitf)

depclean:
+@echo "Cleaning all dependency files."
+@rm -f $(build)*.d
