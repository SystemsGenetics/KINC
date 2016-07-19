.RECIPEPREFIX := +
MAKEFLAGS := $(MAKEFLAGS) --no-print-directory
CXX := g++
AR := ar
GENCL := ../../ACE/bin/mkcl.sh
NAME := kinc
FILES := KINC.files



build := ../build/
run := ../run/

acxxflags := $(CXXFLAGS) -g -std=c++11 -I../../ACE/src
aldflags := $(LDFLAGS)
aldlibs := $(LDLIBS) -lOpenCL -L../../ACE/lib -lace

progf := $(run)$(NAME)

src := $(filter %.cpp,$(shell cat $(FILES)))
dpds := $(addprefix $(build),$(src:%.cpp=%.d))
objs := $(addprefix $(build),$(src:%.cpp=%.o))

.PHONY: clean depclean all



all: $(progf)

include $(dpds)

$(progf): $(objs) $(dpds)
+@echo "Building program."
+@$(CXX) $(objs) $(aldflags) $(aldlibs) -o $@

depend: $(dpds)
+@echo Done.

$(build)%.o : %.cpp
+@echo "Building object $(*F)."
+@$(CXX) $(acxxflags) -c $< -o $@

%.cl.h: %.cl
+@echo "Building kernel $(*F)."
+@$(GENCL) $<

$(build)%.d: %.cpp
+@echo "Building depend $(*F)."
+@echo -n "$@ $(build)" > $@
+@$(CXX) $(acxxflags) -MM $< >> $@

clean:
+@echo "Cleaning all."
+@rm -f $(build)*.o $(unitf)

depclean:
+@echo "Cleaning all dependency files."
+@rm -f $(build)*.d
