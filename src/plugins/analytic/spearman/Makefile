.RECIPEPREFIX := +
MAKEFLAGS := $(MAKEFLAGS) --no-print-directory
CXX := g++
NAME := kinc
FILES := spearman.files
GENCL := ../../../../bin/mkcl.sh



build := ./build/

acxxflags := $(CXXFLAGS) -g -std=c++11
aldflags := $(LDFLAGS)
aldlibs := $(LDLIBS) -lOpenCL

progf := $(run)$(NAME)

src := $(filter %.cpp,$(shell cat $(FILES)))

dpds := $(addprefix $(build),$(src:%.cpp=%.d))
objs := $(addprefix $(build),$(src:%.cpp=%.o))



.PHONY: clean all program test check plugins pclean

all: $(TARGET)

include $(dpds)

$(TARGET): $(objs) $(dpds)
+@echo "Building plugin."
+@ar rc $@ $(objs)

depend: $(dpds)
+@echo Done.

$(build)%.o : %.cpp
+@echo "Building object $@"
+@$(CXX) $(acxxflags) -c $< -o $@

$(build)%.d: %.cpp
+@echo "Building depend $@"
+@echo -n "$@ $(build)" > $@
+@$(CXX) $(acxxflags) -MM $< >> $@

%.cl.h: %.cl
+@echo "Building kernel $@"
+@$(GENCL) $<

clean:
+@echo "Cleaning all."
+@rm -f $(build)*.o $(progf)

depclean:
+@echo "Cleaning all dependency files."
+@rm -f $(build)*.d
