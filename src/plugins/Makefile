.RECIPEPREFIX := +
AR=ar
CXX := g++



acxxflags := $(CXXFLAGS) -g -std=c++11

data := $(shell find data -mindepth 1 -maxdepth 1 -type d)
data := $(patsubst data/%,%.dat.a,$(data))

analytic := $(shell find analytic -mindepth 1 -maxdepth 1 -type d)
analytic := $(patsubst analytic/%,%.alc.a,$(analytic))

binds := $(strip $(data) $(analytic))
binding := analytic/plugin.o data/plugin.o

arflags := rcs
arscript := ar.sh



.PHONY: clean

plugins.a: $(binds) $(binding)
+@echo "Building plugins library."
+@echo "CREATE $@" > $(arscript)
+@for a in $(binds); do (echo "ADDLIB $$a" >> $(arscript)); done
+@echo "ADDMOD $(binding)" >> $(arscript)
+@echo "SAVE" >> $(arscript)
+@echo "END" >> $(arscript)
+@ar -M < $(arscript)

%.dat.a: data/%
+@echo "Building data plugin $*"
+@cd $< && $(MAKE) TARGET=../../$@

%.alc.a: analytic/%
+@echo "Building analytic plugin $*"
+@cd $< && $(MAKE) TARGET=../../$@

analytic/plugin.cpp: analytic/gen.sh
+@echo "Generating analytic binder."
+@cd analytic && ./gen.sh

data/plugin.cpp: analytic/gen.sh
+@echo "Generating data binder."
+@cd data && ./gen.sh

%.o:%.cpp
+@echo "Building binder $@"
+@$(CXX) $(acxxflags) -c $< -o $@

clean:
+@echo "Cleaning all plugins."
+@rm -f *.a data/*.o data/*.cpp analytic/*.o analytic/*.cpp
