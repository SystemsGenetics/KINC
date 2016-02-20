.RECIPEPREFIX := +
AR=ar



data := $(shell find data -mindepth 1 -maxdepth 1 -type d)
analytic := $(shell find analytic -mindepth 1 -maxdepth 1 -type d)
data := $(patsubst data/%,%.dat.a,$(data))
analytic := $(patsubst analytic/%,%.alc.a,$(analytic))

arflags := rcs
arscript := ar.sh



.PHONY: plugclean

plugins.a: $(data) $(analytic)
+@echo "Building plugins library."
+@echo "CREATE $@" > $(arscript)
+@for a in $^; do (echo "ADDLIB $$a" >> $(arscript)); done
+@echo "SAVE" >> $(arscript)
+@echo "END" >> $(arscript)
+@ar -M < $(arscript)

%.dat.a: data/%
+@echo "Building data plugin $*"
+@cd $< && $(MAKE) TARGET=../../$@

%.alc.a: analytic/%
+@echo "Building analytic plugin $*"
+@cd $< && $(MAKE) TARGET=../../$@

plugclean:
+@echo "Cleaning all plugins."
+@rm -f *.a
