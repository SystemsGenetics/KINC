
DEBUG ?= 1

NVCC = nvcc
NVCCFLAGS = -std=c++11

ifeq ($(DEBUG), 1)
NVCCFLAGS += --device-debug
else
NVCCFLAGS += --generate-line-info
endif

kinc-mini: kinc-mini.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -f kinc-mini
