FROM nvidia/cuda:10.0-devel
MAINTAINER Ben Shealy <btsheal@clemson.edu>

ARG NVIDIA_HEADLESS=0
ARG ACE_REVISION="develop"
ARG KINC_REVISION="master"

# install package dependencies
RUN apt-get update -qq \
	&& apt-get install -qq -y \
		qt5-default \
		clinfo git libgsl-dev liblapacke-dev libopenblas-dev libopenmpi-dev ocl-icd-opencl-dev \
		python3-pip

# install headless driver for cpu image
RUN if [ ${NVIDIA_HEADLESS} = 1 ]; then apt-get install -qq -y nvidia-headless-418 ; fi

# add NVIDIA platform to OpenCL
RUN mkdir -p /etc/OpenCL/vendors \
	&& echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# install python dependencies
RUN pip3 install -q argparse matplotlib numpy pandas scipy seaborn

# install StatsLib
WORKDIR /opt

RUN git clone -q https://github.com/kthohr/gcem \
	&& cp -r gcem/include/* /usr/local/include

RUN git clone -q https://github.com/kthohr/stats \
	&& cp -r stats/include/* /usr/local/include

# install ACE
WORKDIR /opt

RUN git clone -q https://github.com/SystemsGenetics/ACE.git \
	&& cd ACE/build \
	&& git checkout -q ${ACE_REVISION} \
	&& qmake ../src/ACE.pro GUI=no \
	&& make -s -j $(nproc) \
	&& make -s qmake_all \
	&& make -s install

ENV LD_LIBRARY_PATH "/usr/local/lib:$LD_LIBRARY_PATH"

# install KINC
WORKDIR /opt

RUN git clone -q https://github.com/SystemsGenetics/KINC.git \
	&& cd KINC/build \
	&& git checkout -q ${KINC_REVISION} \
	&& qmake ../src/KINC.pro GUI=no \
	&& make -s -j $(nproc) \
	&& make -s qmake_all \
	&& make -s install

# initialize default work directory
WORKDIR /workspace
