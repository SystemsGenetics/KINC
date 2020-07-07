FROM nvidia/cuda:10.0-devel
MAINTAINER Ben Shealy <btsheal@clemson.edu>

ARG NVIDIA_HEADLESS=0
ARG ACE_REVISION="develop"
ARG KINC_REVISION="develop"

ENV CUDADIR="/usr/local/cuda"
ENV CPLUS_INCLUDE_PATH="${CUDADIR}/include:${CPLUS_INCLUDE_PATH}"
ENV LIBRARY_PATH="${CUDADIR}/lib64:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${CUDADIR}/lib64:${LD_LIBRARY_PATH}"

# install package dependencies
RUN apt-get update -qq \
	&& apt-get install -qq -y \
		qt5-default \
		clinfo \
		git \
		libboost-dev \
		libgsl-dev \
		liblapacke-dev \
		libopenblas-dev \
		libopenmpi-dev \
		ocl-icd-opencl-dev \
		python3-pip

# install headless driver for cpu image
RUN if [ ${NVIDIA_HEADLESS} = 1 ]; then apt-get install -qq -y nvidia-headless-418 ; fi

# add NVIDIA platform to OpenCL
RUN mkdir -p /etc/OpenCL/vendors \
	&& echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# install StatsLib
WORKDIR /opt

RUN git clone -q https://github.com/kthohr/gcem \
	&& cp -r gcem/include/* /usr/local/include

RUN git clone -q https://github.com/kthohr/stats \
	&& cp -r stats/include/* /usr/local/include

# install ACE
WORKDIR /opt

ADD https://api.github.com/repos/SystemsGenetics/ACE/git/refs/heads/${ACE_REVISION} ace-version.json
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

ADD https://api.github.com/repos/SystemsGenetics/KINC/git/refs/heads/${KINC_REVISION} kinc-version.json
RUN git clone -q https://github.com/SystemsGenetics/KINC.git \
	&& mkdir -p KINC/build \
	&& cd KINC/build \
	&& git checkout -q ${KINC_REVISION} \
	&& qmake ../src/KINC.pro GUI=no \
	&& make -s -j $(nproc) \
	&& make -s qmake_all \
	&& make -s install

# install python dependencies
WORKDIR /opt/KINC

RUN pip3 install -q -r requirements.txt

# initialize default work directory
WORKDIR /workspace
