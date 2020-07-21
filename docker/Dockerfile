FROM nvidia/cuda:10.0-devel
MAINTAINER Ben Shealy <btsheal@clemson.edu>

ARG NVIDIA_HEADLESS=0
ARG ACE_REVISION="develop"
ARG KINC_REVISION="develop"
ARG KINC_R_REVISION="master"

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


# Install headless driver for cpu image
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

RUN git clone -q https://github.com/SystemsGenetics/ACE.git \
	&& cd ACE/build \
	&& git checkout -q ${ACE_REVISION} \
	&& qmake ../src/ACE.pro \
	&& make -s -j $(nproc) \
	&& make -s qmake_all \
	&& make -s install

ENV LD_LIBRARY_PATH "/usr/local/lib:$LD_LIBRARY_PATH"

# install KINC
WORKDIR /opt

# For older version of KINC we need to use qmake. For newer versions after 3.3.x we
# can use the make file in the root directory of KINC.
RUN if echo "${KINC_REVISION}" | grep -Eq 'v3\.[32]\.[0-9]$' ; then \
    git clone -q https://github.com/SystemsGenetics/KINC.git \
    && cd KINC \
	&& git checkout -q ${KINC_REVISION} \
    && cd build \
	&& qmake ../src/KINC.pro \
	&& make -s -j $(nproc) \
	&& make -s qmake_all \
	&& make -s install; \
else \
    git clone -q https://github.com/SystemsGenetics/KINC.git \
    && cd KINC \
	&& git checkout -q ${KINC_REVISION} \
	&& make -s -j $(nproc) \
	&& make -s install; \
fi


# Add in a few additional requirements for the >=3.4 version of KINC. These are needed
# for the R and Python scripts in the bin folder.
RUN if ! echo "${KINC_REVISION}" | grep -Eq 'v3\.[32]\.[0-9]$' ; then \
    apt-get install -qq -y \
        m4 \
        autoconf \
        automake \
        libxml2-dev \
        bison \
        flex \
        build-essential \
        libcurl4-gnutls-dev \
        libssl-dev \
        software-properties-common; \
fi

# Install KINC.R
ENV DEBIAN_FRONTEND=noninteractive
RUN if ! echo "${KINC_REVISION}" | grep -Eq 'v3\.[32]\.[0-9]$' ; then \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
    && apt-get update -qq \
    && apt-get -qq install -y r-base \
    && R -q -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.us.r-project.org')" \
    && R -e "library('devtools'); install_github('SystemsGenetics/KINC.R', ref='${KINC_R_REVISION}')"; \
fi

# install python dependencies
WORKDIR /opt/KINC

RUN if [ -e requirements.txt ]; then pip3 install  -r requirements.txt; fi



# initialize default work directory
WORKDIR /workspace

# Tini for signal processing and zombie killing.
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini


# Define the command and parameters that will be executed when this
# container is first run.
ENTRYPOINT ["/tini", "--"]
