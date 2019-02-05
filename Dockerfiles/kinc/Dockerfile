FROM nvidia/cuda:9.0-devel
MAINTAINER Ben Shealy <btsheal@clemson.edu>

# install package dependencies
RUN apt-get update -qq \
	&& apt-get install -qq -y \
		wget xz-utils \
		build-essential pkg-config libgl1-mesa-dev \
		clinfo git libgsl-dev liblapacke-dev libopenblas-dev libopenmpi-dev ocl-icd-opencl-dev \
		python3-pip

# add NVIDIA platform to OpenCL
RUN mkdir -p /etc/OpenCL/vendors \
	&& echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# install Qt
WORKDIR /opt

ENV QTDIR "/opt/Qt-5.11.1"
ENV PATH "$QTDIR/bin:$PATH"

RUN wget -q http://download.qt.io/official_releases/qt/5.11/5.11.1/single/qt-everywhere-src-5.11.1.tar.xz \
	&& tar -xf qt-everywhere-src-5.11.1.tar.xz \
	&& rm qt-everywhere-src-5.11.1.tar.xz \
	&& cd qt-everywhere-src-5.11.1 \
	&& ./configure -opensource -confirm-license -prefix $QTDIR \
		-skip qt3d \
		-skip qtactiveqt \
		-skip qtandroidextras \
		-skip qtcanvas3d \
		-skip qtcharts \
		-skip qtconnectivity \
		-skip qtdatavis3d \
		-skip qtdeclarative \
		-skip qtdoc \
		-skip qtgamepad \
		-skip qtgraphicaleffects \
		-skip qtimageformats \
		-skip qtlocation \
		-skip qtmacextras \
		-skip qtmultimedia \
		-skip qtnetworkauth \
		-skip qtpurchasing \
		-skip qtquickcontrols \
		-skip qtquickcontrols2 \
		-skip qtremoteobjects \
		-skip qtscript \
		-skip qtscxml \
		-skip qtsensors \
		-skip qtserialbus \
		-skip qtserialport \
		-skip qtspeech \
		-skip qtsvg \
		-skip qttools \
		-skip qttranslations \
		-skip qtvirtualkeyboard \
		-skip qtwayland \
		-skip qtwebchannel \
		-skip qtwebengine \
		-skip qtwebglplugin \
		-skip qtwebsockets \
		-skip qtwebview \
		-skip qtwinextras \
		-skip qtx11extras \
		-skip qtxmlpatterns \
	&& make -s -j $(nproc) \
	&& make -s install \
	&& cd .. \
	&& rm -rf qt-everywhere-src-5.11.1

# install ACE
WORKDIR /opt

RUN git clone https://github.com/SystemsGenetics/ACE.git \
	&& cd ACE/build \
	&& git checkout v3.0.2 \
	&& qmake ../src/ACE.pro \
	&& make -s -j $(nproc) \
	&& make -s qmake_all \
	&& make -s install

ENV LD_LIBRARY_PATH "/usr/local/lib:$LD_LIBRARY_PATH"

# install KINC
WORKDIR /opt

RUN git clone https://github.com/SystemsGenetics/KINC.git \
	&& cd KINC/build \
	&& git checkout v3.2.2 \
	&& qmake ../src/KINC.pro \
	&& make -s -j $(nproc) \
	&& make -s qmake_all \
	&& make -s install

# install python dependencies
RUN pip3 install -q argparse matplotlib numpy pandas scipy seaborn

# return to home directory
WORKDIR /root
