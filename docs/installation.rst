.. _installation:

Installation
------------

Ubuntu
~~~~~~

Use the following steps to install KINC from source on Ubuntu 16.04:

Install Dependencies
====================

Most of the dependencies are available as packages:

.. code:: bash

   sudo apt install build-essential libgsl-dev libopenblas-dev libopenmpi-dev ocl-icd-opencl-dev liblapacke-dev

For device drivers (AMD, Intel, NVIDIA, etc), refer to the manufacturer's website.

Install Qt (>=5.7)
==================

Select a suitable `version of Qt <http://download.qt.io/official_releases/qt>`__ and install Qt:

.. code:: bash

   wget http://download.qt.io/official_releases/qt/5.7/5.7.1/qt-opensource-linux-x64-5.7.1.run
   sh ./qt-opensource-linux-x64-5.7.1.run

If you install Qt locally then you must add Qt to the executable path:

.. code:: bash

   # append to ~/.bashrc
   export QTDIR="$HOME/Qt/5.7.1/gcc_64"
   export PATH="$QTDIR/bin:$PATH"

Install ACE
===========

Select a suitable `version of ACE <https://github.com/SystemsGenetics/ACE/releases>`__ and set the following environment variable:

.. code:: bash

   export ACE_VERSION=v3.0.2

Next, clone the ACE repository:

.. code:: bash

   git clone https://github.com/SystemsGenetics/ACE.git
   cd ACE/build
   git checkout $ACE_VERSION

By default, ACE will try to install itself into ``/usr/local``. To install ACE to a different directory (e.g. ``/local/software``), set the ``INSTALL_PREFIX`` environment variable accordingly:

.. code:: bash

   export INSTALL_PREFIX="/local/software"

Now, within the ``ACE/build`` directory run the following to build the ACE libraries:

.. code:: bash

   qmake ../src/ACE.pro PREFIX=$INSTALL_PREFIX/ACE-$ACE_VERSION
   make qmake_all
   make
   make qmake_all
   make install

This will install ACE into the directory specified by ``INSTALL_PREFIX`` in a directory named with the ACE version.

Install KINC
============

Select a suitable `version of KINC <https://github.com/SystemsGenetics/KINC/releases>`__ and set the environment variable:

.. code:: bash

   export ACE_VERSION=v3.0.2
   export KINC_VERSION=v3.2.2

Next, clone the KINC repository:

.. code:: bash

   git clone https://github.com/SystemsGenetics/KINC.git
   cd KINC/build
   git checkout $KINC_VERSION

By default, KINC will try to install itself into ``/usr/local``. To install KINC to a different directory (e.g. ``/local/software``), set the ``INSTALL_PREFIX`` environment variable accordingly:

.. code:: bash

   export INSTALL_PREFIX="/local/software"

Before you can build KINC, the compiler must be able to find the ACE libraries.  Several environment variables help with this:

.. code:: bash

   export PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/bin:$PATH"
   export LD_LIBRARY_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$LD_LIBRARY_PATH"
   export LIBRARY_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$LIBRARY_PATH"
   export CPATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$CPATH"
   export C_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$C_INCLUDE_PATH"
   export CPLUS_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$CPLUS_INCLUDE_PATH"
   export OBJC_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$OBJC_INCLUDE_PATH"

Now build and install KINC:

.. code:: bash

   qmake ../src/KINC.pro PREFIX=$INSTALL_PREFIX/KINC-$KINC_VERSION
   make qmake_all
   make
   make qmake_all
   make install

To run KINC you must update the ``LD_LIBRARY_PATH`` in your ``~/.bashrc`` file.  Use the following command to get the exact text you need to add.

.. code:: bash

   echo "export LD_LIBRARY_PATH=\"$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$INSTALL_PREFIX/KINC-$KINC_VERSION/lib:\$LD_LIBRARY_PATH\""
   echo "export PATH=\"$INSTALL_PREFIX/ACE-$ACE_VERSION/bin:$INSTALL_PREFIX/KINC-$KINC_VERSION/bin:\$PATH\""

Append the resulting text to your ~/.bashrc file. You should now be able to run KINC

Windows
~~~~~~~

Windows is currently not supported because there is no OpenMPI library for the Windows platform. Future support for Windows will be added when MPI becomes an optional dependency.
