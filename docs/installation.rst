Installation
============

This page provides step-by-step instructions to download, compile and install KINC as well as its dependencies on a stand-alone Ubuntu 18.04 system. However, the same set of dependencies are required for installation on any other platform or on a High Performance Computing (HPC) platform.

Even if you install KINC on an HPC system, you may want to install KINC on a local stand-alone machine so that you can view output files using the `qkinc` graphical interface, export network files more easily or perform other tasks.

Additionally, KINC can also be run in a Docker container. Consult the :doc:`usage` section for more information.

.. _installation_reference_label:

Dependencies
------------

KINC requires the following software packages.

- `NVIDIA CUDA Toolkit <https://developer.nvidia.com/cuda-zone>`_
- `OpenCL <https://www.khronos.org/opencl/>`_
- `OpenMPI <https://www.open-mpi.org/>`_
- `QT <https://www.qt.io/>`_
- `ACE <https://github.com/SystemsGenetics/ACE>`_
- `The GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
- `OpenBLAS <https://www.openblas.net/>`_
- `LAPACK <http://www.netlib.org/lapack/>`_
- `GCEM <https://www.kthohr.com/gcem.html>`_
- `StatsLib <https://www.kthohr.com/statslib.html>`_
- `Boost C++ libraries <https://www.boost.org/>`_
Some functionality of KINC (i.e. condition-specific network construction) require two additional set of dependences:  `KINC.R v1.1 <https://github.com/SystemsGenetics/KINC.R>`_ , Python3 and a variety of Python modules for the 3D visualization.

To install KINC.R please follow the installation instructions on the `KINC.R repository <https://github.com/SystemsGenetics/KINC.R>`_. KINC.R requires a variety of other R modules.

If you desire to use the Python v3 `Plotly Dash <https://plotly.com/dash/>`_ 3D visualization script you must also install the following packages:

- `Numpy <https://numpy.org/>`_
- `Pandas <https://pandas.pydata.org/>`_
- `IGraph for Python <https://igraph.org/python/>`_
- `Plotly <https://plotly.com/>`_
- `Seaborn <https://seaborn.pydata.org/>`_
- `Dash <https://plotly.com/dash/>`_
- `progress <https://github.com/verigak/progress/>`_
- `fa2 <https://github.com/bhargavchippada/forceatlas2>`_

Install these Python v3 packages using your favorite package manager (i.e. pip, Anaconda, etc.)

Installing KINC on Ubuntu 18.04
-------------------------------

Install Dependencies
~~~~~~~~~~~~~~~~~~~~

Most dependencies are available as packages via Ubuntu and can be installed using the `apt` framework:

.. code:: bash

  sudo apt install \
    qt5-default \
    libgsl-dev \
    libopenblas-dev \
    libopenmpi-dev \
    ocl-icd-opencl-dev \
    liblapacke-dev \
    nvidia-cuda-toolkit \
    libboost-dev

Additionally, since KINC uses the NVIDIA Driver API, you must install either the appropriate NVIDIA drivers for your system or the NVIDIA headless driver if you don't have a GPU:

.. code:: bash

  # install NVIDIA driver
  sudo apt install nvidia-driver-435

  # install NVIDIA headless driver
  sudo apt install nvidia-headless-435

For specific device drivers other than those provided by Ubuntu (i.e. AMD, Intel, NVIDIA, etc), please refer to the manufacturer's website for installation instructions.

Install StatsLib and GCEM
~~~~~~~~~~~~~~~~~~~~~~~~~

Both StatsLib and GCEM are header-only libraries. To install them, you only need to download the packages and put them where they can be found.  The easiest location is in ``/usr/local/include`` (which requires root access).  Below are example installation instructions:

To install StatsLib into ``/usr/local/``:

.. code:: bash

  git clone -b master --single-branch https://github.com/kthohr/stats ./stats
  sudo cp -R ./stats/include/* /usr/local/include


To install CGEM into ``/usr/local/``:

.. code:: bash

  git clone https://github.com/kthohr/gcem.git ./gcem
  sudo cp -R ./gcem/include/* /usr/local/include

Install ACE
~~~~~~~~~~~

KINC v3.4 requires ACE v3.2. ACE requires some of the same dependencies as KINC (such as QT, CUDA, OpenMPI, OpenCL, etc).  Therefore, if all dependencies above are installed, ACE should compile. To start, set the following environment variable:

.. code:: bash

  export ACE_VERSION=v3.2.0

Next, clone the ACE repository:

.. code:: bash

  git clone https://github.com/SystemsGenetics/ACE.git
  cd ACE/build
  git checkout $ACE_VERSION

Default installation location
*****************************

Next compile:

.. code:: bash

  qmake ../src/ACE.pro
  make qmake_all
  make
  make qmake_all
  make install

Alternative installation location
*********************************

By default, ACE will try to install into ``/usr/local``. To install ACE to a different directory (e.g. ``/local/software``), set the ``INSTALL_PREFIX`` environment variable accordingly:

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
~~~~~~~~~~~~

Select a suitable `version of KINC <https://github.com/SystemsGenetics/KINC/releases>`__ and set the environment variable:

.. code:: bash

  export ACE_VERSION=v3.2.0
  export KINC_VERSION=v3.4.2

Next, clone the KINC repository:

.. code:: bash

  git clone https://github.com/SystemsGenetics/KINC.git
  cd KINC
  git checkout $KINC_VERSION

Default installation location
*****************************

Next compile:

.. code:: bash

  make
  make install

Alternative installation location
*********************************

By default, KINC will try to install itself into ``/usr/local``. To install KINC to a different directory (e.g. ``/local/software``), set the ``INSTALL_PREFIX`` environment variable accordingly:

.. code:: bash

  export INSTALL_PREFIX="/local/software"

Now build and install KINC:

.. code:: bash

  make
  make install

If ACE is not in /usr/local
***************************

If ACE was not installed into an alternative location other than the default ``/usr/local`` then should set several environment variables help the compiler find ACE libraries and headers:

.. code:: bash

  export PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/bin:$PATH"
  export LD_LIBRARY_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$LD_LIBRARY_PATH"
  export LIBRARY_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$LIBRARY_PATH"
  export CPATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$CPATH"
  export C_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$C_INCLUDE_PATH"
  export CPLUS_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$CPLUS_INCLUDE_PATH"
  export OBJC_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$OBJC_INCLUDE_PATH"


Preparing to Run KINC
~~~~~~~~~~~~~~~~~~~~~

If KINC was installed in the default location you can skip the :doc:`usage` page for futher instructions, otherwise, if you installed KINC in an alternative location, you must update the ``LD_LIBRARY_PATH`` in your ``~/.bashrc`` file.  Use the following command to get the exact text you need to add.

.. code:: bash

  echo "export LD_LIBRARY_PATH=\"$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$INSTALL_PREFIX/KINC-$KINC_VERSION/lib:\$LD_LIBRARY_PATH\""
  echo "export PATH=\"$INSTALL_PREFIX/ACE-$ACE_VERSION/bin:$INSTALL_PREFIX/KINC-$KINC_VERSION/bin:\$PATH\""

Append the resulting text to your ``~/.bashrc`` file. You should now be able to run KINC


Installing on Windows
---------------------

Windows is currently not supported because there is no OpenMPI library for the Windows platform. Future support for Windows will be added when MPI becomes an optional dependency.

Installing on an HPC System
---------------------------

Usage of KINC on high-performance computing (HPC) systems will require assistance of the cluster's systems admin to ensure all dependencies are installed and available.  Software management on clusters is specific to each cluster, although there are often commonalities.  Regardless, it is not possible to provide comprehensive instructions that would apply to every cluster.

Palmetto
~~~~~~~~

The following instructions are specific to the Palmetto cluster at Clemson University, however they can be adapted with some effort to other HPC clusters.

If you have previously used any version of KINC or ACE, be sure to remove the modules from your libraries. Furthermore, check to make sure that your ``.bashrc`` is clear of any designated paths for ACE or KINC.

Obtain an interactive node with at least 8 cores. Run the command:

.. code:: bash

  qsub -I -l select=1:ncpus=8:ngpus=2:gpu_model=p100

Once you have obtained an interactive node, run the following commands from your home directory:

.. code:: bash

  git clone https://github.com/bentsherman/pbs-toolkit.git
  ./pbs-toolkit/modules/install-ace.sh v3.2.0
  ./pbs-toolkit/modules/install-statslib.sh
  ./pbs-toolkit/modules/install-kinc.sh v3.4.2 v3.2.0

These scripts will install ACE and KINC into your home directory, establishing them as modules that can be run from anywhere. It will also update your environment so that the modules can be called when necessary. It uses a module called ``use.own``, which when added will make KINC and ACE available to be used interactively. You should now be able to load KINC as a module:

.. code:: bash

  module use ${HOME}/modules
  module load kinc/v3.4.2
