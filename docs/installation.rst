Installation
============

Currently, there is no stand-alone binary for KINC. It must be compiled and installed. The instructions on this page provide step-by-step instructions to compile and install KINC.

KINC can also be run in a Docker container. Consult the Usage section for more information.

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

The instructions on this page provides details for compiling KINC.

Ubuntu 18.04
------------

Use the following steps to install KINC from source on Ubuntu 18.04:

Install Dependencies
~~~~~~~~~~~~~~~~~~~~

Most of these dependencies are available as packages on Ubuntu and can be installed using the `apt` framework:

.. code:: bash

   sudo apt install \
     qt5-default \
     libgsl-dev \
     libopenblas-dev \
     libopenmpi-dev \
     ocl-icd-opencl-dev \
     liblapacke-dev \
     nvidia-cuda-dev

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

KINC v3.3.0 requires ACE v3.1.0. ACE requires some of the same dependencies as KINC (such as QT, CUDA, OpenMPI, OpenCL, etc).  Therefore, if all dependencies above are installed, ACE should compile. To start, set the following environment variable:

.. code:: bash

   export ACE_VERSION=v3.1.0

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

   export ACE_VERSION=v3.1.0
   export KINC_VERSION=v3.3.0

Next, clone the KINC repository:

.. code:: bash

   git clone https://github.com/SystemsGenetics/KINC.git
   cd KINC/build
   git checkout $KINC_VERSION

Default installation location
*****************************

Next compile:

.. code:: bash

   qmake ../src/KINC.pro
   make qmake_all
   make
   make qmake_all
   make install

Alternative installation location
*********************************

By default, KINC will try to install itself into ``/usr/local``. To install KINC to a different directory (e.g. ``/local/software``), set the ``INSTALL_PREFIX`` environment variable accordingly:

.. code:: bash

   export INSTALL_PREFIX="/local/software"

Now build and install KINC:

   .. code:: bash

      qmake ../src/KINC.pro PREFIX=$INSTALL_PREFIX/KINC-$KINC_VERSION
      make qmake_all
      make
      make qmake_all
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

Windows
-------

Windows is currently not supported because there is no OpenMPI library for the Windows platform. Future support for Windows will be added when MPI becomes an optional dependency.

HPC Systems
-----------

Usage of KINC on high-performance computing (HPC) systems will require assistance of the cluster's systems admin to ensure all dependencies are installed and available.  Software management on clusters is specific to each cluster, although there are often commonalities.  Regardless, it is not possible to provide comprehensive instructions that would apply to every cluster.

Palmetto
~~~~~~~~

The following instructions are specific to the Palmetto cluster at Clemson University, however they can be adapted with some effort to other HPC clusters.

If you have previously used any version of KINC or ACE, be sure to remove the modules from your libraries. Furthermore, check to make sure that your ``.bashrc`` is clear of any designated paths for ACE or KINC.

Obtain an interactive node with at least 8 cores. Run the command:

.. code:: bash

   qsub -I -l select=1:ncpus=8

Once you have obtained an interactive node, run the following commands from your home directory:

.. code:: bash

   /zfs/feltus/btsheal/install-ace.sh
   /zfs/feltus/btsheal/install-kinc.sh

These scripts will install ACE and KINC into your home directory, establishing them as modules that can be run from anywhere. It will also update your environment so that the modules can be called when necessary. It uses a module called ``use.own``, which when added will make KINC and ACE available to be used interactively. You should now be able to load KINC as a module:

.. code:: bash

   module add use.own
   module add KINC
