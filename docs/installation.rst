Installation
============

This page provides instructions for installing KINC-Python on various platforms.

Requirements
------------

**Python Version:**

- Python 3.8 or later

**Core Dependencies:**

- NumPy >= 1.20.0
- SciPy >= 1.7.0
- pandas >= 1.3.0
- h5py >= 3.0.0
- scikit-learn >= 1.0.0
- tqdm >= 4.60.0

**Optional Dependencies:**

- **GPU support**: CuPy (for CUDA GPUs)
- **MPI support**: mpi4py (for distributed computing)
- **Parallel CPU**: Dask[distributed] (for multi-core parallelization)

Basic Installation
------------------

Using pip
~~~~~~~~~

The simplest way to install KINC-Python is via pip:

.. code-block:: bash

   pip install kinc-python

This will install KINC-Python and all required dependencies.

From Source
~~~~~~~~~~~

To install from source (for development or latest features):

.. code-block:: bash

   git clone https://github.com/SystemsGenetics/KINC-Python.git
   cd KINC-Python
   pip install -e .

The ``-e`` flag installs in "editable" mode, so changes to the code take effect immediately.

Installation with Optional Features
------------------------------------

GPU Support
~~~~~~~~~~~

For GPU-accelerated computation using NVIDIA GPUs:

.. code-block:: bash

   # Requires CUDA Toolkit installed on your system
   pip install kinc-python[gpu]

Or install CuPy separately:

.. code-block:: bash

   # For CUDA 12.x
   pip install cupy-cuda12x

   # For CUDA 11.x
   pip install cupy-cuda11x

**Note**: CuPy requires NVIDIA CUDA Toolkit to be installed. See the `CuPy installation guide <https://docs.cupy.dev/en/stable/install.html>`_ for details.

MPI Support
~~~~~~~~~~~

For distributed computing on HPC clusters:

.. code-block:: bash

   pip install kinc-python[mpi]

This requires an MPI implementation (OpenMPI, MPICH, etc.) to be installed on your system.

On Ubuntu/Debian:

.. code-block:: bash

   sudo apt install libopenmpi-dev openmpi-bin
   pip install kinc-python[mpi]

On macOS with Homebrew:

.. code-block:: bash

   brew install open-mpi
   pip install kinc-python[mpi]

Dask for Multi-Core Parallelization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For parallel processing on a single machine:

.. code-block:: bash

   pip install kinc-python[dask]

Full Installation
~~~~~~~~~~~~~~~~~

To install all optional dependencies:

.. code-block:: bash

   pip install kinc-python[all]

Development Installation
------------------------

For contributors and developers:

.. code-block:: bash

   git clone https://github.com/SystemsGenetics/KINC-Python.git
   cd KINC-Python
   pip install -e ".[dev]"

This installs additional development tools:

- pytest (testing)
- pytest-cov (code coverage)
- black (code formatting)
- flake8 (linting)

Platform-Specific Installation
-------------------------------

Ubuntu/Debian Linux
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Install system dependencies
   sudo apt update
   sudo apt install python3 python3-pip python3-dev
   sudo apt install libhdf5-dev  # For h5py

   # Install KINC-Python
   pip install kinc-python

For GPU support on Ubuntu:

.. code-block:: bash

   # Install NVIDIA drivers and CUDA
   sudo apt install nvidia-driver-535
   sudo apt install nvidia-cuda-toolkit

   # Install KINC-Python with GPU support
   pip install kinc-python[gpu]

macOS
~~~~~

.. code-block:: bash

   # Install Homebrew if not already installed
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

   # Install Python and HDF5
   brew install python hdf5

   # Install KINC-Python
   pip install kinc-python

Windows
~~~~~~~

.. code-block:: powershell

   # Install Python from python.org
   # Then use pip in PowerShell or Command Prompt

   pip install kinc-python

**Note**: GPU support on Windows requires NVIDIA CUDA Toolkit and appropriate CuPy wheel.

HPC Cluster Installation
-------------------------

Many HPC systems use module systems to manage software. Here's an example workflow:

.. code-block:: bash

   # Load required modules
   module load python/3.10
   module load gcc/11.2
   module load openmpi/4.1

   # Create a virtual environment
   python3 -m venv ~/kinc-env
   source ~/kinc-env/bin/activate

   # Install KINC-Python with MPI support
   pip install kinc-python[mpi]

For GPU nodes with CUDA:

.. code-block:: bash

   module load cuda/12.0
   pip install kinc-python[gpu,mpi]

Conda Installation
------------------

If you prefer conda/mamba:

.. code-block:: bash

   # Create conda environment
   conda create -n kinc python=3.10
   conda activate kinc

   # Install dependencies
   conda install numpy scipy pandas h5py scikit-learn

   # Install KINC-Python
   pip install kinc-python

Docker Installation
-------------------

A Docker image is available for easy deployment:

.. code-block:: bash

   docker pull systemsgenetics/kinc-python:latest

   # Run KINC in container
   docker run -v $(pwd):/data systemsgenetics/kinc-python kinc --help

Verifying Installation
-----------------------

To verify KINC-Python is installed correctly:

.. code-block:: bash

   # Check version
   kinc --version

   # Run help command
   kinc --help

   # Test import in Python
   python -c "import kinc; print(kinc.__version__)"

Troubleshooting
---------------

**Import Error: No module named 'kinc'**

Make sure you're in the correct Python environment:

.. code-block:: bash

   which python
   pip list | grep kinc

**HDF5 Installation Issues**

If h5py fails to install, you may need to install HDF5 system libraries:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt install libhdf5-dev

   # macOS
   brew install hdf5

   # Then reinstall h5py
   pip install --no-cache-dir h5py

**CuPy Installation Issues**

CuPy requires CUDA Toolkit. Verify CUDA is installed:

.. code-block:: bash

   nvcc --version
   nvidia-smi

See the `CuPy documentation <https://docs.cupy.dev/en/stable/install.html>`_ for platform-specific issues.

**MPI4Py Installation Issues**

Ensure MPI is installed:

.. code-block:: bash

   # Check for mpicc
   which mpicc

   # If not found, install OpenMPI
   sudo apt install libopenmpi-dev

Updating KINC-Python
--------------------

To update to the latest version:

.. code-block:: bash

   pip install --upgrade kinc-python

For development installations:

.. code-block:: bash

   cd KINC-Python
   git pull
   pip install -e .

Uninstallation
--------------

To remove KINC-Python:

.. code-block:: bash

   pip uninstall kinc-python
