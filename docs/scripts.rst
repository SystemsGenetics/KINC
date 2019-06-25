Auxiliary Scripts
=================

Included in the KINC repository are several Bash and Python scripts for a variety of helpful tasks. These scripts are documented here. For each of these scripts, you can run the script with no arguments (or with ``-h`` for Python scripts) to see the help text.

kinc.sh
-------

The ``kinc.sh`` script can run each KINC analytic in a customizable sequence on any input GEM, using any runner (``serial``, ``cuda``, ``opencl``) and any number of processes via MPI. These three things are provided as command-line arguments. However, you can also set which analytics to run by setting the ``DO_*`` variables in the script source.

.. code:: bash

   scripts/kinc.sh serial 1 GEM.txt

kinc-py.sh
----------

The ``kinc-py.sh`` script does the same thing as ``kinc.sh`` but with Python scripts: ``kinc-similarity.py``, ``kinc-threshold.py``, and ``kinc-extract.py``. These scripts can also be used invidiually.

.. code:: bash

   scripts/kinc-py.sh GEM.txt
