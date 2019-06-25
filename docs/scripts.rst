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

make-input-data.py
------------------

The ``make-input-data.py`` script creates a "fake" GEM given the number of genes, samples, and classes, as well as other customization options. This script is useful for testing KINC when you don't have any real data on hand.

.. code:: bash

   python scripts/make-input-data.py

validate.py
-----------

The ``validate.py`` script attempts to measure the difference between two similarity matrices by comparing the features of each edge (number of clusters, sample statistics, correlation value, sample mask, etc.). This script is useful for comparing similarity matrices from different runs to confirm that they are equivalent or to simply quantify their differences. You must provide the plain-text similarity matrix file, not the network file. This file can be generated using the ``export-cmx`` analytic.

.. code:: bash

   python scripts/validate.py GEM1.cmx.txt GEM2.cmx.txt

visualize.py
------------

The ``visualize.py`` script can create several useful visualizations of a network, such as scatter plots of the gene pairs in the network. Consult the help text to see all visualization options.

.. code:: bash

   python scripts/visualize.py --emx GEM.txt --netlist GEM.coexpnet.txt --corrdist
