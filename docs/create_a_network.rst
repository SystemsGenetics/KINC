Create a Network
================
To execute the functions described in the :doc:`functional_overview` section you will always

The easiest way to learn how to use KINC is to study the ``kinc.sh`` script, which can run the entire KINC workflow. Additionally, you can use the ``make-input-data.py`` script to generate a "fake" GEM with which to test KINC quickly:

.. code:: bash

   # generate fake GEM
   python scripts/make-input-data.py

   # run KINC
   scripts/kinc.sh serial 1 GEM.txt
