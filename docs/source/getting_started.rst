Getting started
=================
This repository provides importable functions that can be used to
assemble a pipeline as needed. An example is in the helper script
template_pipeline.py, which you can follow in the next
pages.

In the config.ini file, replace the DIR_HOME string with the path to your copy
of the repository. (The other directories will be made by the pipeline.)

You can also specify other directories in the config.ini file if you want to change them
from the defaults.

Run

.. code-block:: bash

  python template_pipeline.py

which will proceed to make directories with the function make_dirs(). Put the pinhole
image into its directory, and rerun the above command.

Making the directories: comparison.angOffset_plateScale()
^^^^^^^^^

This function will make CDF plots with 1-sigma-equivalent boundaries,
like these:

.. _label: This is an example picture
.. figure:: images/plate_scale_sx_2018dec.pdf
	   :scale: 70 %
           :align: center
	   :alt: Alternative text
