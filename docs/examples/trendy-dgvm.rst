TRENDY DGVM Emulation
=====================

This workflow emulates 20 TRENDY-v11 dynamic global vegetation models for the
historical ``S1``–``S3`` experiments. It compares NPP and land carbon pools and
fluxes, and demonstrates an additional doubled-warming experiment.

The notebook source is ``notebooks/CNit_example_TRENDY_DGVMs_emulation.ipynb`` in
the CNit project repository. It is available in a source checkout, including when
CNit is installed in editable mode.

To run the notebook locally, download ``cali_trendy_best_fitting.csv``,
``TRENDY_cali_forcing.nc``, and ``TRENDY_cali_target.nc`` from the `CNit Model
Calibration and Results <https://doi.org/10.5281/zenodo.18116573>`_ and put them
in the project-level ``data/`` directory. See :doc:`../project_info/calibration-data`
for more about the archive.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   TRENDY notebook <../CNit_example_TRENDY_DGVMs_emulation>
