CMIP6 ESM Emulation
===================

This workflow emulates 16 CMIP6 Earth system models across ``1pctCO2``,
``hist_ssp245``, and ``hist_ssp585``. It compares NPP, land carbon pools and
fluxes, and nitrogen pools and fluxes against the model outputs.

The notebook source is ``notebooks/CNit_example_CMIP6_ESMs_emulation.ipynb`` in
the CNit project repository. It is available in a source checkout, including when
CNit is installed in editable mode.

To run the notebook locally, download ``cali_cmip6_best_fitting.csv``,
``CMIP6_cali_forcing.nc``, and ``CMIP6_cali_target.nc`` from the `CNit Model
Calibration and Results <https://doi.org/10.5281/zenodo.18116573>`_ and put them
in the project-level ``data/`` directory. See :doc:`../project_info/calibration-data`
for more about the archive.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   CMIP6 notebook <../CNit_example_CMIP6_ESMs_emulation>
