CMIP6 ESM Emulation
===================

This workflow emulates 16 CMIP6 Earth system models across ``1pctCO2``,
``hist_ssp245``, and ``hist_ssp585``. It compares NPP, land carbon pools and
fluxes, and nitrogen pools and fluxes against the model outputs.

The `CMIP6 ESM emulation notebook
<https://github.com/gangtang-au/cnit/blob/main/notebooks/CNit_example_CMIP6_ESMs_emulation.ipynb>`_
is in the project-level ``notebooks/`` directory. It is available in a source
checkout, including when CNit is installed in editable mode.

To run the notebook locally, download ``cali_cmip6_best_fitting.csv``,
``CMIP6_cali_forcing.nc``, and ``CMIP6_cali_target.nc`` from the `CNit Model
Calibration and Results <https://doi.org/10.5281/zenodo.18116573>`_ and put them
in the project-level ``data/`` directory. See :doc:`../project_info/calibration-data`
for more about the archive.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   CMIP6 notebook <../CNit_example_CMIP6_ESMs_emulation>
