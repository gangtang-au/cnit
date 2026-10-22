Examples
========

CNit's emulation examples are rendered from the current Jupyter notebooks. They
show the complete workflows, code, and saved plots. Both notebooks use calibrated
parameters loaded from a model-specific CSV column and forcing data passed to
``CNitExpConfig`` objects.

CMIP6 ESM Emulation
-------------------

This notebook emulates 16 CMIP6 Earth system models across ``1pctCO2``,
``hist_ssp245``, and ``hist_ssp585``. It compares NPP, land carbon pools and
fluxes, and nitrogen pools and fluxes against the model outputs. To run it
locally, download ``cali_cmip6_best_fitting.csv``, ``CMIP6_cali_forcing.nc``,
and ``CMIP6_cali_target.nc`` from the `CNit Model Calibration and Results
<https://doi.org/10.5281/zenodo.18116573>`_ and put them in the project-level
``data/`` directory.

.. toctree::
   :maxdepth: 1

   CNit_example_CMIP6_ESMs_emulation

TRENDY DGVM Emulation
---------------------

This notebook emulates 20 TRENDY-v11 dynamic global vegetation models for the
historical ``S1``–``S3`` experiments, compares NPP and land carbon pools and
fluxes, and demonstrates an additional doubled-warming experiment. To run it
locally, download ``cali_trendy_best_fitting.csv``,
``TRENDY_cali_forcing.nc``, and ``TRENDY_cali_target.nc`` from the `CNit Model Calibration and Results <https://doi.org/10.5281/zenodo.18116573>`_ and put them in the
project-level ``data/`` directory.

.. toctree::
   :maxdepth: 1

   CNit_example_TRENDY_DGVMs_emulation
