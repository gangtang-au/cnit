Changelog
=========

All notable changes to CNit will be documented here.

The format follows `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Version 2.4.0 (2026-10-31)
--------------------------

**Changed**

* Replaced the mineral nitrogen pool's ODE solve with an analytical interval update
  using linearly interpolated net input and a harmonic-mean effective turnover time.
* Improved interpolation performance by packing synchronized numeric series into a
  shared linear interpolator.
* Significantly improved model computational performance. For example, running the
  TRENDY S0–S3 experiments (1700–2020) for 20 DGVMs now takes approximately 1 second,
  while running the 1pctCO2 (1850–2000), hist_ssp245 (1850–2100), and hist_ssp585
  (1850–2100) experiments for 16 CMIP6 ESMs takes approximately 1.2 seconds (see the
  example notebooks in the documentation).
* ``run_experiments()`` now applies each experiment's nitrogen switch to a copy of the
  model, leaving the original model unchanged.
* Pool switches are now validated as fixed-length binary sequences and stored as
  immutable tuples.
* Model parameters and experiment input series are checked for finite values during
  configuration validation.

Version 2.3.0 (2026-10-25)
--------------------------

**Changed**

* Optimized regrowth and afforestation-decay pulse scheduling by recording rate
  changes at interval boundaries instead of updating every affected future interval.
  This preserves the pulse rates for a given duration while reducing repeated work.
* Improved subannual time-axis support and documentation across CNit. In land-use
  flux calculations, regrowth and afforestation-decay durations are rounded to the
  nearest whole model interval (minimum one), rather than truncated with ``int()`` as
  before.
* Strengthened validation of model and experiment configurations and calibration CSV
  files.
* Reorganized and expanded the documentation, including the configuration API.

Version 2.2.0 (2026-10-22)
--------------------------

**Added**

* Model configuration constructors from dictionaries and calibration CSV files.
* Quickstart, CMIP6 ESM, and TRENDY DGVM notebooks, together with example calibration
  tables and forcing and target data.

**Changed**

* Refactored the carbon and nitrogen cycle calculations to use numeric values
  inside the ODE solver, while preserving unit-aware model inputs and outputs.
  The model equations and parameterizations are unchanged. The benchmarked
  simulation ran approximately 6 times faster than in CNit v2.1.0.
* Updated the quickstart and emulation examples for the CNit v2.2.0 model interface.

Version 2.1.0 (2026-01-01)
--------------------------

First public release of CNit as a Python package, with code improvements and
comprehensive documentation.

**Added**

* Python package structure with pip and conda installation and dependency management.
* Pint-based unit handling for physical quantities.
* Sphinx documentation with scientific background and API reference.
* Quickstart guide and a TRENDY DGVM emulation workflow with code, outputs, and
  visualizations.

**Model Features**

* Coupled carbon-nitrogen cycle for terrestrial ecosystems
* Four nitrogen pools (plant, litter, soil, mineral) and three carbon pools
* CO₂ fertilization, temperature sensitivity, and land use change effects
* Nitrogen limitation on carbon uptake
* Annual timesteps with subannual process representation through within-year flux and
  pool-turnover allocation fractions. Evenly spaced subannual timesteps are also
  supported, with allocation fractions applying within each timestep.
* Separation of deforestation and afforestation for land-use emissions.
* Computational efficiency for integrated assessment modeling.
* Readiness for MAGICC integration.

See the :doc:`model overview <../overview>` for a detailed model description.

Earlier Source Releases
-----------------------

Source code for earlier CNit versions is archived in these Zenodo releases:

* Tang, G., Nicholls, Z., Norton, A., Zaehle, S., & Meinshausen, M. (2024).
  *Model Source Code for "Coupled Carbon-Nitrogen Cycle in MAGICC: Model Description and Calibration"*.
  Zenodo. https://doi.org/10.5281/zenodo.12204421

* Tang, G., Zaehle, S., Nicholls, Z., Norton, A., Ziehn, T., & Meinshausen, M. (2025).
  *Model Source Code for the MAGICC Coupled Carbon–Nitrogen Cycle Model – CNit*.
  Zenodo. https://doi.org/10.5281/zenodo.15569386

.. note::

   Version numbering follows `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_
   (MAJOR.MINOR.PATCH):

   * **MAJOR:** Incompatible API changes
   * **MINOR:** New features, backward-compatible
   * **PATCH:** Bug fixes, backward-compatible
