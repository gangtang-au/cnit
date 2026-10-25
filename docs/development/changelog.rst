Changelog
=========

All notable changes to CNit will be documented here.

The format follows `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Version 2.2.0 (2026-10-22)
--------------------------

**Added**

* Model configuration helpers and support for loading calibration parameters
  from user-supplied CSV files.

**Changed**

* Refactored the carbon and nitrogen cycle calculations to use numeric values
  inside the ODE solver, while preserving unit-aware model inputs and outputs.
  The model equations and parameterizations are unchanged. The benchmarked
  simulation ran approximately 6 times faster than in CNit 2.1.0.
* Updated the quick start and examples to match the current model interface.

Version 2.1.0 (2026-01-01)
--------------------------

Initial public release of CNit as a Python package with enhanced code and full
documentation.

**Added**

* Python package structure with pip/conda installation and dependency management
* Unit handling with Pint integration for robust physical quantities
* Comprehensive Sphinx documentation with scientific background and API reference
* Quick start guide and complete example notebook with code, outputs, and visualizations

**Model Features**

* Coupled carbon-nitrogen cycle for terrestrial ecosystems
* Four nitrogen pools (plant, litter, soil, mineral) and three carbon pools
* CO₂ fertilization, temperature sensitivity, and land use change effects
* Nitrogen limitation on carbon uptake
* Annual timesteps with subannual process representation through within-year flux and
  pool-turnover allocation fractions. Evenly spaced subannual timesteps are also
  supported, with allocation fractions applying within each timestep.
* Separation of deforestation/afforestation for land use emissions
* Computationally efficient for integrated assessment modeling
* MAGICC integration ready

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
