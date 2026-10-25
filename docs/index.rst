CNit Documentation
==================

Welcome to the documentation for **CNit** (Carbon-Nitrogen Interactions in Terrestrial
ecosystems), the coupled carbon-nitrogen cycle model for **MAGICC** (Model for the
Assessment of Greenhouse Gas Induced Climate Change).

.. image:: _static/CNit_logo.png
   :alt: CNit logo
   :height: 150px
   :align: right

CNit is a process-based terrestrial biogeochemistry model that simulates the coupled
dynamics of carbon and nitrogen cycles in terrestrial ecosystems. The model tracks carbon
and nitrogen through plant, litter, soil, and mineral pools, with explicit representation
of carbon-nitrogen feedbacks, environmental effects (CO₂ fertilization, temperature, land
use change), and nitrogen limitation on carbon uptake.

**Key Features:**

* Explicit carbon-nitrogen coupling with nitrogen limitation feedbacks
* Four nitrogen pools (plant, litter, soil, mineral) and three carbon pools
* Environmental modifiers: CO₂, temperature, land use, nitrogen availability
* Separation of deforestation and afforestation for CDR scenarios
* Annual timesteps use flux-allocation and pool-turnover fractions to represent
  within-year transfers
* Evenly spaced subannual timesteps are also supported; allocation fractions apply
  within each timestep
* Computationally efficient for integrated assessment modeling

Quick Navigation
----------------

**🔬 For Scientists**
   Understand the model's scientific basis, processes, and when to use CNit.

   → Start with :doc:`overview`, then see published papers in References.

**📚 For Model Users**
   Install CNit and run your first simulation.

   → See :doc:`quickstart` for installation and basic usage.

   → See :doc:`api/cnit.configs` for model parameters and experiment inputs.

   → See :doc:`example` for complete CMIP6 ESM and TRENDY DGVM notebooks with code and outputs.

**🛠️ For Developers**
   Contribute to CNit development.

   → See :doc:`api/cnit.physics` for implementation details with scientific descriptions.

   → See :doc:`contributing` for development setup and coding guidelines.

Getting Started
---------------

New to CNit? Start here to understand the model and run your first simulation.

.. toctree::
   :caption: Getting started
   :maxdepth: 2

   quickstart

API Reference
-------------

The overview introduces CNit's structure and processes. These API pages provide the
detailed process parameterizations, model calculations, configuration parameters, and
class methods.

.. toctree::
   :caption: API reference
   :maxdepth: 2

   api/cnit.physics
   api/cnit.configs

Examples
--------

See complete CMIP6 ESM and TRENDY DGVM emulation workflows.

.. toctree::
   :caption: Examples
   :maxdepth: 2

   example

Model Overview
--------------

Read about CNit's structure, processes, and scientific basis.

.. toctree::
   :caption: Model overview
   :maxdepth: 2

   overview

Development
-----------

.. toctree::
   :caption: Development
   :maxdepth: 2

   changelog
   contributing


Contact
-------

CNit is developed and maintained by **Gang Tang**.

**Questions or feedback?**

* Open an issue: `GitHub Issues <https://github.com/gangtang-au/cnit/issues>`_
* Email: gang.tang.au@gmail.com

Interested in contributing or collaboration? Please get in touch!

License & Citation
------------------

CNit is released under the **BSD 3-Clause License**.

Copyright (c) 2026, Gang Tang and contributors.

See the `LICENSE file <https://github.com/gangtang-au/cnit/blob/main/LICENSE>`_
for full details.

**If you use CNit in your research, please cite our papers.**

* Tang, G., Nicholls, Z., Norton, A., Zaehle, S., and Meinshausen, M.: Synthesizing
  global carbon–nitrogen coupling effects – the MAGICC coupled carbon–nitrogen
  cycle model v1.0, Geosci. Model Dev., 18, 2193–2230,
  https://doi.org/10.5194/gmd-18-2193-2025, 2025.

* Tang, G., Zaehle, S., Nicholls, Z., Norton, A., Ziehn, T., & Meinshausen, M. (2026).
  Understanding the drivers of carbon–nitrogen cycle variability in CMIP6 ESMs
  with MAGICC CNit v2.0: Model and calibration updates. Journal of Advances in
  Modeling Earth Systems, 18, e2025MS005270. https://doi.org/10.1029/2025MS005270

Related Projects
----------------

MAGICC - Model for the Assessment of Greenhouse Gas Induced Climate Change

* https://magicc.org/
* https://gitlab.com/magicc/magicc

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
