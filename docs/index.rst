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
* Environmental modifiers:  CO₂, temperature, land use, nitrogen availability
* Separation of deforestation and afforestation for CDR scenarios
* Annual timestep with sub-annual process representation
* Computationally efficient for integrated assessment modeling

Quick Navigation
----------------

**🔬 For Scientists**
   Understand the model's scientific basis, processes, and when to use CNit.

   → Start with :doc:`overview`, then see published papers in References.

**📚 For Model Users**
   Install CNit and run your first simulation.

   → See :doc:`quickstart` for installation and basic usage.

   → See :doc:`example` for a complete example with code and outputs.

**🛠️ For Developers**
   Contribute to CNit development.

   → See :doc:`api/cnit.physics` for implementation details with scientific descriptions.

   → See :doc:`contributing` for development setup and coding guidelines.

Understanding CNit
------------------

Learn about the model's structure, processes, and scientific basis before installing.

.. toctree::
   :maxdepth: 2

   overview

Getting Started
---------------

New to CNit? Start here to understand the model and run your first simulation.

.. toctree::
   :maxdepth: 2

   quickstart

User Guide
----------

Learn how to use CNit effectively with examples and best practices.

.. toctree::
   :maxdepth: 2

   example

API Reference
-------------

Complete documentation of all classes, methods, and modules.  The API documentation
includes detailed scientific descriptions, equations, and implementation details.

.. toctree::
   :maxdepth: 2

   api/cnit.physics

Development
-----------

..  toctree::
   :maxdepth: 2

   changelog
   contributing

Contact
-------

CNit is developed and maintained by **Gang Tang**.

**Questions or feedback?**

* Open an issue:  `GitHub Issues <https://github.com/gangtang-au/cnit/issues>`_
* Email: gang.tang.au@gmail.com

Interested in contributing or collaboration?  Please get in touch!

License & Citation
------------------

CNit is released under the **BSD 3-Clause License**.

Copyright (c) 2026, Gang Tang and contributors.

See the `LICENSE file <https://github.com/yourusername/cnit/blob/main/LICENSE>`_
for full details.

**If you use CNit in your research, please cite our papers.**

* Tang, G., Nicholls, Z., Norton, A., Zaehle, S., and Meinshausen, M.: Synthesizing
  global carbon–nitrogen coupling effects – the MAGICC coupled carbon–nitrogen
  cycle model v1.0, Geosci. Model Dev., 18, 2193–2230,
  https://doi.org/10.5194/gmd-18-2193-2025, 2025.

* Tang, G., Zaehle, S., Nicholls, Z., Norton, A., Ziehn, T., & Meinshausen, M.
  Understanding the Drivers of Carbon-Nitrogen Cycle Variability in CMIP6 ESMs
  with MAGICC CNit v2.0: Model and Calibration Updates. ESS Open Archive.
  June 16, 2025. https://doi.org/10.22541/essoar.175008280.09297369/v1 [accepted by
  Journal of Advances in Modeling Earth Systems (JAMES)]

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