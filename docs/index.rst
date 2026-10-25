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

.. image:: https://static.pepy.tech/personalized-badge/cnit?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads
   :target: https://pepy.tech/projects/cnit
   :alt: PyPI downloads

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

   → See :doc:`CNit Model Overview <overview>` for the model’s scientific basis and :doc:`References <project_info/references>` for recommended publications.

   → See :doc:`CNit Physics API <api/cnit.physics>` for detailed scientific process descriptions and parameterizations.

**📚 For Model Users**
   Install CNit and run your first simulation.

   → See :doc:`quickstart` for installation and basic usage.

   → See :doc:`CMIP6 ESM emulation <examples/cmip6-esm>` and
   :doc:`TRENDY DGVM emulation <examples/trendy-dgvm>` for complete workflows,
   code, and outputs.

   → See :doc:`CNit Physics API <api/cnit.physics>` for detailed process equations and parameterizations.

   → See :doc:`api/cnit.configs` for model parameters and experiment inputs.

**🛠️ For Developers**
   Contribute to CNit development.

   → See :doc:`api/cnit.physics` for implementation details.

   → See :doc:`contributing` for development setup and coding guidelines.

Model Overview
--------------

Learn about CNit’s structure, processes, scientific basis, and intended use.

.. toctree::
   :caption: Model overview
   :maxdepth: 2

   overview

Getting Started
---------------

Install CNit, then follow the basic usage guide to set up and run a simulation.

.. toctree::
   :caption: Getting started
   :maxdepth: 2

   quickstart

API Reference
-------------

Look up detailed process parameterizations, model calculations, configuration
parameters, and class methods.

.. toctree::
   :caption: API reference
   :maxdepth: 2

   api/cnit.physics
   api/cnit.configs

Emulation Workflows
-------------------

Explore complete CMIP6 ESM and TRENDY DGVM emulation workflows, with notebooks and
model outputs.

.. toctree::
   :caption: Emulation Workflows
   :maxdepth: 2

   examples/cmip6-esm
   examples/trendy-dgvm

Development
-----------

Read the release history and find guidance for contributing to CNit.

.. toctree::
   :caption: Development
   :maxdepth: 2

   development/changelog
   contributing


Project Information
-------------------

Find maintainer contact details, licensing and citation information, and related projects.

.. toctree::
   :caption: Project information
   :maxdepth: 2

   project_info/references
   project_info/calibration-data
   project_info/license
   project_info/related-projects
   project_info/contact
