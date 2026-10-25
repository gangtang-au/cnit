CNit Physics API
================

This section documents the main physics modules and submodules of the CNit package.
Each module contains classes and functions for modeling various aspects of the carbon and nitrogen cycles.

For model and experiment configuration classes, see :doc:`cnit.configs`.


Land Use
--------

Handles land use and land cover change processes.

.. automodule:: cnit.physics.landuse
   :members:
   :undoc-members:
   :show-inheritance:


Effects
-------

Implements physical and biogeochemical effect calculations.

.. automodule:: cnit.physics.effects
   :members:
   :undoc-members:
   :show-inheritance:


Carbon Cycle
------------

Classes and functions for carbon cycle modeling.

.. automodule:: cnit.physics.carbon_cycle
   :members:
   :undoc-members:
   :show-inheritance:


Nitrogen Cycle
--------------

Classes and functions for modeling nitrogen cycling.

.. automodule:: cnit.physics.nitrogen_cycle
   :members:
   :undoc-members:
   :show-inheritance:


Coupled Carbon-Nitrogen Cycle
-----------------------------

The coupled model combines the carbon and nitrogen cycle components. Its
implementation class, ``CarbonNitrogenCycleModel``, is exported from the package root
as ``CNitModel`` (``from cnit import CNitModel``). Use it with ``CNitModelConfig`` and
provide scenario forcings through ``CNitExpConfig``. See :doc:`cnit.configs` for
detailed parameter and input documentation. See :doc:`../quickstart` for a
configuration workflow.

.. autoclass:: cnit.physics.carbon_nitrogen_cycle.CarbonNitrogenCycleModel
   :exclude-members: from_config, run_experiments
   :show-inheritance:
