Coupled Carbon-Nitrogen Cycle
=============================

The coupled model combines the carbon and nitrogen cycle components. Its
implementation class, ``CarbonNitrogenCycleModel``, is exported from the package
root as ``CNitModel`` (``from cnit import CNitModel``). Use it with
``CNitModelConfig`` and provide scenario input series through ``CNitExpConfig``.
See :doc:`../../quickstart` for a configuration workflow and
:doc:`../cnit.configs` for detailed parameter and input documentation.

.. autoclass:: cnit.physics.carbon_nitrogen_cycle.CarbonNitrogenCycleModel
   :exclude-members: from_config, run_experiments
   :show-inheritance:
