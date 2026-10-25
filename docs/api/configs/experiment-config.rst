Experiment config
=================

Import the experiment configuration as ``CNitExpConfig`` from the package root:
``from cnit import CNitExpConfig``. This public alias refers to the implementation
class ``CarbonNitrogenCycleExperimentConfig`` documented below.

An experiment configuration contains one named scenario's forcing data. Every forcing
series must have the same one-dimensional shape as ``time_axis``. On creation, the config
checks compatible units and verifies that the axis has at least two finite, strictly
increasing, evenly spaced time coordinates in years.

Evenly spaced subannual time axes are supported. Allocation fractions apply within each
model timestep, and the supplied calibration files use an annual timescale.

.. list-table:: Experiment inputs
   :class: config-table
   :header-rows: 1
   :widths: 28 17 55

   * - Input
     - Unit
     - Description
   * - ``name``
     - —
     - Scenario label stored in the experiment results.
   * - ``switch_N``
     - dimensionless
     - Enables nitrogen-cycle coupling for this experiment when 1; disables it when 0.
   * - ``time_axis``
     - yr
     - Evenly spaced time coordinates.
   * - ``dT_s``
     - K
     - Temperature anomaly relative to the model reference temperature.
   * - ``CO2_s``
     - ppm
     - Atmospheric CO₂ concentration.
   * - ``CemsLUnet_s``
     - GtC/yr
     - Net carbon land-use emission rate; used when ``method_cLU`` is 0.
   * - ``CemsLUgrs_s``
     - GtC/yr
     - Gross carbon land-use emission rate; used when ``method_cLU`` is 1.
   * - ``NflxAD_s``
     - GtN/yr
     - Nitrogen input rate from atmospheric deposition.
   * - ``NflxFT_s``
     - GtN/yr
     - Nitrogen input rate from fertilizer.
   * - ``NemsLUnet_s``
     - GtN/yr
     - Net organic nitrogen land-use emission rate; used when ``method_nLU`` is 0.
   * - ``NemsLUgrs_s``
     - GtN/yr
     - Gross organic nitrogen land-use emission rate; used when ``method_nLU`` is 1.
   * - ``NemsLUmin_s``
     - GtN/yr
     - Direct mineral nitrogen loss associated with land use.
   * - ``metadata``
     - —
     - Optional user metadata stored on the config and omitted from model inputs.



Net and gross land-use series are both required by the config, although the model uses
only the series selected by the corresponding land-use method.

Use an experiment with ``run_experiments``:

.. code-block:: python

   import numpy as np

   from cnit import (
       CNitExpConfig,
       CNitModel,
       CNitModelConfig,
       Q,
   )

   model = CNitModel.from_config(CNitModelConfig())
   time_axis = Q(np.arange(1850, 1853), "yr")
   n = time_axis.size
   experiment = CNitExpConfig(
       name="baseline",
       switch_N=Q(1, "1"),
       time_axis=time_axis,
       dT_s=Q(np.zeros(n), "K"),
       CO2_s=Q(np.linspace(280, 282, n), "ppm"),
       CemsLUnet_s=Q(np.zeros(n), "GtC/yr"),
       CemsLUgrs_s=Q(np.zeros(n), "GtC/yr"),
       NflxAD_s=Q(np.zeros(n), "GtN/yr"),
       NflxFT_s=Q(np.zeros(n), "GtN/yr"),
       NemsLUnet_s=Q(np.zeros(n), "GtN/yr"),
       NemsLUgrs_s=Q(np.zeros(n), "GtN/yr"),
       NemsLUmin_s=Q(np.zeros(n), "GtN/yr"),
   )
   results = model.run_experiments([experiment])



.. automodule:: cnit.configs.experiment_config
   :members:
   :undoc-members:
   :show-inheritance:
