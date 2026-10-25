Experiment config
=================

.. automodule:: cnit.configs.experiment_config
   :members:
   :undoc-members:
   :show-inheritance:



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
