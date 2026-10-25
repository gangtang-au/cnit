Model config
============

.. automodule:: cnit.configs.model_config
   :members:
   :undoc-members:
   :show-inheritance:

.. code-block:: python

   from cnit import CNitModel, CNitModelConfig, Q

   config = CNitModelConfig.from_dict(
       {
           "NPP0": Q(50, "GtC/yr"),
           "tau_CplsP": Q(10, "yr"),
       }
   )
   model = CNitModel.from_config(config)

**Loading a calibration table**

The CSV loader uses the requested model name as a column heading. It reads each row's
``parameter``, ``unit``, and selected model value; blank values and parameters omitted
from the table keep their defaults.

.. code-block:: python

   from cnit import CNitModel, CNitModelConfig

   config = CNitModelConfig.from_parameter_csv(
       "data/cali_cmip6_best_fitting.csv",
       model_name="ACCESS-ESM1-5",
   )
   model = CNitModel.from_config(config)
