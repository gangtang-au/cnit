Model and Experiment Configuration
==================================

CNit separates **model parameters** from **experiment inputs**. A model
configuration sets parameter values and initial pools. An experiment configuration
provides the time axis and forcing series for one scenario. The same model instance
can run several experiments with different forcing data.

Model Parameters
----------------

Create a model configuration with defaults, override selected parameters with
``from_dict``, or load values from a calibration table with ``from_parameter_csv``.
Unspecified parameters retain their defaults. Quantity values must use compatible
units; the unit shown in each table is the recommended unit for reading and setting
that parameter. The config validates fraction ranges, allocation sums, turnover
times, and initial pool sizes when it is created.

.. code-block:: python

   from cnit import CNitModel, CNitModelConfig, Q

   config = CNitModelConfig.from_dict(
       {
           "NPP0": Q(50, "GtC/yr"),
           "tau_CplsP": Q(10, "yr"),
       }
   )
   model = CNitModel.from_config(config)

Land-use parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description
   * - ``frc_rgr``
     - 0.9
     - dimensionless
     - Fraction of deforestation that can regrow naturally; the remainder is permanent land
       conversion.
   * - ``tau_rgr``
     - 100
     - yr
     - Duration over which regrowth is distributed after deforestation; afforestation decay
       uses the same rounded duration.
   * - ``method_cLU``
     - 0
     - dimensionless
     - Selects net (0) or gross (1) carbon land-use emissions as input.
   * - ``method_nLU``
     - 1
     - dimensionless
     - Selects net (0) or gross (1) nitrogen land-use emissions as input.

CO₂ response parameters
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description
   * - ``CO2ref``
     - 284.875
     - ppm
     - Reference atmospheric CO₂ concentration.
   * - ``CO2b``
     - 31
     - ppm
     - CO₂ compensation point for the rectangular-hyperbolic NPP response.
   * - ``sns_CO2_log``
     - 1
     - dimensionless
     - Sensitivity of NPP to CO₂ in the logarithmic response.
   * - ``sns_CO2_sig``
     - 1
     - dimensionless
     - Sensitivity of NPP to CO₂ in the sigmoid response.
   * - ``eff_CO2_sig_max``
     - 2
     - dimensionless
     - Maximum effect factor in the sigmoid CO₂ response.
   * - ``method_CO2_NPP``
     - 0
     - dimensionless
     - Selects or blends the logarithmic, rectangular-hyperbolic, and sigmoid CO₂ responses;
       values outside 0–2 give an effect factor of 1.

Temperature response parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description
   * - ``sns_NPP2dT``
     - 0
     - 1/K
     - Sensitivity of NPP to temperature in the exponential response.
   * - ``sns_NPP2dT_sig``
     - 0
     - 1/K
     - Sensitivity of NPP to temperature in the sigmoid response.
   * - ``method_dT_NPP``
     - 0
     - dimensionless
     - Blends the exponential (0) and sigmoid (1) NPP temperature responses.
   * - ``sns_LPR2dT``
     - 0
     - 1/K
     - Sensitivity of litter production respiration to temperature.
   * - ``sns_cLP2dT``
     - 0
     - 1/K
     - Sensitivity of carbon litter production to temperature.
   * - ``sns_cLD2dT``
     - 0
     - 1/K
     - Sensitivity of carbon litter decomposition to temperature.
   * - ``sns_cSR2dT``
     - 0
     - 1/K
     - Sensitivity of carbon soil respiration to temperature.
   * - ``sns_PU2dT``
     - 0
     - 1/K
     - Sensitivity of plant nitrogen uptake to temperature.
   * - ``sns_BNF2dT``
     - 0
     - 1/K
     - Sensitivity of biological nitrogen fixation to temperature.
   * - ``sns_nLP2dT``
     - 0
     - 1/K
     - Sensitivity of nitrogen litter production to temperature.
   * - ``sns_nLD2dT``
     - 0
     - 1/K
     - Sensitivity of nitrogen litter decomposition to temperature.
   * - ``sns_nSR2dT``
     - 0
     - 1/K
     - Sensitivity of nitrogen soil respiration to temperature.
   * - ``sns_nLSgas2dT``
     - 0
     - 1/K
     - Sensitivity of gaseous mineral nitrogen loss to temperature.

Carbon–nitrogen coupling parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters describe how nitrogen availability and carbon supply
influence carbon-cycle processes, plant nitrogen uptake, and biological fixation.

Responses to plant nitrogen deficit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``sns_NPP2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of NPP to plant nitrogen uptake deficit.
   * - ``sns_cLP2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of carbon litter production to plant nitrogen uptake deficit.
   * - ``sns_cLD2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of carbon litter decomposition to plant nitrogen uptake deficit.
   * - ``sns_cSR2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of carbon soil respiration to plant nitrogen uptake deficit.
   * - ``sns_PU2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of plant nitrogen uptake to nitrogen deficit.
   * - ``sns_BNF2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of biological nitrogen fixation to nitrogen deficit.
   * - ``sns_nLP2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of nitrogen litter production to plant nitrogen uptake deficit.
   * - ``sns_nLD2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of nitrogen litter decomposition to plant nitrogen uptake deficit.
   * - ``sns_nSR2PUdef``
     - 0
     - yr/GtN
     - Sensitivity of nitrogen soil respiration to plant nitrogen uptake deficit.

NPP effects on plant nitrogen acquisition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``eff_C_PU_max``
     - 2
     - dimensionless
     - Maximum NPP-dependent effect on plant nitrogen uptake.
   * - ``sns_PU2NPPrd``
     - 2.5
     - dimensionless
     - Sensitivity of plant nitrogen uptake to relative NPP changes.
   * - ``eff_C_BNF_max``
     - 2
     - dimensionless
     - Maximum NPP-dependent effect on biological nitrogen fixation.
   * - ``sns_BNF2NPPrd``
     - 2.5
     - dimensionless
     - Sensitivity of biological nitrogen fixation to relative NPP changes.

Land-use response parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description
   * - ``sns_cLP2LUrgr``
     - 0
     - yr/GtC
     - Sensitivity of carbon litter production to regrowth flux.
   * - ``sns_cLD2LUrgr``
     - 0
     - yr/GtC
     - Sensitivity of carbon litter decomposition to regrowth flux.
   * - ``sns_cSR2LUrgr``
     - 0
     - yr/GtC
     - Sensitivity of carbon soil respiration to regrowth flux.
   * - ``sns_nLP2LUrgr``
     - 0
     - yr/GtC
     - Sensitivity of nitrogen litter production to regrowth flux.
   * - ``sns_nLD2LUrgr``
     - 0
     - yr/GtC
     - Sensitivity of nitrogen litter decomposition to regrowth flux.
   * - ``sns_nSR2LUrgr``
     - 0
     - yr/GtC
     - Sensitivity of nitrogen soil respiration to regrowth flux.

Carbon-cycle parameters
~~~~~~~~~~~~~~~~~~~~~~~

These parameters define reference carbon fluxes, pool turnover, transfers
between pools, and initial pool sizes.

Carbon process rates and turnover times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``NPP0``
     - 50
     - GtC/yr
     - Reference net primary production before environmental effects.

   * - ``LPR0``
     - 0
     - GtC/yr
     - Reference litter production respiration before environmental effects.

   * - ``tau_CplsP``
     - 10
     - yr
     - Plant carbon pool turnover time under reference conditions.

   * - ``tau_CplsL``
     - 1
     - yr
     - Litter carbon pool turnover time under reference conditions.

   * - ``tau_CplsS``
     - 100
     - yr
     - Soil carbon pool turnover time under reference conditions.

Carbon allocation fractions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``frc_NPP2P``
     - 0.55
     - dimensionless
     - Fraction of NPP retained in the plant pool at the end of a timestep.

   * - ``frc_NPP2L``
     - 0.35
     - dimensionless
     - Fraction of NPP allocated to litter through within-timestep plant turnover; the
       remainder goes to soil.

   * - ``frc_cLP2L``
     - 0.7
     - dimensionless
     - Fraction of litter production retained in litter; the remainder goes directly to soil.

   * - ``frc_cLD2S``
     - 0.3
     - dimensionless
     - Fraction of litter decomposition transferred to soil; the remainder is respired.

   * - ``frc_cLUgrs2P``
     - 0.6
     - dimensionless
     - Fraction of gross deforestation carbon loss from the plant pool.

   * - ``frc_cLUgrs2L``
     - 0.3
     - dimensionless
     - Fraction of gross deforestation carbon loss from the litter pool; the remainder comes
       from soil.

Initial carbon pools
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``CplsP0``
     - 400
     - GtC
     - Initial plant carbon pool size.

   * - ``CplsL0``
     - 50
     - GtC
     - Initial litter carbon pool size.

   * - ``CplsS0``
     - 1550
     - GtC
     - Initial soil carbon pool size.

Nitrogen-cycle parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters define nitrogen supply and mineralization, pool turnover,
transfers and losses, and initial pool sizes.

Nitrogen inputs and mineralization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``PU0``
     - 0.5
     - GtN/yr
     - Reference plant nitrogen uptake rate before environmental effects.

   * - ``BNF0``
     - 0.1
     - GtN/yr
     - Reference biological nitrogen fixation rate before environmental effects.

   * - ``NCratio_PUreq``
     - 0.02
     - GtN/GtC
     - Required nitrogen-to-carbon ratio used to calculate plant nitrogen demand.

   * - ``NetMIN0``
     - 0.5
     - GtN/yr
     - Reference net mineralization rate; balances organic nitrogen release and microbial
       uptake.

   * - ``sns_NetMIN2dT``
     - 0
     - 1/K
     - Temperature sensitivity of net mineralization.

   * - ``NetMINlt0``
     - 0.5
     - GtN/yr
     - Reference long-term net mineralization rate, also affected by deposition and NPP.

   * - ``sns_NetMINlt2AD``
     - 0
     - yr/GtN
     - Sensitivity of long-term net mineralization to atmospheric nitrogen deposition.

Nitrogen pool turnover times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``tau_NplsP``
     - 10
     - yr
     - Plant nitrogen pool turnover time under reference conditions.

   * - ``tau_NplsL``
     - 1
     - yr
     - Litter nitrogen pool turnover time under reference conditions.

   * - ``tau_NplsS``
     - 100
     - yr
     - Soil nitrogen pool turnover time under reference conditions.

   * - ``tau_NplsM``
     - 1
     - yr
     - Mineral nitrogen pool turnover time under reference conditions.

Nitrogen allocation and loss fractions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``frc_BNF2P``
     - 0.55
     - dimensionless
     - Fraction of biological nitrogen fixation allocated to the plant pool.

   * - ``frc_BNF2L``
     - 0.35
     - dimensionless
     - Fraction of biological nitrogen fixation allocated to litter through within-timestep
       plant turnover; the remainder goes to soil.

   * - ``frc_PU2P``
     - 0.55
     - dimensionless
     - Fraction of plant nitrogen uptake retained in plant biomass at the end of a timestep.

   * - ``frc_PU2L``
     - 0.35
     - dimensionless
     - Fraction of plant nitrogen uptake allocated to litter through within-timestep
       turnover; the remainder goes to soil.

   * - ``frc_nLP2L``
     - 0.7
     - dimensionless
     - Fraction of litter production nitrogen retained in litter; the remainder goes directly
       to soil.

   * - ``frc_nLD2S``
     - 0.3
     - dimensionless
     - Fraction of litter decomposition nitrogen transferred to soil organic matter; the
       remainder is mineralized.

   * - ``frc_nLSgas``
     - 1
     - dimensionless
     - Fraction of mineral nitrogen loss occurring as gaseous loss; the remainder is leached.

   * - ``frc_nLUgrs2P``
     - 0.6
     - dimensionless
     - Fraction of gross deforestation nitrogen loss from the plant pool.

   * - ``frc_nLUgrs2L``
     - 0.3
     - dimensionless
     - Fraction of gross deforestation nitrogen loss from the litter pool; the remainder
       comes from soil.

Initial nitrogen pools
^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description

   * - ``NplsP0``
     - 5
     - GtN
     - Initial plant nitrogen pool size.

   * - ``NplsL0``
     - 1
     - GtN
     - Initial litter nitrogen pool size.

   * - ``NplsS0``
     - 150
     - GtN
     - Initial soil nitrogen pool size.

   * - ``NplsM0``
     - 1
     - GtN
     - Initial mineral nitrogen pool size; this pool turns over faster than the organic
       pools.

Initialization and pool switches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 12 14 54

   * - Parameter
     - Default
     - Unit
     - Description
   * - ``time0``
     - 1850
     - yr
     - Time at which initial pool sizes apply; the time axis must start at or after this
       value.
   * - ``switch_N``
     - 1
     - dimensionless
     - Zero disables nitrogen-cycle coupling; any non-zero value enables it.
   * - ``switch_Cpls``
     - [1, 1, 1]
     - dimensionless
     - For plant, litter, and soil carbon pools, 1 enables and 0 freezes the pool; other
       values scale its derivative.
   * - ``switch_Npls``
     - [1, 1, 1, 1]
     - dimensionless
     - For plant, litter, and soil nitrogen pools, other values scale derivatives; the
       mineral pool is disabled at 0 and enabled otherwise.

Land-use method selectors accept only 0 or 1. An invalid selector raises an error
when non-zero land-use inputs are processed; all-zero flux inputs return zero fluxes.

Loading a calibration table
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Experiment Inputs
-----------------

An experiment configuration contains one named scenario's forcing data. Every forcing
series must have the same one-dimensional shape as ``time_axis``. The axis must contain
at least two finite, strictly increasing, evenly spaced time coordinates in years. The
input quantities are checked for compatible units when the config is created.

The model accepts evenly spaced subannual time axes. Forcing rates are treated as
constant over each interval, and allocation fractions apply within each model timestep.
The calibration files provided with CNit are based on an annual timescale.

.. list-table:: Experiment inputs
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
     - Evenly spaced time coordinates. Each value marks the start of an input interval.
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



Forcing values themselves are not checked for finiteness by this config. Net and gross
land-use series are both required by the config, although the model uses only the series
selected by the corresponding land-use method.

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

Configuration API
-----------------

``CNitModel``
   The main coupled model. Create it from a configuration with
   ``CNitModel.from_config``. Use ``run`` for one experiment or ``run_experiments``
   to run several experiments.

``CNitModelConfig``
   Stores model parameters and initial pool sizes. Use ``from_dict`` to override
   selected parameters or ``from_parameter_csv`` to load a calibration.

``CNitExpConfig``
   Stores one experiment's time axis and forcing series. Use ``to_input_dict`` to
   obtain the inputs for ``CNitModel.run``.
