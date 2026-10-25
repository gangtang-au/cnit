"""Configuration of CNit model parameters and initial pool sizes.

The detailed parameter reference and attribute documentation are in
:class:`CarbonNitrogenCycleModelConfig`.
"""

import csv
import math
from pathlib import Path

from attrs import define, field, fields
import numpy as np
import pint

from ..utils.units import Q, check_units


@define
class CarbonNitrogenCycleModelConfig:
    """Stores shared model parameters and initial pool sizes for CNit runs.

    Use this configuration with
    :meth:`~cnit.physics.carbon_nitrogen_cycle.CarbonNitrogenCycleModel.from_config`.
    The public alias is ``CNitModelConfig``, available from the package root::

        from cnit import CNitModelConfig

    Create a config with defaults, override selected values with :meth:`from_dict`, or
    load values from a calibration table with :meth:`from_parameter_csv`. Unspecified
    parameters retain their defaults.

    Quantity fields require dimensionally compatible units. The units shown in the
    parameter tables are recommended units; equivalent units are accepted. On creation,
    the config checks fraction ranges, allocation sums, positive turnover times, and
    non-negative initial pool sizes.

    The tables summarize defaults and parameter meanings. The class attributes below
    provide the detailed attribute documentation, and the physics API documents the
    corresponding process equations.

    **Land-use parameters**

    These parameters select the land-use input method and control the duration and fraction of deforestation-related regrowth.

    .. list-table:: Land-use parameters
       :class: config-table
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

    .. list-table:: Land-use response parameters
       :class: config-table
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

    **Environmental response parameters**

    These parameters define how atmospheric CO₂ and temperature affect carbon and nitrogen process rates.

    .. list-table:: CO₂ response parameters
       :class: config-table
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

    .. list-table:: Temperature response parameters
       :class: config-table
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

    **Carbon–nitrogen coupling parameters**

    These parameters describe how plant nitrogen deficit and NPP-dependent effects influence carbon processes, nitrogen uptake, and biological fixation.

    .. list-table:: Responses to plant nitrogen deficit
       :class: config-table
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

    .. list-table:: NPP effects on plant nitrogen acquisition
       :class: config-table
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

    **Carbon-cycle parameters**

    These parameters set reference carbon fluxes, pool turnover times, within-timestep allocation, and initial carbon pool sizes.

    .. list-table:: Carbon process rates and turnover times
       :class: config-table
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

    .. list-table:: Carbon allocation fractions
       :class: config-table
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

    .. list-table:: Initial carbon pools
       :class: config-table
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

    **Nitrogen-cycle parameters**

    These parameters define nitrogen supply and mineralization, pool turnover,
    transfers and losses, and initial pool sizes.

    .. list-table:: Nitrogen inputs and mineralization
       :class: config-table
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

    .. list-table:: Nitrogen pool turnover times
       :class: config-table
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

    .. list-table:: Nitrogen allocation and loss fractions
       :class: config-table
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

    .. list-table:: Initial nitrogen pools
       :class: config-table
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

    **Initialization and pool-switch parameters**

    These parameters set the initial-condition year and enable or scale individual carbon and nitrogen pools.

    .. list-table:: Initialization and pool switches
       :class: config-table
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
    """

    frc_rgr: pint.Quantity = field(default=Q(0.9, "1"), validator=check_units("1"))
    """Fraction of deforestation that can regrow [dimensionless].

    Represents the portion of deforested area that can regrow naturally.
    The remaining fraction (1 - frc_rgr) represents permanent land conversion with no
    natural regrowth.
    """

    tau_rgr: pint.Quantity = field(default=Q(100, "yr"), validator=check_units("yr"))
    """Regrowth timescale [yr].

    Time period over which regrowth occurs following deforestation.
    Regrowth is distributed as a constant flux over this period.
    """

    method_cLU: pint.Quantity = field(default=Q(0, "1"), validator=check_units("1"))
    """Carbon land use flux calculation method [dimensionless].

    Selects the input data type and calculation approach:

    - 0: Use net emissions (CemsLUnet) as input, suitable for AFOLU emissions data
      where regrowth must be calculated from net emissions
    - 1: Use gross emissions (CemsLUgrs) as input, suitable for harvest/deforestation
      data where gross emissions are directly observed
    """

    method_nLU: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    """Nitrogen land use flux calculation method [dimensionless].

    Selects the input data type and calculation approach:

    - 0: Use net emissions (NemsLUnet) as input (less commonly used)
    - 1: Use gross emissions (NemsLUgrs) as input, suitable for gross deforestation
      data (typical for nitrogen)
    """

    CO2ref: pint.Quantity = field(
        default=Q(284.875, "ppm"), validator=check_units("ppm")
    )
    """Atmospheric CO₂ concentration at pre-industrial level [ppm].
    """

    CO2b: pint.Quantity = field(default=Q(31, "ppm"), validator=check_units("ppm"))
    """Atmospheric CO₂ concentration when NPP = 0 (for rectangular hyperbolic effect,
    default 31 ppm) [ppm].
    """

    sns_CO2_log: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    """Sensitivity of NPP to atmospheric CO₂ concentration changes
    (logarithmic form) [dimensionless].
    """

    sns_CO2_sig: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    """Sensitivity of NPP to atmospheric CO₂ concentration changes
    (sigmoid form) [dimensionless].
    """

    eff_CO2_sig_max: pint.Quantity = field(
        default=Q(2, "1"), validator=check_units("1")
    )
    """Maximum CO₂ effect factor for the sigmoid form [dimensionless].
    """

    method_CO2_NPP: pint.Quantity = field(default=Q(0, "1"), validator=check_units("1"))
    """CO₂ effect factor calculation method (value from 0 to 2),
    otherwise eff_CO2 is set equal to 1.
    """

    sns_NPP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of NPP to temperature changes, exponential form [K⁻¹]
    """

    sns_NPP2dT_sig: pint.Quantity = field(
        default=Q(0, "1/K"), validator=check_units("1/K")
    )
    """Sensitivity of NPP to temperature changes, sigmoid form [K⁻¹]
    """

    method_dT_NPP: pint.Quantity = field(default=Q(0, "1"), validator=check_units("1"))
    """Blending parameter between exponential (0) and sigmoid (1) NPP responses
    [dimensionless]
    """

    sns_LPR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of litter production respiration to temperature changes
    [K⁻¹]
    """

    sns_cLP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of carbon litter production to temperature changes [K⁻¹]
    """

    sns_cLD2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of carbon litter decomposition to temperature changes [K⁻¹]
    """

    sns_cSR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of carbon soil respiration (soil organic matter decomposition) to
    temperature changes [K⁻¹]
    """

    sns_PU2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen plant uptake to temperature changes [K⁻¹]
    """

    sns_BNF2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of biological nitrogen fixation to temperature changes [K⁻¹]
    """

    sns_nLP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen litter production to temperature changes [K⁻¹]
    """

    sns_nLD2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen litter decomposition to temperature changes [K⁻¹]
    """

    sns_nSR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen soil respiration (soil organic matter decomposition) to
    temperature changes [K⁻¹]
    """

    sns_nLSgas2dT: pint.Quantity = field(
        default=Q(0, "1/K"), validator=check_units("1/K")
    )
    """Sensitivity of mineral nitrogen gaseous loss to temperature changes [K⁻¹]
    """

    sns_NPP2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of net primary production to nitrogen plant uptake deficit [yr/GtN]
    """

    sns_cLP2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of litter production (carbon) to nitrogen plant uptake
    deficit [yr/GtN]
    """

    sns_cLD2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of litter decomposition (carbon) to nitrogen plant uptake
    deficit [yr/GtN]
    """

    sns_cSR2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of soil respiration (carbon) to nitrogen plant uptake
    deficit [yr/GtN]
    """

    sns_PU2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of plant uptake to N deficit.

    Controls how strongly plant uptake responds to nitrogen deficiency [yr/GtN].
    Higher values indicate stronger uptake response to N limitation.
    """

    sns_BNF2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of biological nitrogen fixation to N deficit.

    Controls how strongly fixation responds to nitrogen deficiency [yr/GtN].
    Higher values indicate stronger fixation response to N limitation.
    """

    sns_nLP2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of litter production (nitrogen) to nitrogen plant uptake
    deficit [yr/GtN]
    """

    sns_nLD2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of litter decomposition (nitrogen) to nitrogen plant uptake
    deficit [yr/GtN]
    """

    sns_nSR2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of soil respiration (nitrogen) to nitrogen plant uptake
    deficit [yr/GtN]
    """

    eff_C_PU_max: pint.Quantity = field(default=Q(2, "1"), validator=check_units("1"))
    """Maximum effect of NPP on plant uptake.

    Upper limit of the NPP-dependent multiplier on plant uptake rate [dimensionless].
    Values > 1 allow enhanced uptake under high NPP conditions.
    """

    sns_PU2NPPrd: pint.Quantity = field(default=Q(2.5, "1"), validator=check_units("1"))
    """Sensitivity of plant uptake to NPP.

    Controls how strongly plant uptake responds to changes in NPP [dimensionless].
    Higher values indicate stronger coupling between NPP and uptake.
    """

    eff_C_BNF_max: pint.Quantity = field(default=Q(2, "1"), validator=check_units("1"))
    """Maximum effect of NPP on biological nitrogen fixation.

    Upper limit of the NPP-dependent multiplier on fixation rate [dimensionless].
    Values > 1 allow enhanced fixation under high NPP conditions.
    """

    sns_BNF2NPPrd: pint.Quantity = field(
        default=Q(2.5, "1"), validator=check_units("1")
    )
    """Sensitivity of biological nitrogen fixation to NPP.

    Controls how strongly fixation responds to changes in NPP [dimensionless].
    Higher values indicate stronger coupling between NPP and fixation.
    """

    sns_cLP2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"), validator=check_units("yr/GtC")
    )
    """Sensitivity of carbon litter production to regrowth flux [yr/GtC]."""

    sns_cLD2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"), validator=check_units("yr/GtC")
    )
    """Sensitivity of carbon litter decomposition to regrowth flux [yr/GtC]."""

    sns_cSR2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"), validator=check_units("yr/GtC")
    )
    """Sensitivity of carbon soil respiration to regrowth flux [yr/GtC]"""

    sns_nLP2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"), validator=check_units("yr/GtC")
    )
    """Sensitivity of nitrogen litter production to regrowth flux [yr/GtC]"""

    sns_nLD2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"), validator=check_units("yr/GtC")
    )
    """Sensitivity of nitrogen litter decomposition to regrowth flux [yr/GtC]"""

    sns_nSR2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"), validator=check_units("yr/GtC")
    )
    """Sensitivity of nitrogen soil respiration to regrowth flux [yr/GtC]"""

    NPP0: pint.Quantity = field(
        default=Q(50, "GtC/yr"), validator=check_units("GtC/yr")
    )
    """Base net primary production (NPP) rate under reference conditions, before
    environmental effects [GtC/yr].
    """

    LPR0: pint.Quantity = field(default=Q(0, "GtC/yr"), validator=check_units("GtC/yr"))
    """Base litter production respiration (LPR) rate under reference conditions,
    before environmental effects [GtC/yr].
    """

    frc_NPP2P: pint.Quantity = field(default=Q(0.55, "1"), validator=check_units("1"))
    """Fraction of net primary production allocated to plant pool [dimensionless].

    Represents the portion of newly fixed carbon remaining in plant biomass at the
    end of a model timestep, primarily in long-lived tissues (wood, structural
    tissues, storage organs) that do not turn over within the model timestep.
    """

    frc_NPP2L: pint.Quantity = field(default=Q(0.35, "1"), validator=check_units("1"))
    """Fraction of net primary production allocated to litter pool [dimensionless].

    Represents newly fixed carbon that enters the litter pool within the model
    timestep through rapid plant-to-litter turnover (e.g., fine root turnover, leaf
    senescence). The remainder (1 - frc_NPP2P - frc_NPP2L) goes to the soil pool,
    representing carbon that cycles through both plant and litter pools and enters
    soil within the model timestep.
    """

    frc_cLP2L: pint.Quantity = field(default=Q(0.7, "1"), validator=check_units("1"))
    """Fraction of litter production carbon retained in litter pool [dimensionless].

    Represents carbon in plant materials that remain in the litter pool at the end
    of a model timestep as recognizable plant structures.  The remainder (1 - frc_cLP2L)
    goes directly to soil, representing materials that decompose rapidly and enter
    the soil pool within the model timestep, bypassing extended residence in the
    litter stage.
    """

    frc_cLD2S: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))
    """Fraction of litter decomposition carbon transferred to soil pool [dimensionless].

    Represents carbon incorporated into stable soil organic matter during decomposition.
    The remainder (1 - frc_cLD2S) is respired to the atmosphere (cLD2A), representing
    carbon released as CO₂ during litter breakdown.
    """

    frc_cLUgrs2P: pint.Quantity = field(default=Q(0.6, "1"), validator=check_units("1"))
    """Fraction of gross deforestation carbon loss from plant pool [dimensionless].

    Represents the distribution of carbon lost during land conversion, with this
    fraction coming from plant biomass removal or burning.
    """

    frc_cLUgrs2L: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))
    """Fraction of gross deforestation carbon loss from litter pool [dimensionless].

    Represents carbon lost from litter during land conversion. The remainder
    (1 - frc_cLUgrs2P - frc_cLUgrs2L) comes from the soil pool.
    """

    tau_CplsP: pint.Quantity = field(default=Q(10, "yr"), validator=check_units("yr"))
    """Plant carbon pool turnover time under reference conditions [yr].

    Determines the baseline rate at which carbon is transferred from plant
    biomass to litter through senescence and mortality.
    """

    tau_CplsL: pint.Quantity = field(default=Q(1, "yr"), validator=check_units("yr"))
    """Litter carbon pool turnover time under reference conditions [yr].

    Determines the baseline rate at which carbon is released from litter
    and transferred to soil organic matter through decomposition.
    """

    tau_CplsS: pint.Quantity = field(default=Q(100, "yr"), validator=check_units("yr"))
    """Soil carbon pool turnover time under reference conditions [yr].

    Determines the baseline rate at which carbon is released from soil
    organic matter through decomposition and respiration processes.
    """

    CplsP0: pint.Quantity = field(default=Q(400, "GtC"), validator=check_units("GtC"))
    """Initial plant carbon pool size [GtC].

    Carbon contained in living plant biomass at the start of the simulation
    (time0).
    """

    CplsL0: pint.Quantity = field(default=Q(50, "GtC"), validator=check_units("GtC"))
    """Initial litter carbon pool size [GtC].

    Carbon contained in dead plant material (litter) at the start of the simulation
    (time0).
    """

    CplsS0: pint.Quantity = field(default=Q(1550, "GtC"), validator=check_units("GtC"))
    """Initial soil carbon pool size [GtC].

    Carbon contained in soil organic matter at the start of the simulation (time0).
    """

    PU0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"), validator=check_units("GtN/yr")
    )
    """Base plant uptake rate under reference conditions, before environmental
    effects [GtN/yr].
    """

    BNF0: pint.Quantity = field(
        default=Q(0.1, "GtN/yr"), validator=check_units("GtN/yr")
    )
    """Base biological nitrogen fixation rate under reference conditions, before
    environmental effects [GtN/yr].
    """

    NCratio_PUreq: pint.Quantity = field(
        default=Q(0.02, "GtN/GtC"), validator=check_units("GtN/GtC")
    )
    """Required nitrogen-to-carbon ratio for plant uptake [GtN/GtC].  Stoichiometric
    ratio determining nitrogen demand based on carbon uptake.  This represents the
    N:C ratio needed to support new plant biomass production.
    """

    NetMIN0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"), validator=check_units("GtN/yr")
    )
    """Initial net mineralization rate under reference conditions [GtN/yr].
    Net mineralization is the balance between nitrogen release from organic matter
    decomposition (mineralization) and nitrogen uptake by microbes (immobilization).
    This flux is temperature-sensitive.
    """

    sns_NetMIN2dT: pint.Quantity = field(
        default=Q(0, "1/K"), validator=check_units("1/K")
    )
    """Temperature sensitivity of net mineralization [1/K].  Controls how strongly
    mineralization responds to temperature changes. Higher values indicate stronger
    temperature effect on mineralization rates.  Positive values indicate increased
    mineralization with warming.
    """

    NetMINlt0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"), validator=check_units("GtN/yr")
    )
    """Long-term net mineralization rate under reference conditions [GtN/yr].
    This component is sensitive to nitrogen deposition and changes in net primary
    productivity, representing slower-turnover mineral nitrogen pools.
    """

    sns_NetMINlt2AD: pint.Quantity = field(
        default=Q(0, "yr/GtN"), validator=check_units("yr/GtN")
    )
    """Sensitivity of long-term net mineralization to nitrogen deposition [yr/GtN].
    Controls how strongly the long-term mineralization rate responds to the
    atmospheric nitrogen deposition rate. This represents priming effects where
    added nitrogen can stimulate or suppress organic matter decomposition.
    """

    tau_NplsP: pint.Quantity = field(default=Q(10, "yr"), validator=check_units("yr"))
    """Plant nitrogen pool turnover time under reference conditions [yr].

    Determines the baseline rate at which nitrogen is transferred from plant
    biomass to litter through senescence and mortality.
    """

    tau_NplsL: pint.Quantity = field(default=Q(1, "yr"), validator=check_units("yr"))
    """Litter nitrogen pool turnover time under reference conditions [yr].

    Determines the baseline rate at which nitrogen is released from litter
    and transferred to soil organic matter through decomposition.
    """

    tau_NplsS: pint.Quantity = field(default=Q(100, "yr"), validator=check_units("yr"))
    """Soil nitrogen pool turnover time under reference conditions [yr].

    Determines the baseline rate at which nitrogen is released from soil
    organic matter through decomposition and mineralization processes.
    """

    tau_NplsM: pint.Quantity = field(default=Q(1, "yr"), validator=check_units("yr"))
    """Mineral nitrogen pool turnover time under reference conditions [yr].

    Determines the baseline rate at which mineral nitrogen is lost from the
    system through gaseous emissions (denitrification, volatilization) and
    leaching.
    """

    frc_BNF2P: pint.Quantity = field(default=Q(0.55, "1"), validator=check_units("1"))
    """Fraction of biological nitrogen fixation allocated to plant pool [dimensionless].

    Represents the portion of newly fixed nitrogen directly incorporated into plant
    biomass through symbiotic associations (e.g., legume-rhizobia symbiosis) or
    free-living fixers in the rhizosphere.
    """

    frc_BNF2L: pint.Quantity = field(default=Q(0.35, "1"), validator=check_units("1"))
    """Fraction of biological nitrogen fixation allocated to litter pool
    [dimensionless].

    Represents newly fixed nitrogen that enters the litter pool within the model
    timestep through rapid plant-to-litter turnover (e.g., fine root turnover, leaf
    senescence of N-fixing plants). The remainder (1 - frc_BNF2P - frc_BNF2L) goes
    to the soil pool, representing nitrogen that cycles through both plant and litter
    pools and enters soil within the model timestep.
    """

    frc_PU2P: pint.Quantity = field(default=Q(0.55, "1"), validator=check_units("1"))
    """Fraction of plant nitrogen uptake allocated to plant pool [dimensionless].

    Represents nitrogen remaining in plant biomass at the end of a model timestep,
    primarily in long-lived tissues (wood, structural tissues, storage organs) that
    do not turn over within the model timestep.
    """

    frc_PU2L: pint.Quantity = field(default=Q(0.35, "1"), validator=check_units("1"))
    """Fraction of plant nitrogen uptake allocated to litter pool [dimensionless].

    Represents nitrogen that cycles from mineral pool to plant and then to litter
    within the model timestep through rapid turnover of short-lived tissues (fine
    roots, senescing leaves, reproductive structures). The remainder (1 - frc_PU2P -
    frc_PU2L) goes to the soil pool, representing nitrogen that cascades through
    plant, litter, and into soil within the model timestep.
    """

    frc_nLP2L: pint.Quantity = field(default=Q(0.7, "1"), validator=check_units("1"))
    """Fraction of litter production nitrogen retained in litter pool [dimensionless].

    Represents nitrogen in plant materials that remain in the litter pool at the end
    of a model timestep as recognizable plant structures. The remainder (1 - frc_nLP2L)
    goes directly to soil, representing materials that decompose rapidly and enter
    the soil pool within the model timestep, bypassing extended residence in the
    litter stage.
    """

    frc_nLD2S: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))
    """Fraction of litter decomposition nitrogen transferred to soil pool
    [dimensionless].

    Represents nitrogen incorporated into stable soil organic matter during
    decomposition. The remainder (1 - frc_nLD2S) is mineralized to the mineral pool
    (nLD2M), representing nitrogen released as plant-available forms (NH₄⁺,
    NO₃⁻) during litter breakdown.
    """

    frc_nLSgas: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    """Fraction of mineral nitrogen loss occurring as gaseous loss [dimensionless].
    Partitions total mineral nitrogen loss between gaseous pathways (denitrification,
    volatilization) and leaching. Values range from 0 (all leaching) to 1 (all
    gaseous).
    Only gaseous losses are temperature-sensitive in this formulation.
    """

    frc_nLUgrs2P: pint.Quantity = field(default=Q(0.6, "1"), validator=check_units("1"))
    """Fraction of gross deforestation nitrogen loss from plant pool [dimensionless].

    Represents the distribution of nitrogen lost during land conversion, with this
    fraction coming from plant biomass removal or burning.
    """

    frc_nLUgrs2L: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))
    """Fraction of gross deforestation nitrogen loss from litter pool [dimensionless].

    Represents nitrogen lost from litter during land conversion. The remainder
    (1 - frc_nLUgrs2P - frc_nLUgrs2L) comes from the soil pool.
    """

    NplsP0: pint.Quantity = field(default=Q(5, "GtN"), validator=check_units("GtN"))
    """Initial plant nitrogen pool size [GtN].

    Nitrogen contained in living plant biomass at the start of the simulation
    (time0).
    """

    NplsL0: pint.Quantity = field(default=Q(1, "GtN"), validator=check_units("GtN"))
    """Initial litter nitrogen pool size [GtN].

    Nitrogen contained in dead plant material (litter) at the start of the simulation
    (time0).
    """

    NplsS0: pint.Quantity = field(default=Q(150, "GtN"), validator=check_units("GtN"))
    """Initial soil nitrogen pool size [GtN].

    Nitrogen contained in soil organic matter at the start of the simulation (time0).
    """

    NplsM0: pint.Quantity = field(default=Q(1, "GtN"), validator=check_units("GtN"))
    """Initial mineral nitrogen pool size [GtN].

    Plant-available inorganic nitrogen (primarily NH₄⁺ and NO₃⁻) at the start of
    the simulation (time0). This pool has much faster turnover than organic pools.
    """

    time0: pint.Quantity = field(default=Q(1850, "yr"), validator=check_units("yr"))
    """Initialization time [yr].

    Time point at which the initial pool sizes apply. The simulation time_axis must
    start at or after this time. The default is the pre-industrial baseline year
    (1850).
    """

    switch_N: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    """Switch to enable/disable nitrogen cycle [dimensionless].

    Controls whether nitrogen limitation effects are applied to carbon cycle:
    - 1: Full carbon-nitrogen coupling (default)
    - 0: Carbon cycle only (no nitrogen limitation)

    When set to 0, the model runs as a carbon-only model with no nitrogen feedbacks.
    """

    switch_Cpls: list = field(factory=lambda: [1, 1, 1])
    """Carbon pool switches for [plant, litter, soil] [dimensionless].

    Use 1 to enable or 0 to freeze a pool at its initial value. Other numeric values
    scale that pool's derivative. These switches are useful for diagnostic purposes
    or simplified model configurations.

    - switch_Cpls[0]: Plant pool switch
    - switch_Cpls[1]: Litter pool switch
    - switch_Cpls[2]: Soil pool switch
    """

    switch_Npls: list = field(factory=lambda: [1, 1, 1, 1])
    """Nitrogen pool switches for [plant, litter, soil, mineral] [dimensionless].

    For plant, litter, and soil pools, use 1 to enable or 0 to freeze the pool at its
    initial value; other numeric values scale its derivative. For the mineral pool,
    0 disables its dynamics and any non-zero value enables them. These switches are
    useful for diagnostic purposes or simplified model configurations.

    - switch_Npls[0]: Plant pool switch
    - switch_Npls[1]: Litter pool switch
    - switch_Npls[2]: Soil pool switch
    - switch_Npls[3]: Mineral pool switch
    """

    def __attrs_post_init__(self) -> None:
        """Validate physical parameter ranges when the config is created."""
        for name, expected_length in (("switch_Cpls", 3), ("switch_Npls", 4)):
            switches = getattr(self, name)
            if not isinstance(switches, (list, tuple, np.ndarray)):
                raise TypeError(f"{name} must be a sequence of numeric switches.")
            try:
                values = np.asarray(switches)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{name} must contain {expected_length} finite numeric values."
                ) from exc

            if values.ndim != 1 or values.size != expected_length:
                raise ValueError(
                    f"{name} must contain exactly {expected_length} values."
                )
            if values.dtype.kind not in "biuf" or not np.all(np.isfinite(values)):
                raise ValueError(f"{name} must contain only finite numeric values.")

        fraction_names = (
            "frc_rgr",
            "frc_NPP2P",
            "frc_NPP2L",
            "frc_cLP2L",
            "frc_cLD2S",
            "frc_cLUgrs2P",
            "frc_cLUgrs2L",
            "frc_BNF2P",
            "frc_BNF2L",
            "frc_PU2P",
            "frc_PU2L",
            "frc_nLP2L",
            "frc_nLD2S",
            "frc_nLSgas",
            "frc_nLUgrs2P",
            "frc_nLUgrs2L",
        )
        fractions = {}
        for name in fraction_names:
            value = float(getattr(self, name).to("1").m)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"{name} must be finite and between 0 and 1.")
            fractions[name] = value

        partition_pairs = (
            ("frc_NPP2P", "frc_NPP2L"),
            ("frc_cLUgrs2P", "frc_cLUgrs2L"),
            ("frc_BNF2P", "frc_BNF2L"),
            ("frc_PU2P", "frc_PU2L"),
            ("frc_nLUgrs2P", "frc_nLUgrs2L"),
        )
        for first, second in partition_pairs:
            if fractions[first] + fractions[second] > 1 + 1e-12:
                raise ValueError(f"{first} + {second} must not exceed 1.")

        time_names = (
            "tau_rgr",
            "tau_CplsP",
            "tau_CplsL",
            "tau_CplsS",
            "tau_NplsP",
            "tau_NplsL",
            "tau_NplsS",
            "tau_NplsM",
        )
        for name in time_names:
            value = float(getattr(self, name).to("yr").m)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a positive finite duration in years.")

        pool_names = (
            "CplsP0",
            "CplsL0",
            "CplsS0",
            "NplsP0",
            "NplsL0",
            "NplsS0",
            "NplsM0",
        )
        for name in pool_names:
            unit = "GtC" if name.startswith("C") else "GtN"
            value = float(getattr(self, name).to(unit).m)
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative.")

    @classmethod
    def from_dict(cls, params: dict) -> "CarbonNitrogenCycleModelConfig":
        """Create a config from parameter overrides.

        Parameters omitted from ``params`` retain their default values. Quantity
        values must include units, for example ``{"NPP0": Q(50, "GtC/yr")}``.
        Unknown parameter names and values that fail the config validators raise an
        error.

        Parameters
        ----------
        params
            Mapping from model parameter names to their values.

        Returns
        -------
        CarbonNitrogenCycleModelConfig
            Validated model configuration.
        """
        valid_fields = {item.name for item in fields(cls)}
        unknown = set(params) - valid_fields
        if unknown:
            raise ValueError(f"Unknown parameter(s) in config: {sorted(unknown)}")
        return cls(**params)

    @classmethod
    def from_parameter_csv(
        cls,
        csv_file,
        model_name: str,
    ) -> "CarbonNitrogenCycleModelConfig":
        """Build a model config from a user-selected calibration table.

        Values in the selected model column override defaults. Blank cells and
        parameters not listed in the table retain their defaults. The table must have
        ``parameter`` and ``unit`` columns, with one column for each model selection.
        Values are converted to quantities and validated when the config is created.

        Parameters
        ----------
        csv_file
            Path to the calibration CSV file, such as
            ``"data/cali_cmip6_best_fitting.csv"``.
        model_name
            Name of the model column to load, such as ``"ACCESS-ESM1-5"``.

        Returns
        -------
        CarbonNitrogenCycleModelConfig
            Validated model configuration with the selected calibration overrides.

        Raises
        ------
        ValueError
            If required columns are missing or duplicated, a parameter row is
            duplicated or invalid, the model column is unavailable, or a value
            cannot be parsed or validated.
        """
        parameter_path = Path(csv_file)
        csv_context = parameter_path.open("r", encoding="utf-8", newline="")
        config_fields = {item.name for item in fields(cls)}

        with csv_context:
            reader = csv.DictReader(csv_context)
            if not reader.fieldnames:
                raise ValueError(
                    f"Calibration file {parameter_path.name} has no header row."
                )

            duplicate_headers = {
                name for name in reader.fieldnames if reader.fieldnames.count(name) > 1
            }
            if duplicate_headers:
                raise ValueError(
                    f"Calibration file {parameter_path.name} has duplicate column "
                    f"names: {sorted(duplicate_headers)}"
                )

            required_headers = {"parameter", "unit"}
            missing_headers = required_headers - set(reader.fieldnames)
            if missing_headers:
                raise ValueError(
                    f"Calibration file {parameter_path.name} is missing required "
                    f"column(s): {', '.join(sorted(missing_headers))}"
                )

            if model_name not in reader.fieldnames:
                available = [
                    name for name in reader.fieldnames if name not in required_headers
                ]
                raise ValueError(
                    f"Unknown model {model_name!r}. Available models: "
                    f"{', '.join(available)}"
                )

            values = {}
            seen_parameters = set()
            for row in reader:
                name = (row.get("parameter") or "").strip()
                if not name:
                    raise ValueError(
                        f"Calibration file {parameter_path.name} has a row with an "
                        "empty parameter name."
                    )
                if name in seen_parameters:
                    raise ValueError(
                        f"Parameter {name!r} appears more than once in "
                        f"{parameter_path.name}."
                    )
                seen_parameters.add(name)

                if name not in config_fields:
                    raise ValueError(
                        f"Parameter {name!r} in {parameter_path.name} is not a "
                        "field of CarbonNitrogenCycleModelConfig"
                    )
                raw_value = (row.get(model_name) or "").strip()
                if not raw_value:
                    continue
                unit = (row.get("unit") or "").strip()
                if unit.lower() == "dimensionless":
                    unit = "1"
                try:
                    values[name] = Q(float(raw_value), unit)
                except (TypeError, ValueError, pint.errors.PintError) as exc:
                    raise ValueError(
                        f"Invalid value for parameter {name!r} and model "
                        f"{model_name!r}: {raw_value!r} {unit}"
                    ) from exc

        # Resolve table values against config defaults before deriving pool switches.
        config = cls.from_dict(values)
        carbon_active = config.CplsP0.m != 0
        config.switch_Cpls = [1, 1, 1] if carbon_active else [0, 0, 0]

        # Freeze nitrogen pools when nitrogen is disabled or NplsP0 is absent/zero.
        nitrogen_active = (
            "NplsP0" in values and config.switch_N.m != 0 and config.NplsP0.m != 0
        )
        config.switch_Npls = [1, 1, 1, 1] if nitrogen_active else [0, 0, 0, 0]

        return config
