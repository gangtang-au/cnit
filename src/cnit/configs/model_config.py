import csv
import math
from pathlib import Path

from attrs import define, field, fields
import pint

from ..utils.units import Q, check_units


@define
class CarbonNitrogenCycleModelConfig:
    """Parameters and initial pool values for a CNit model.

    Pass this configuration to
    :meth:`~cnit.physics.carbon_nitrogen_cycle.CarbonNitrogenCycleModel.from_config`,
    or create it from overrides with :meth:`from_dict` or from a calibration table
    with :meth:`from_parameter_csv`.

    Quantity fields require units with the correct physical dimensionality. For
    example, a turnover time must have time units and a pool size must have carbon- or
    nitrogen-mass units. The supplied unit does not need to use the exact spelling
    shown in the parameter descriptions, provided it is dimensionally compatible.

    Parameters are grouped by process:

    - Land use: ``frc_rgr``, ``tau_rgr``, ``method_cLU``, and ``method_nLU``.
    - CO₂ response: ``CO2ref``, ``CO2b``, the ``sns_CO2_*`` sensitivities,
      ``eff_CO2_sig_max``, and ``method_CO2_NPP``.
    - Temperature response: ``sns_*2dT`` sensitivities and ``method_dT_NPP``.
    - Nitrogen-deficit response: ``sns_*2PUdef`` sensitivities and the ``eff_*_max``
      and ``sns_*2NPPrd`` parameters controlling carbon effects on nitrogen acquisition.
    - Land-use response: ``sns_*2LUrgr`` sensitivities.
    - Carbon and nitrogen fluxes: baseline rates, stoichiometry, and allocation
      fractions beginning with ``NPP0``, ``LPR0``, ``PU0``, ``BNF0``, ``NetMIN`` or
      ``frc_``.
    - Pool dynamics: ``tau_Cpls*`` and ``tau_Npls*`` are turnover times; ``Cpls*0``
      and ``Npls*0`` are initial pool sizes.
    - Initialization and switches: ``time0``, ``switch_N``, ``switch_Cpls``, and
      ``switch_Npls``.

    Fractions are checked to lie between zero and one, and allocation pairs may not
    sum to more than one. Turnover times must be positive and initial pool sizes must
    be non-negative. The model equations and additional parameter meanings are
    documented with their calculators in :mod:`cnit.physics`.
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
    """Sensitivity of NPP to temperature changes, exponential form [K\ :sup:`-1`]
    """

    sns_NPP2dT_sig: pint.Quantity = field(
        default=Q(0, "1/K"), validator=check_units("1/K")
    )
    """Sensitivity of NPP to temperature changes, sigmoid form [K\ :sup:`-1`]
    """

    method_dT_NPP: pint.Quantity = field(default=Q(0, "1"), validator=check_units("1"))
    """Blending parameter between exponential (0) and sigmoid (1) NPP responses
    [dimensionless, 0-1] [dimensionless]
    """

    sns_LPR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of litter production respiration to temperature changes
    [K\ :sup:`-1`]
    """

    sns_cLP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of carbon litter production to temperature changes [K\ :sup:`-1`]
    """

    sns_cLD2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of carbon litter decomposition to temperature changes [K\ :sup:`-1`]
    """

    sns_cSR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of carbon soil respiration (soil organic matter decomposition) to
    temperature changes [K\ :sup:`-1`]
    """

    sns_PU2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen plant uptake to temperature changes [K\ :sup:`-1`]
    """

    sns_BNF2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of biological nitrogen fixation to temperature changes [K\ :sup:`-1`]
    """

    sns_nLP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen litter production to temperature changes [K\ :sup:`-1`]
    """

    sns_nLD2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen litter decomposition to temperature changes
    [K\ :sup:`-1`]
    """

    sns_nSR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    """Sensitivity of nitrogen soil respiration (soil organic matter decomposition) to
    temperature changes [K\ :sup:`-1`]
    """

    sns_nLSgas2dT: pint.Quantity = field(
        default=Q(0, "1/K"), validator=check_units("1/K")
    )
    """Sensitivity of mineral nitrogen gaseous loss to temperature changes [K\
    :sup:`-1`]
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
            If the model column is unavailable, a table parameter is not a config
            field, or a value cannot be parsed or validated.
        """
        parameter_path = Path(csv_file)
        csv_context = parameter_path.open("r", encoding="utf-8", newline="")
        config_fields = {item.name for item in fields(cls)}

        with csv_context:
            reader = csv.DictReader(csv_context)
            if not reader.fieldnames or model_name not in reader.fieldnames:
                available = reader.fieldnames[2:] if reader.fieldnames else []
                raise ValueError(
                    f"Unknown model {model_name!r}. Available models: "
                    f"{', '.join(available)}"
                )

            values = {}
            for row in reader:
                name = row["parameter"].strip()
                if name not in config_fields:
                    raise ValueError(
                        f"Parameter {name!r} in {parameter_path.name} is not a "
                        "field of CarbonNitrogenCycleModelConfig"
                    )
                raw_value = row[model_name].strip()
                if not raw_value:
                    continue
                unit = row["unit"].strip()
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
