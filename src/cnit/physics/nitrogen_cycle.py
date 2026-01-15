"""
MAGICC Terrestrial Nitrogen Cycle Model Implementation

This module provides the implementation of the nitrogen cycle component of MAGICC's
coupled carbon-nitrogen cycle model (CNit). It simulates the dynamics of terrestrial
nitrogen pools and their interactions with carbon cycling, climate change, atmospheric
deposition, fertilizer application, and land use change.

**Nitrogen Cycle Components:**

- :class:`NitrogenPUBNFCalculator`: Calculates plant nitrogen uptake (PU) and biological
  nitrogen fixation (BNF) fluxes with environmental modifiers and nitrogen deficit feedbacks
- :class:`NitrogenTurnoverCalculator`: Computes nitrogen turnover fluxes (litter production,
  litter decomposition, soil respiration, mineral losses) based on pool sizes and turnover times
- :class:`NitrogenCycleCalculator`: Handles nitrogen mass balance and flux partitioning among
  plant, litter, soil, and mineral pools using empirical allocation fractions
- :class:`NitrogenCycleModel`: Solves the nitrogen cycle ODEs in two stages using
  scipy.integrate.solve_ivp
- :class:`NitrogenCycleModelResult`: Stores and processes nitrogen cycle simulation results

**Key Model Features:**

1. **Four-Pool Structure**: The model tracks nitrogen in plant, litter, soil organic matter,
   and mineral (plant-available inorganic) pools, representing both organic and inorganic
   nitrogen reservoirs.

2. **Two-Stage Solution**: The nitrogen cycle is solved in two sequential stages:

   - Stage 1: Organic nitrogen pools (plant, litter, soil)
   - Stage 2: Mineral nitrogen pool using net mineralization from Stage 1

   This approach reflects the coupling between organic matter decomposition and mineral
   nitrogen availability.

3. **Nitrogen Deficit Feedback**: A key feature is the explicit calculation of nitrogen
   deficit for plant uptake:

   ..  math::

       PUdef = PUreq - PUavail

   where nitrogen requirement (PUreq) is based on NPP and stoichiometric N: C ratios, and
   nitrogen availability (PUavail) comes from mineralization and external inputs.  This
   deficit drives compensatory responses in plant uptake efficiency and biological fixation.

4. **Flux-Based Mineralization**: Net mineralization is parameterized as a flux-based process
   rather than relying on explicit mineral pool dynamics.  This approach is necessary for
   annual timestep integration, as the mineral nitrogen pool turnover time (days to weeks)
   is much shorter than the annual timestep.  The flux-based approach captures essential
   nitrogen limitation dynamics without modeling sub-annual mineral pool fluctuations.

5. **Annual Timestep with Cascading**: The model operates on annual timesteps but accounts
   for rapid sub-annual cycling through empirical allocation fractions that represent nitrogen
   cascading through multiple pools within a year (e.g., nitrogen taken up by plants, then
   entering litter through turnover, then entering soil through decomposition—all within
   one year).

6. **Environmental Modifiers**: All process rates are modified by dimensionless effect factors
   representing:

   - Temperature effects on turnover and mineralization
   - Nitrogen availability effects (carbon-nitrogen coupling)
   - Land use change effects on turnover rates

   (see :mod:`effect_factor_calculators`)

7. **Multiple Nitrogen Inputs**: The model includes both natural and anthropogenic nitrogen
   inputs:

   - Biological nitrogen fixation (BNF): Natural N₂ fixation by microorganisms
   - Atmospheric deposition (AD): Wet and dry deposition from atmospheric sources
   - Fertilizer application (FT): Direct anthropogenic nitrogen addition

8. **Nitrogen Losses**: The model represents nitrogen losses through:

   - Gaseous emissions:  Denitrification (N₂, N₂O, NO) and ammonia volatilization
   - Leaching: Nitrate leaching and dissolved organic nitrogen export
   - Land use change: Direct losses during deforestation and land conversion

**Model Structure:**

The nitrogen cycle includes four terrestrial pools:

- **Plant pool (NplsP)**: Nitrogen in living plant biomass (structural and storage compounds)
- **Litter pool (NplsL)**: Nitrogen in dead plant material undergoing decomposition
- **Soil pool (NplsS)**: Nitrogen in soil organic matter (proteins, amino acids, humus)
- **Mineral pool (NplsM)**: Plant-available inorganic nitrogen (primarily NH₄⁺ and NO₃⁻)

Key processes:

- **External inputs**: Biological nitrogen fixation (BNF), atmospheric deposition (AD),
  fertilizer application (FT)
- **Plant uptake (PU)**: Transfer from mineral to organic pools
- **Organic matter turnover**: Litter production (nLP), litter decomposition (nLD),
  soil respiration/turnover (nSR)
- **Mineralization**: Release of plant-available nitrogen from organic matter decomposition
  (NetMIN = nLD2M + nSR)
- **Nitrogen losses**: Gaseous emissions and leaching from mineral pool (nLS)
- **Land use emissions**: Nitrogen lost from organic pools (nLUgrs) and mineral pool (nLUmin)
  during land conversion
"""

from typing import Callable, Dict, List
from attrs import define, field
import numpy as np
import scipy
import xarray as xr
import pint

from ..utils.units import check_units, Q
from ..utils.maths import _get_interp, SolveError
from ..utils.data_utils import make_dataset_from_var_dict


@define
class NitrogenPUBNFCalculator:
    """Calculator for nitrogen plant uptake and biological nitrogen fixation.

    This class computes nitrogen fluxes in the terrestrial biosphere, including:
    - Plant uptake (PU)
    - Biological nitrogen fixation (BNF)

    The calculator accounts for multiple environmental controls including temperature,
    carbon cycle dynamics, nitrogen availability, nitrogen deposition,
    and anthropogenic inputs (fertilization). A key feature is the plant uptake deficit feedback, where
    insufficient nitrogen availability for plant uptake (positive deficit) can stimulate
    compensatory responses in both uptake efficiency and biological fixation.
    """

    PU0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Initial plant uptake rate, base rate of nitrogen uptake by plants 
    under reference conditions without any environmental effects [GtN/yr]. 
    """

    BNF0: pint.Quantity = field(
        default=Q(0.1, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Initial biological nitrogen fixation rate, base rate of nitrogen 
    fixation by symbiotic and free-living organisms under reference conditions 
    without any environmental effects [GtN/yr].
    """

    # Stoichiometric Ratios
    NCratio_PUreq: pint.Quantity = field(
        default=Q(0.02, "GtN/GtC"),
        validator=check_units("GtN/GtC"),
    )
    """Required nitrogen-to-carbon ratio for plant uptake [GtN/GtC].  Stoichiometric 
    ratio determining nitrogen demand based on carbon uptake.  This represents the 
    N:C ratio needed to support new plant biomass production.
    """

    NetMIN0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Initial net mineralization rate under reference conditions [GtN/yr]. 
    Net mineralization is the balance between nitrogen release from organic matter 
    decomposition (mineralization) and nitrogen uptake by microbes (immobilization). 
    This flux is temperature-sensitive. 
    """

    sns_NetMIN2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """Temperature sensitivity of net mineralization [1/K].  Controls how strongly 
    mineralization responds to temperature changes. Higher values indicate stronger 
    temperature effect on mineralization rates.  Positive values indicate increased 
    mineralization with warming.
    """

    NetMINlt0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Long-term net mineralization rate under reference conditions [GtN/yr]. 
    This component is sensitive to nitrogen deposition and changes in net primary 
    productivity, representing slower-turnover mineral nitrogen pools.
    """

    sns_NetMINlt2AD: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """Sensitivity of long-term net mineralization to nitrogen deposition [yr/GtN]. 
    Controls how strongly long-term mineralization responds to cumulative nitrogen 
    deposition.  This represents priming effects where added nitrogen can stimulate 
    or suppress organic matter decomposition.
    """

    def calculate(
            self,
            eff_dT_PU: pint.Quantity,
            eff_C_PU: pint.Quantity,
            eff_N_PU: pint.Quantity,
            eff_dT_BNF: pint.Quantity,
            eff_C_BNF: pint.Quantity,
            eff_N_BNF: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate nitrogen fluxes for plant uptake (PU) and biological nitrogen fixation (BNF).

        Parameters
        ----------
        eff_dT_PU
            Temperature effect on plant uptake [dimensionless].  Values > 1 indicate
            enhanced uptake with warming, < 1 indicate reduced uptake.
        eff_C_PU
            Carbon cycle effect on plant uptake [dimensionless].  Represents how
            changes in NPP affect nitrogen uptake capacity. This is a sigmoid function
            of relative NPP change (see ``calculate_eff_C_PUBNF``), where increased
            NPP enhances the ability to take up nitrogen up to a maximum value.
        eff_N_PU
            Nitrogen cycle effect on plant uptake [dimensionless]. Represents
            the response of plant nitrogen uptake to nitrogen limitation. This is
            an exponential function of nitrogen deficit (see ``calculate_eff_N_PUBNF``),
            where positive deficit (N limitation) can enhance uptake efficiency through
            mechanisms like increased root allocation, upregulation of N transporters,
            or mycorrhizal associations. With positive sensitivity (``sns_PU2PUdef > 0``),
            values > 1 indicate enhanced uptake under N limitation.
        eff_dT_BNF
            Temperature effect on biological nitrogen fixation [dimensionless].
            Values > 1 indicate enhanced fixation with warming.
        eff_C_BNF
            Carbon cycle effect on biological nitrogen fixation [dimensionless].
            Represents how changes in NPP affect BNF rates. BNF is energetically
            costly and requires carbon from photosynthesis. This is a sigmoid function
            of relative NPP change (see ``calculate_eff_C_PUBNF``), where increased
            NPP provides more carbon energy to support fixation up to a maximum value.
        eff_N_BNF
            Nitrogen cycle effect on biological nitrogen fixation [dimensionless].
            Represents the response of BNF to nitrogen limitation. This is an
            exponential function of nitrogen deficit (see ``calculate_eff_N_PUBNF``).
            With positive sensitivity (``sns_BNF2PUdef > 0``), N limitation (positive
            deficit) stimulates BNF as it becomes favorable to invest energy in fixation
            when soil N is scarce. With negative sensitivity, high N availability
            (negative deficit) can down-regulate fixation since it is energetically
            expensive and unnecessary when soil N is abundant.

        Returns
        -------
            Dictionary with keys:

            - NflxPU: Plant uptake nitrogen flux [GtN/yr]
            - NflxBNF: Biological nitrogen fixation flux [GtN/yr]

        Notes
        -----
        The following formulas are used:

        ..  math::

            PU = PU_0 \\times \\epsilon_{T(PU)} \\times \\epsilon_{C(PU)} \\times \\epsilon_{N(PU)}

            BNF = BNF_0 \\times \\epsilon_{T(BNF)} \\times \\epsilon_{C(BNF)} \\times \\epsilon_{N(BNF)}

        where:

        - :math:`PU_0` is the initial plant uptake rate
        - :math:`BNF_0` is the initial biological nitrogen fixation rate
        - :math:`\\epsilon_{T}` represents temperature effects
        - :math:`\\epsilon_{C}` represents carbon cycle effects (NPP-driven enhancement via sigmoid response)
        - :math:`\\epsilon_{N}` represents nitrogen deficit effects (exponential response to N limitation)

        The carbon effects (:math:`\\epsilon_{C(PU)}` and :math:`\\epsilon_{C(BNF)}`)
        are computed as sigmoid functions of relative NPP change, ensuring smooth
        transitions from reference conditions and saturation at high productivity levels.

        The nitrogen effects (:math:`\\epsilon_{N(PU)}` and :math:`\\epsilon_{N(BNF)}`)
        are exponential functions of the plant uptake deficit, representing compensatory
        responses when nitrogen availability for plant uptake is insufficient. When the
        deficit is positive (required > available for PU), these effects can enhance
        nitrogen acquisition, creating negative feedbacks that help maintain plant
        nitrogen supply.
        """

        return {
            "NflxPU": self.PU0 * eff_dT_PU * eff_C_PU * eff_N_PU,
            "NflxBNF": self.BNF0 * eff_dT_BNF * eff_C_BNF * eff_N_BNF
        }

    def calculate_NflxPU_req_avail_def(
            self,
            CflxNPP: pint.Quantity,
            CflxNPP0: pint.Quantity,
            NflxAD: pint.Quantity,
            NflxFT: pint.Quantity,
            dT: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate required, available, and deficit nitrogen for plant uptake.

        This method diagnoses the nitrogen budget specifically for plant uptake by
        comparing the nitrogen required for plant uptake (based on stoichiometric
        requirements for NPP) with the nitrogen available for plant uptake (from
        mineralization and external inputs). The plant uptake deficit determines
        whether plants have sufficient nitrogen to support their carbon uptake and
        is used to calculate feedback effects on various carbon-nitrogen processes.

        The calculation uses potential NPP (NPP without nitrogen limitation) to
        determine nitrogen requirement, ensuring that the deficit calculation is
        independent of the current nitrogen limitation state. This allows for
        proper feedback dynamics in carbon-nitrogen coupling.

        Parameters
        ----------
        CflxNPP
            Net primary production carbon flux [GtC/yr].  This should be potential
            NPP (NPPpot) calculated with no nitrogen limitation effect (ε_N(NPP) = 1)
            to properly determine nitrogen requirement independent of current limitation.
        CflxNPP0
            Initial/reference net primary production carbon flux [GtC/yr]. Used
            to normalize NPP changes for long-term mineralization calculation.
        NflxAD
            Nitrogen deposition flux [GtN/yr].  Atmospheric nitrogen deposition
            from anthropogenic and natural sources. All deposition is assumed to be
            immediately available for plant uptake.
        NflxFT
            Nitrogen fertilizer flux [GtN/yr].  Direct fertilizer application
            (currently not used in calculation but included for future extension).
        dT
            Temperature change relative to reference [K]. Used to calculate
            temperature-dependent mineralization.

        Returns
        -------
            Dictionary with keys:

            - NflxPUreq: Required nitrogen for plant uptake [GtN/yr]
            - NflxPUavail: Total available nitrogen for plant uptake [GtN/yr]
            - NflxPUavail_netMIN: Available nitrogen from net mineralization [GtN/yr]
            - NflxPUdef:  Nitrogen deficit (positive means demand exceeds supply) [GtN/yr]

        Notes
        -----
        The following formulas are used:

        **Nitrogen Required for Plant Uptake:**

        ..  math::

            PUreq = NPPpot \\times NCratio

        where : math:`NPPpot` is the potential NPP (calculated without nitrogen
        limitation) and :math:`NCratio` is the required nitrogen-to-carbon ratio
        (``NCratio_PUreq``). This reflects the requirement for plants to maintain
        specific carbon-to-nitrogen stoichiometry during carbon assimilation.

        **Nitrogen Available for Plant Uptake:**

        Nitrogen availability comes from two sources:  atmospheric deposition and
        net mineralization. The net mineralization available for plant uptake
        (NetMIN) is separated into two components:

        .. math::

            PUavail = AD + NetMIN0 \\times e^{dT \\times s_{netmin2dT}} +
            NetMINlt_0 \\times e^{AD \\times s_{netmin2AD}} \\times
            \\frac{NPPpot}{NPP0}

        where:

        - :math:`AD` is atmospheric nitrogen deposition (assumed immediately available)
        - :math:`NetMIN0` is the baseline short-term net mineralization rate
        - :math:`s_(netmin2dT)` is the temperature sensitivity (``sns_NetMIN2dT``)
        - :math:`dT` is the temperature change
        - :math:`NetMINlt0` is the baseline long-term net mineralization rate
        - :math:`s_{netmin2AD}` is the sensitivity to nitrogen deposition (``sns_NetMINlt2AD``)

        **Components of Net Mineralization:**

        1. **Short-term component** (:math:`NetMIN_0 \\times e^{dT \\times sns_T}`):
           Represents the portion of current net mineralization available for plant
           uptake from fast-turnover organic matter pools. Temperature-dependent
           because decomposition rates increase with warming.

        2. **Long-term component** (:math:`NetMINlt0 \\times e^{AD \\times s_{netmin2AD}}
           \\times \\frac{NPPpot}{NPP0}`): Represents nitrogen that accumulates in
           the mineral nitrogen pool (from the long-term past) and becomes available
           for plant uptake in the current timestep. Scaled by:

           - :math:`NPPpot/NPP0` ratio:  Higher NPP produces more organic matter
             substrate for decomposition/mineralization, controlling the baseline
             rate for long-term net mineralization
           - Exponential response to AD: Incorporates that net mineralization
             (mineralization - immobilization) is nitrogen-demanding, with deposition
             affecting the balance through priming effects (Cheng et al., 2019)

        **Plant Uptake Deficit:**

        ..  math::

            PUdef = PUreq - PUavail

        The plant uptake deficit represents the balance between nitrogen required
        and nitrogen available specifically for the plant uptake process:

        - When :math:`PUdef > 0`: Nitrogen required exceeds availability, indicating
          nitrogen limitation on plant growth
        - When :math:`PUdef < 0`: Nitrogen available exceeds requirement, indicating
          nitrogen surplus

        This deficit is a key diagnostic variable used in nitrogen effect calculation
        to determine compensatory responses in various carbon-nitrogen processes (see
        methods in ``EffectCarbonNitrogenCouplingCalculator``). The nitrogen effect
        on NPP follows an exponential response:  :math:`\\epsilon_{N(NPP)} = e^{s_{
        NPP2PUdef} \\times PUdef}`, where the resulting effect is used to calculate actual NPP
        from potential NPP.

        This formulation does not explicitly represent the full ecosystem mineral
        nitrogen budget, which would include additional fluxes such as microbial
        immobilization, denitrification, and leaching. Net mineralization implicitly
        accounts for the balance between gross mineralization and microbial immobilization,
        while other nitrogen losses are not directly modeled in this simplified framework.

        **Methodological Rationale for the Flux-Based Approximation:**

        The parameterization of PUavail effectively approximates mineral nitrogen
        availability for plant uptake without explicitly relying on the current
        mineral nitrogen pool size. This flux-based approximation is necessary for
        representing nitrogen deficit at an annual timestep without modeling sub-annual
        dynamics of mineral nitrogen processes.

        In complex carbon-nitrogen models, nitrogen availability is typically based on
        the current mineral nitrogen pool size (mass unit), and nitrogen requirement is
        computed from integrated fluxes in a given timestep (mass unit) (Thornton et al.,
        2007; Wiltshire et al., 2021; Zaehle et al., 2014). Competition from microbial
        immobilization is also considered in some complex models. However, in a model
        with an annual timestep like CNit, such a pool-based system would be inherently
        unstable:  the mineral nitrogen pool size would be orders of magnitude smaller
        than the annual nitrogen demand, since the turnover time of the mineral nitrogen
        pool is substantially shorter than the annual timestep.

        While the full ecosystem mineral nitrogen budget is usually balanced over long
        timescales (annual in CNit), short-term imbalances are the key that
        determines nutrient limitation on processes.  To represent this effect, we
        transform the mineral nitrogen availability from a pool size basis to a flux
        basis by estimating the net mineralization fluxes that would be available for
        plant uptake over a year. This flux-based approach allows us to compare nitrogen
        supply and demand on comparable annual timescales, capturing the essential
        dynamics of nitrogen limitation without explicitly modeling the rapid turnover
        of mineral nitrogen pools.
        """
        # Calculate required nitrogen
        NflxPUreq = CflxNPP * self.NCratio_PUreq

        # Calculate available nitrogen from net mineralization
        NflxPUavail_netMIN = self.NetMIN0 * np.exp(
            dT * self.sns_NetMIN2dT
        ) + self.NetMINlt0 * np.exp(NflxAD * self.sns_NetMINlt2AD) * (
                                     CflxNPP / CflxNPP0
                             )

        # Calculate total available nitrogen
        NflxPUavail = NflxPUavail_netMIN + NflxAD

        return {
            "NflxPUreq": NflxPUreq,
            "NflxPUavail": NflxPUavail,
            "NflxPUavail_netMIN": NflxPUavail_netMIN,
            "NflxPUdef": NflxPUreq - NflxPUavail,
        }


@define
class NitrogenTurnoverCalculator:
    """Calculator for nitrogen turnover processes in terrestrial pools.

    This calculator determines nitrogen fluxes between different pools in the
    terrestrial nitrogen cycle based on pool sizes, turnover times, and
    environmental effects.  The turnover processes include:

    - Plant nitrogen turnover (nLP): Transfer from plant to litter pool
    - Litter nitrogen decomposition (nLD): Transfer from litter to soil pool
    - Soil nitrogen turnover (nSR): Release from soil organic matter
    - Mineral nitrogen loss (nLS): Gaseous and leaching losses from mineral pool

    Each flux is controlled by pool size, baseline turnover time, and environmental
    effects including temperature, nitrogen availability, and land use change.
    """

    # Pool turnover times
    tau_NplsP: pint.Quantity = field(
        default=Q(10, "yr"),
        validator=check_units("yr"),
    )
    """Plant nitrogen pool turnover time under reference conditions [yr].  
    
    Determines the baseline rate at which nitrogen is transferred from plant 
    biomass to litter through senescence and mortality. 
    """

    tau_NplsL: pint.Quantity = field(
        default=Q(1, "yr"),
        validator=check_units("yr"),
    )
    """Litter nitrogen pool turnover time under reference conditions [yr]. 
    
    Determines the baseline rate at which nitrogen is released from litter 
    and transferred to soil organic matter through decomposition.
    """

    tau_NplsS: pint.Quantity = field(
        default=Q(100, "yr"),
        validator=check_units("yr"),
    )
    """Soil nitrogen pool turnover time under reference conditions [yr]. 
    
    Determines the baseline rate at which nitrogen is released from soil 
    organic matter through decomposition and mineralization processes.
    """

    tau_NplsM: pint.Quantity = field(
        default=Q(1, "yr"),
        validator=check_units("yr"),
    )
    """Mineral nitrogen pool turnover time under reference conditions [yr]. 
    
    Determines the baseline rate at which mineral nitrogen is lost from the 
    system through gaseous emissions (denitrification, volatilization) and 
    leaching.
    """

    # Fractionation parameters
    frc_nLSgas: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """Fraction of mineral nitrogen loss occurring as gaseous loss [dimensionless]. 
    Partitions total mineral nitrogen loss between gaseous pathways (denitrification, 
    volatilization) and leaching.  Values range from 0 (all leaching) to 1 (all gaseous). 
    Only gaseous losses are temperature-sensitive in this formulation.
    """

    def calculate_nLPLDSR(
            self,
            NplsP: pint.Quantity,
            NplsL: pint.Quantity,
            NplsS: pint.Quantity,
            eff_dT_nLP: pint.Quantity,
            eff_dT_nLD: pint.Quantity,
            eff_dT_nSR: pint.Quantity,
            eff_N_nLP: pint.Quantity,
            eff_N_nLD: pint.Quantity,
            eff_N_nSR: pint.Quantity,
            eff_LU_nLP: pint.Quantity,
            eff_LU_nLD: pint.Quantity,
            eff_LU_nSR: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate nitrogen turnover fluxes from plant, litter, and soil nitrogen
        pools based on pool sizes, turnover times, and environmental effect factors.

        Parameters
        ----------
        NplsP
            Plant nitrogen pool size [GtN]. Total nitrogen contained in living
            plant biomass.
        NplsL
            Litter nitrogen pool size [GtN]. Total nitrogen contained in dead
            plant material (litter).
        NplsS
            Soil nitrogen pool size [GtN]. Total nitrogen contained in soil
            organic matter.
        eff_dT_nLP
            Temperature effect on litter production nitrogen turnover [dimensionless].
        eff_dT_nLD
            Temperature effect on litter decomposition nitrogen turnover [dimensionless].
        eff_dT_nSR
            Temperature effect on soil nitrogen turnover [dimensionless].
        eff_N_nLP
            Carbon-nitrogen coupling effect on plant-to-litter turnover [dimensionless].
        eff_N_nLD
            Carbon-nitrogen coupling effect on litter decomposition [dimensionless].
        eff_N_nSR
            Carbon-nitrogen coupling effect on soil nitrogen turnover [dimensionless].
        eff_LU_nLP
            Land use effect on litter production nitrogen turnover [dimensionless].
        eff_LU_nLD
            Land use effect on litter decomposition nitrogen turnover [dimensionless].
        eff_LU_nSR
            Land use effect on soil respiration nitrogen turnover [dimensionless].

        Returns
        -------
            Dictionary with keys:

            - NflxLP: Nitrogen flux from litter production [GtN/yr]
            - NflxLD: Nitrogen flux from litter decomposition [GtN/yr]
            - NflxSR: Nitrogen flux from soil respiration [GtN/yr]

        Notes
        -----
        The following formulas are used:

        .. math::

            nLP = \\frac{N_P}{\\tau_{N_P}} \\times \\epsilon_{T(nLP)} \\times
            \\epsilon_{N(nLP)} \\times \\epsilon_{LU(nLP)}

            nLD = \\frac{N_L}{\\tau_{N_L}} \\times \\epsilon_{T(nLD)} \\times
            \\epsilon_{N(nLD)} \\times \\epsilon_{LU(nLD)}

            nSR = \\frac{N_S}{\\tau_{N_S}} \\times \\epsilon_{T(nSR)} \\times
            \\epsilon_{N(nSR)} \\times \\epsilon_{LU(nSR)}

        where:

        - :math:`N_P`, :math:`N_L`, :math:`N_S` are the plant, litter, and soil
          nitrogen pool sizes
        - :math:`\\tau_{N_P}`, :math:`\\tau_{N_L}`, :math:`\\tau_{N_S}` are the
          respective turnover times
        - :math:`\\epsilon_{T}` represents temperature effects
        - :math:`\\epsilon_{N}` represents carbon-nitrogen coupling effects
        - : math:`\\epsilon_{LU}` represents land use effects

        Each flux is calculated using first-order kinetics where the turnover rate
        equals pool size divided by turnover time, modified by environmental effects.
        This formulation assumes that turnover processes follow exponential decay
        dynamics at baseline conditions, with multiplicative effects from environmental
        factors.
        """
        return {
            "NflxLP": NplsP / self.tau_NplsP * eff_dT_nLP * eff_N_nLP * eff_LU_nLP,
            "NflxLD": NplsL / self.tau_NplsL * eff_dT_nLD * eff_N_nLD * eff_LU_nLD,
            "NflxSR": NplsS / self.tau_NplsS * eff_dT_nSR * eff_N_nSR * eff_LU_nSR,
        }

    def calculate_nLS(
            self,
            NplsM: pint.Quantity,
            eff_dT_nLSgas: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate nitrogen loss flux from the mineral nitrogen pool.

        This method calculates total nitrogen loss from the mineral pool through
        two pathways: gaseous losses (denitrification and volatilization) and
        leaching.  Only gaseous losses are assumed to be temperature-sensitive.

        Parameters
        ----------
        NplsM
            Mineral nitrogen pool size [GtN].  Total plant-available mineral
            nitrogen in the soil.
        eff_dT_nLSgas
            Temperature effect on gaseous nitrogen loss from mineral pool [dimensionless].
            Values > 1 indicate enhanced gaseous losses with warming.  Only affects
            the gaseous component of total losses.

        Returns
        -------
            Dictionary with keys:

            - NflxLS: Total nitrogen loss from mineral pool [GtN/yr]

        Notes
        -----
        The following formula is used:

        .. math::

            nLS = \\frac{N_M}{\\tau_{N_M}} \\times \\left( f_{gas} \\times
            \\epsilon_{T(nLSgas)} + (1 - f_{gas}) \\right)

        where:

        - :math:`N_M` is the mineral nitrogen pool size
        - :math:`\\tau_{N_M}` is the mineral nitrogen pool turnover time
        - :math:`f_{gas}` is the fraction of losses occurring as gaseous emissions
          (``frc_nLSgas``)
        - :math:`\\epsilon_{T(nLSgas)}` is the temperature effect on gaseous losses
        - :math:`(1 - f_{gas})` represents the leaching fraction, which is assumed
          temperature-insensitive

        The total loss is partitioned between:

        1. **Gaseous losses** (:math:`f_{gas}` fraction): Temperature-sensitive
           losses through denitrification (conversion to N₂, N₂O, NO) and ammonia
           volatilization.  These processes are enhanced by warming.

        2. **Leaching losses** (:math:`1 - f_{gas}` fraction): Temperature-insensitive
           losses through nitrate leaching and dissolved organic nitrogen export.
           These are primarily controlled by hydrology rather than temperature.

        This formulation uses first-order kinetics for baseline losses with differential
        temperature sensitivity for the two loss pathways, reflecting their different
        controlling mechanisms.
        """
        return {
            "NflxLS": (
                    NplsM
                    / self.tau_NplsM
                    * (self.frc_nLSgas * eff_dT_nLSgas + (1 - self.frc_nLSgas))
            ),
        }


@define
class NitrogenCycleCalculator:
    """Calculator for terrestrial nitrogen cycle dynamics.

    This class handles the complete nitrogen cycle by tracking nitrogen transfers
    between pools and calculating pool dynamics. The nitrogen cycle includes four
    main pools:

    - Plant nitrogen pool (NplsP): Nitrogen in living plant biomass
    - Litter nitrogen pool (NplsL): Nitrogen in dead plant material
    - Soil nitrogen pool (NplsS): Nitrogen in soil organic matter
    - Mineral nitrogen pool (NplsM): Plant-available inorganic nitrogen

    Key processes include:

    - External inputs: Biological nitrogen fixation (BNF), atmospheric deposition
      (AD), and fertilizer application (FT)
    - Plant uptake (PU): Transfer from mineral to organic pools
    - Organic matter turnover: Litter production (nLP), litter decomposition (nLD),
      and soil respiration (nSR)
    - Mineralization: Release of mineral nitrogen from organic matter decomposition
    - Losses: Gaseous emissions and leaching (nLS)
    - Land use changes:  Gross deforestation (nLUgrs) and direct mineral loss (nLUmin)

    Each input flux is partitioned among pools using empirically-derived fractions
    that represent the distribution of nitrogen through different pathways in the
    terrestrial biosphere.  These fractions account for the cascade of nitrogen through
    multiple pools within the annual timestep, capturing both immediate allocation
    and subsequent rapid turnover processes.
    """

    # Biological Nitrogen Fixation (BNF) fractions
    frc_BNF2P: pint.Quantity = field(
        default=Q(0.55, "1"),
        validator=check_units("1"),
    )
    """Fraction of biological nitrogen fixation allocated to plant pool [dimensionless]. 
    
    Represents the portion of newly fixed nitrogen directly incorporated into plant 
    biomass through symbiotic associations (e.g., legume-rhizobia symbiosis) or 
    free-living fixers in the rhizosphere. 
    """

    frc_BNF2L: pint.Quantity = field(
        default=Q(0.35, "1"),
        validator=check_units("1"),
    )
    """Fraction of biological nitrogen fixation allocated to litter pool [dimensionless].
     
    Represents newly fixed nitrogen that enters the litter pool within the annual 
    timestep through rapid plant-to-litter turnover (e.g., fine root turnover, leaf 
    senescence of N-fixing plants). The remainder (1 - frc_BNF2P - frc_BNF2L) goes 
    to the soil pool, representing nitrogen that cycles through both plant and litter 
    pools and enters soil within the annual timestep.
    """

    # Plant Uptake (PU) fractions
    frc_PU2P: pint.Quantity = field(
        default=Q(0.55, "1"),
        validator=check_units("1"),
    )
    """Fraction of plant nitrogen uptake allocated to plant pool [dimensionless]. 
    
    Represents nitrogen remaining in plant biomass at the end of the annual timestep, 
    primarily in long-lived tissues (wood, structural tissues, storage organs) that 
    do not turn over within the year. 
    """

    frc_PU2L: pint.Quantity = field(
        default=Q(0.35, "1"),
        validator=check_units("1"),
    )
    """Fraction of plant nitrogen uptake allocated to litter pool [dimensionless]. 
    
    Represents nitrogen that cycles from mineral pool to plant and then to litter 
    within the annual timestep through rapid turnover of short-lived tissues (fine 
    roots, senescing leaves, reproductive structures). The remainder (1 - frc_PU2P - 
    frc_PU2L) goes to the soil pool, representing nitrogen that cascades through 
    plant, litter, and into soil within the annual timestep.
    """

    # Litter fractions
    frc_nLP2L: pint.Quantity = field(
        default=Q(0.7, "1"),
        validator=check_units("1"),
    )
    """Fraction of litter production nitrogen retained in litter pool [dimensionless]. 
    
    Represents nitrogen in plant materials that remain in the litter pool at the end 
    of the annual timestep as recognizable plant structures. The remainder (1 - frc_nLP2L) 
    goes directly to soil, representing materials that decompose rapidly and enter 
    the soil pool within the annual timestep, bypassing extended residence in the 
    litter stage.
    """

    frc_nLD2S: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """Fraction of litter decomposition nitrogen transferred to soil pool [dimensionless].
    
    Represents nitrogen incorporated into stable soil organic matter during decomposition. 
    The remainder (1 - frc_nLD2S) is mineralized to the mineral pool (nLD2M), representing 
    nitrogen released as plant-available forms (NH₄⁺, NO₃⁻) during litter breakdown.
    """

    frc_nLSgas: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """Fraction of mineral nitrogen loss occurring as gaseous emissions [dimensionless]. 
     
    Partitions mineral nitrogen losses between gaseous pathways (denitrification 
    producing N₂, N₂O, NO; ammonia volatilization) and leaching (nitrate leaching, 
    dissolved organic nitrogen export). Values range from 0 (all leaching) to 1 (all 
    gaseous).
    """

    # Grazing fractions
    frc_nLUgrs2P: pint.Quantity = field(
        default=Q(0.6, "1"),
        validator=check_units("1"),
    )
    """Fraction of gross deforestation nitrogen loss from plant pool [dimensionless]. 
    
    Represents the distribution of nitrogen lost during land conversion, with this 
    fraction coming from plant biomass removal or burning. 
    """

    frc_nLUgrs2L: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """Fraction of gross deforestation nitrogen loss from litter pool [dimensionless]. 
    
    Represents nitrogen lost from litter during land conversion. The remainder 
    (1 - frc_nLUgrs2P - frc_nLUgrs2L) comes from the soil pool. 
    """

    def calculate_dNplsPLS_dt(
            self,
            NflxPU: pint.Quantity,
            NflxBNF: pint.Quantity,
            NflxLP: pint.Quantity,
            NflxLD: pint.Quantity,
            NflxSR: pint.Quantity,
            NflxLUgrs: pint.Quantity,
    ) -> List[pint.Quantity]:
        """
        Calculate rate of change for plant, litter, and soil nitrogen pools.

        This method computes the time derivatives of organic nitrogen pools based
        on inputs (BNF, PU), internal transfers (nLP, nLD), outputs (nSR), and land
        use changes (nLUgrs). Each input flux is partitioned among pools according
        to the fractional allocation parameters, which account for nitrogen cascading
        through multiple pools within the annual timestep.

        Parameters
        ----------
        NflxPU
            Plant nitrogen uptake flux [GtN/yr].  Total uptake from mineral pool.
        NflxBNF
            Biological nitrogen fixation flux [GtN/yr]. New nitrogen entering system
            from atmospheric N₂.
        NflxLP
            Litter production flux [GtN/yr]. Nitrogen transfer from plant to litter
            through senescence and mortality.
        NflxLD
            Litter decomposition flux [GtN/yr]. Nitrogen release from litter through
            decomposition.
        NflxSR
            Soil respiration flux [GtN/yr]. Nitrogen release from soil organic matter
            through decomposition (also called soil turnover).
        NflxLUgrs
            Gross deforestation nitrogen flux [GtN/yr].  Nitrogen lost from system
            due to land conversion.

        Returns
        -------
            List containing:

            - dNplsP_dt: Rate of change of plant nitrogen pool [GtN/yr]
            - dNplsL_dt: Rate of change of litter nitrogen pool [GtN/yr]
            - dNplsS_dt: Rate of change of soil nitrogen pool [GtN/yr]

        Notes
        -----
        The pool dynamics follow mass balance equations:

        **Plant Pool:**

        ..  math::

            \\frac{dN_P}{dt} = BNF \\times f_{BNF2P} + PU \\times f_{PU2P} - nLP -
            nLUgrs \\times f_{nLUgrs2P}

        **Litter Pool:**

        .. math::

            \\frac{dN_L}{dt} = BNF \\times f_{BNF2L} + PU \\times f_{PU2L} +
            nLP \\times f_{nLP2L} - nLD - nLUgrs \\times f_{nLUgrs2L}

        **Soil Pool:**

        .. math::

            \\frac{dN_S}{dt} = BNF \\times f_{BNF2S} + PU \\times f_{PU2S} +
            nLP \\times f_{nLP2S} + nLD \\times f_{nLD2S} - nSR - nLUgrs \\times f_{nLUgrs2S}

        where:

        - :math:`f_{BNF2S} = 1 - f_{BNF2P} - f_{BNF2L}` represents BNF nitrogen
          cascading through plant→litter→soil within the annual timestep
        - :math:`f_{PU2S} = 1 - f_{PU2P} - f_{PU2L}` represents uptake nitrogen
          cascading through plant→litter→soil within the annual timestep
        - :math:`f_{nLP2S} = 1 - f_{nLP2L}` represents rapid litter decomposition
          directly entering soil within the annual timestep
        - :math:`f_{nLUgrs2S} = 1 - f_{nLUgrs2P} - f_{nLUgrs2L}` (remainder from soil)

        **Important Note on Annual Timestep:**

        All allocation fractions account for the cascade of nitrogen through multiple
        pools within the annual timestep. For example, :math:`f_{PU2S}` does not
        represent direct allocation from mineral to soil, but rather nitrogen that
        is taken up by plants, then enters litter through turnover, and subsequently
        enters soil through decomposition—all within a single year.  This approach
        effectively captures rapid cycling processes without explicitly modeling
        sub-annual dynamics.
        """
        # Calculate plant uptake partitioning
        NflxPU2P = NflxPU * self.frc_PU2P
        NflxPU2L = NflxPU * self.frc_PU2L
        NflxPU2S = NflxPU - NflxPU2P - NflxPU2L

        # Calculate biological nitrogen fixation partitioning
        NflxBNF2P = NflxBNF * self.frc_BNF2P
        NflxBNF2L = NflxBNF * self.frc_BNF2L
        NflxBNF2S = NflxBNF - NflxBNF2P - NflxBNF2L

        # Calculate grazing uptake partitioning
        NflxLUgrs2P = NflxLUgrs * self.frc_nLUgrs2P
        NflxLUgrs2L = NflxLUgrs * self.frc_nLUgrs2L
        NflxLUgrs2S = NflxLUgrs - NflxLUgrs2P - NflxLUgrs2L

        # Calculate litter production partitioning
        NflxLP2L = NflxLP * self.frc_nLP2L
        NflxLP2S = NflxLP - NflxLP2L

        # Calculate litter decomposition and mineralization
        NflxLD2S = NflxLD * self.frc_nLD2S

        # Calculate changes in nitrogen pools
        dNplsP_dt = NflxBNF2P + NflxPU2P - NflxLP - NflxLUgrs2P
        dNplsL_dt = NflxBNF2L + NflxPU2L + NflxLP2L - NflxLD - NflxLUgrs2L
        dNplsS_dt = NflxBNF2S + NflxPU2S + NflxLP2S + NflxLD2S - NflxSR - NflxLUgrs2S

        return [dNplsP_dt, dNplsL_dt, dNplsS_dt]

    def calculate_dNplsM_dt(
            self,
            NflxPU: pint.Quantity,
            NflxNetMIN: pint.Quantity,
            NflxLS: pint.Quantity,
            NflxAD: pint.Quantity,
            NflxFT: pint.Quantity,
            NflxLUmin: pint.Quantity,
    ) -> List[pint.Quantity]:
        """
        Calculate rate of change for mineral nitrogen pool.

        This method computes the time derivative of the mineral nitrogen pool based
        on inputs (AD, FT, NetMIN), outputs (PU, nLS), and land use changes (nLUmin).

        Parameters
        ----------
        NflxPU
            Plant nitrogen uptake flux [GtN/yr].  Removal of mineral N by plants.
        NflxNetMIN
            Net mineralization flux [GtN/yr]. Release of mineral nitrogen from
            organic matter decomposition (mineralization minus microbial immobilization).
        NflxLS
            Nitrogen loss flux from mineral pool [GtN/yr].  Combined gaseous and
            leaching losses.
        NflxAD
            Atmospheric deposition flux [GtN/yr]. External input from wet and dry
            deposition.
        NflxFT
            Fertilizer application flux [GtN/yr].  Anthropogenic input from agricultural
            fertilization.
        NflxLUmin
            Direct mineral nitrogen loss from land use change [GtN/yr].  Mineral N
            lost during land conversion events.

        Returns
        -------
        List[pint. Quantity]
            List containing [dNplsM_dt] in units of [GtN/yr].

        Notes
        -----
        The mineral pool dynamics follow a mass balance equation:

        .. math::

            \\frac{dN_M}{dt} = AD + FT + NetMIN - PU - nLS - nLUmin

        where:

        - :math:`AD` is atmospheric deposition
        - :math:`FT` is fertilizer application
        - :math:`NetMIN` is net mineralization (gross mineralization - immobilization)
        - :math:`PU` is plant uptake
        - :math:`nLS` is total loss (gaseous + leaching)
        - :math:`nLUmin` is land use-induced loss

        The mineral nitrogen pool represents plant-available inorganic nitrogen
        (primarily NH₄⁺ and NO₃⁻) and typically has rapid turnover relative to
        organic pools. This pool mediates the coupling between organic nitrogen
        cycling and plant nitrogen demand.
        """

        # Calculate changes in nitrogen pools
        dNplsM_dt = NflxAD + NflxFT + NflxNetMIN - NflxPU - NflxLS - NflxLUmin

        return [dNplsM_dt]

    def calculate_Nflx_all(
            self,
            NflxPU: pint.Quantity,
            NflxBNF: pint.Quantity,
            NflxLP: pint.Quantity,
            NflxLD: pint.Quantity,
            NflxSR: pint.Quantity,
            NflxLS: pint.Quantity,
            NflxAD: pint.Quantity,
            NflxFT: pint.Quantity,
            NflxLUgrs: pint.Quantity,
            NflxLUmin: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate all nitrogen fluxes and their partitioning in the system.

        This method computes the complete set of nitrogen fluxes including all
        partitioned flows between pools.  Useful for detailed budget analysis and
        diagnostic output.

        Parameters
        ----------
        NflxPU
            Plant nitrogen uptake flux [GtN/yr].
        NflxBNF
            Biological nitrogen fixation flux [GtN/yr].
        NflxLP
            Litter production flux [GtN/yr].
        NflxLD
            Litter decomposition flux [GtN/yr].
        NflxSR
            Soil respiration flux [GtN/yr] (also called soil turnover).
        NflxLS
            Nitrogen loss flux [GtN/yr].
        NflxAD
            Atmospheric deposition flux [GtN/yr].
        NflxFT
            Fertilizer application flux [GtN/yr].
        NflxLUgrs
            Gross deforestation nitrogen flux [GtN/yr].
        NflxLUmin
            Direct mineral nitrogen loss from land use [GtN/yr].

        Returns
        -------
            Dictionary containing all nitrogen fluxes and their partitioned components:

            - Primary fluxes:  NflxPU, NflxBNF, NflxLP, NflxLD, NflxSR, NflxLS,
              NflxAD, NflxFT, NflxLUgrs, NflxLUmin
            - PU partitioning: NflxPU2P, NflxPU2L, NflxPU2S
            - BNF partitioning: NflxBNF2P, NflxBNF2L, NflxBNF2S
            - LP partitioning: NflxLP2L, NflxLP2S
            - LD partitioning:  NflxLD2S, NflxLD2M
            - LUgrs partitioning: NflxLUgrs2P, NflxLUgrs2L, NflxLUgrs2S
            - Derived flux: NflxNetMIN (total net mineralization = nLD2M + nSR)

        Notes
        -----
        Net mineralization (NetMIN) is calculated as:

        .. math::

            NetMIN = nLD_{2M} + nSR

        where:

        - :math:`nLD_{2M} = nLD \\times (1 - f_{nLD2S})` is litter decomposition to mineral
        - :math:`nSR` is soil respiration (all assumed to mineralize)

        This represents the total flux of nitrogen from organic to mineral forms,
        which is the internal source of plant-available nitrogen in the system.
        """
        NflxPU2P = NflxPU * self.frc_PU2P
        NflxPU2L = NflxPU * self.frc_PU2L
        NflxPU2S = NflxPU - NflxPU2P - NflxPU2L

        NflxBNF2P = NflxBNF * self.frc_BNF2P
        NflxBNF2L = NflxBNF * self.frc_BNF2L
        NflxBNF2S = NflxBNF - NflxBNF2P - NflxBNF2L

        NflxLUgrs2P = NflxLUgrs * self.frc_nLUgrs2P
        NflxLUgrs2L = NflxLUgrs * self.frc_nLUgrs2L
        NflxLUgrs2S = NflxLUgrs - NflxLUgrs2P - NflxLUgrs2L

        NflxLP2L = NflxLP * self.frc_nLP2L
        NflxLP2S = NflxLP - NflxLP2L

        NflxLD2S = NflxLD * self.frc_nLD2S
        NflxLD2M = NflxLD - NflxLD2S

        NflxNetMIN = NflxLD2M + NflxSR  # total nitrogen mineralization flux

        return {
            "NflxPU": NflxPU,
            "NflxPU2P": NflxPU2P,
            "NflxPU2L": NflxPU2L,
            "NflxPU2S": NflxPU2S,
            "NflxBNF": NflxBNF,
            "NflxBNF2P": NflxBNF2P,
            "NflxBNF2L": NflxBNF2L,
            "NflxBNF2S": NflxBNF2S,
            "NflxLP": NflxLP,
            "NflxLP2L": NflxLP2L,
            "NflxLP2S": NflxLP2S,
            "NflxLD": NflxLD,
            "NflxLD2S": NflxLD2S,
            "NflxLD2M": NflxLD2M,
            "NflxSR": NflxSR,
            "NflxNetMIN": NflxNetMIN,
            "NflxLS": NflxLS,
            "NflxLUgrs": NflxLUgrs,
            "NflxLUgrs2P": NflxLUgrs2P,
            "NflxLUgrs2L": NflxLUgrs2L,
            "NflxLUgrs2S": NflxLUgrs2S,
            "NflxLUmin": NflxLUmin,
            "NflxAD": NflxAD,
            "NflxFT": NflxFT,
        }


@define
class NitrogenCycleModelResult:
    """Results from running the nitrogen cycle model.

    This class stores the complete output from a nitrogen cycle model simulation,
    including nitrogen pool states, fluxes, and all environmental effect modifiers.
    Results can be converted to an xarray Dataset for analysis and visualization.

    The nitrogen cycle model tracks four terrestrial nitrogen pools:

    - Plant nitrogen pool (NplsP): Nitrogen in living plant biomass
    - Litter nitrogen pool (NplsL): Nitrogen in dead plant material
    - Soil nitrogen pool (NplsS): Nitrogen in soil organic matter
    - Mineral nitrogen pool (NplsM): Plant-available inorganic nitrogen

    Key processes represented include:

    - External inputs: Biological nitrogen fixation (BNF), atmospheric deposition (AD),
      and fertilizer application (FT)
    - Plant uptake (PU): Transfer from mineral to organic pools
    - Organic matter turnover (nLP, nLD, nSR)
    - Mineral nitrogen losses (nLS): Gaseous emissions and leaching
    - Land use changes:  Gross deforestation (nLUgrs) and mineral losses (nLUmin)
    - Environmental modifiers: temperature, nitrogen, and land use effects

    Note:
        All `_t` attributes are time-continuous functions created through linear
        interpolation of discrete time points.  These functions map from time [yr]
        to their respective units, allowing for smooth temporal interpolation of
        forcing data and effect modifiers.
    """

    n_state: Dict[str, pint.Quantity]
    """Nitrogen pool states over time [GtN].   

    Dictionary containing time series of nitrogen pools:  
    - 'NplsP': Plant nitrogen pool
    - 'NplsL': Litter nitrogen pool
    - 'NplsS':  Soil nitrogen pool
    - 'NplsM': Mineral nitrogen pool
    """

    calc_n_turnover: NitrogenTurnoverCalculator
    """Calculator for nitrogen turnover processes.  

    Instance of NitrogenTurnoverCalculator containing turnover time parameters
    and methods for calculating organic matter breakdown and mineral N loss fluxes 
    (nLP, nLD, nSR, nLS).
    """

    calc_n_cycle: NitrogenCycleCalculator
    """Calculator for nitrogen cycle mass balance. 

    Instance of NitrogenCycleCalculator containing allocation fractions and
    methods for calculating flux partitioning and pool dynamics. 
    """

    time_axis: pint.Quantity
    """Time points for simulation [yr]. 

    Array of discrete time points at which the model state is evaluated,
    typically representing annual timesteps.
    """

    # Time-interpolated flux functions
    NflxPU_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated plant nitrogen uptake flux function [GtN/yr].  

    Continuous function mapping time to plant uptake flux.  Interpolated
    from discrete forcing data to allow smooth temporal variation.
    """

    NflxBNF_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated biological nitrogen fixation flux function [GtN/yr].

    Continuous function mapping time to BNF flux. Represents new nitrogen entering
    the system from atmospheric N₂.  Interpolated from discrete forcing data. 
    """

    NflxAD_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated atmospheric nitrogen deposition flux function [GtN/yr].

    Continuous function mapping time to deposition flux. Represents external nitrogen
    input from wet and dry deposition. Interpolated from discrete forcing data.
    """

    NflxFT_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated fertilizer application flux function [GtN/yr].

    Continuous function mapping time to fertilizer flux.  Represents anthropogenic
    nitrogen input from agricultural fertilization. Interpolated from discrete forcing data.
    """

    NflxLUgrs_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated gross land use change nitrogen flux function [GtN/yr].

    Continuous function mapping time to land use emissions. Represents nitrogen
    lost from organic pools during land conversion (deforestation, agricultural expansion).
    Interpolated from discrete forcing data.
    """

    NflxLUmin_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated land use mineral nitrogen loss flux function [GtN/yr].

    Continuous function mapping time to direct mineral nitrogen losses. Represents
    mineral N lost during land conversion events. Interpolated from discrete forcing data.
    """

    # Time-interpolated effect functions
    eff_dT_nLP_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated temperature effect on litter production function [dimensionless].

    Continuous function mapping time to temperature modifier for plant-to-litter
    nitrogen turnover. Interpolated from calculated temperature effects.
    """

    eff_dT_nLD_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated temperature effect on litter decomposition function [dimensionless].

    Continuous function mapping time to temperature modifier for litter-to-soil
    nitrogen turnover. Interpolated from calculated temperature effects.
    """

    eff_dT_nSR_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated temperature effect on soil respiration function [dimensionless].

    Continuous function mapping time to temperature modifier for soil organic
    nitrogen turnover. Interpolated from calculated temperature effects.
    """

    eff_dT_nLSgas_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated temperature effect on gaseous nitrogen losses function [dimensionless]. 

    Continuous function mapping time to temperature modifier for gaseous emissions
    (denitrification, volatilization) from the mineral nitrogen pool. Interpolated
    from calculated temperature effects.
    """

    eff_N_nLP_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated nitrogen effect on litter production function [dimensionless].

    Continuous function mapping time to nitrogen modifier for plant-to-litter nitrogen
    turnover. Represents carbon-nitrogen coupling effects on plant tissue turnover.
    Interpolated from calculated nitrogen effects.
    """

    eff_N_nLD_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated nitrogen effect on litter decomposition function [dimensionless].

    Continuous function mapping time to nitrogen modifier for litter decomposition.
    Represents how nitrogen availability affects decomposer activity. Interpolated
    from calculated nitrogen effects. 
    """

    eff_N_nSR_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated nitrogen effect on soil respiration function [dimensionless].  

    Continuous function mapping time to nitrogen modifier for soil organic nitrogen
    turnover. Represents carbon-nitrogen coupling in soil processes. Interpolated
    from calculated nitrogen effects.
    """

    eff_LU_nLP_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated land use effect on litter production function [dimensionless].

    Continuous function mapping time to land use modifier for plant-to-litter
    nitrogen turnover. Represents how land use change affects plant mortality rates.
    Interpolated from calculated land use effects.
    """

    eff_LU_nLD_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated land use effect on litter decomposition function [dimensionless].  

    Continuous function mapping time to land use modifier for litter decomposition.
    Represents how land use change affects litter breakdown rates. Interpolated
    from calculated land use effects.
    """

    eff_LU_nSR_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated land use effect on soil respiration function [dimensionless]. 

    Continuous function mapping time to land use modifier for soil organic nitrogen
    turnover. Represents how land use change affects soil nitrogen dynamics. 
    Interpolated from calculated land use effects.
    """

    def add_non_state_variables(self) -> xr.Dataset:
        """
        Calculate derived variables and create xarray Dataset with complete results.

        This method computes all diagnostic fluxes and derived quantities from the
        model state variables and returns a complete xarray Dataset containing:

        - State variables:  Nitrogen pools (NplsP, NplsL, NplsS, NplsM)
        - Turnover fluxes: Calculated from pools and environmental effects (nLP, nLD, nSR, nLS)
        - All partitioned fluxes: PU, BNF, and land use change distributed among pools
        - Derived quantities:  Total organic nitrogen (NplsPLS), total terrestrial nitrogen
          (NplsPLSM), net terrestrial nitrogen flux (NflxNetPLSM)
        - Net mineralization: Total nitrogen release from organic matter (NflxNetMIN)

        Returns
        -------
        xr.Dataset
            Complete model results as xarray Dataset with time dimension and all
            variables labeled with descriptive names and units.  The Dataset includes:

            **State Variables:**

            - NplsP, NplsL, NplsS, NplsM: Individual pool sizes [GtN]
            - NplsPLS: Total organic nitrogen (plant + litter + soil) [GtN]
            - NplsPLSM: Total terrestrial nitrogen (organic + mineral) [GtN]

            **Primary Fluxes:**

            - NflxPU: Plant nitrogen uptake [GtN/yr]
            - NflxBNF:  Biological nitrogen fixation [GtN/yr]
            - NflxAD: Atmospheric deposition [GtN/yr]
            - NflxFT: Fertilizer application [GtN/yr]
            - NflxLP:  Litter production [GtN/yr]
            - NflxLD:  Litter decomposition [GtN/yr]
            - NflxSR: Soil respiration (soil turnover) [GtN/yr]
            - NflxLS: Mineral nitrogen losses [GtN/yr]
            - NflxLUgrs: Gross land use emissions [GtN/yr]
            - NflxLUmin: Land use mineral nitrogen losses [GtN/yr]

            **Partitioned Fluxes:**

            - NflxPU2P, NflxPU2L, NflxPU2S:  Plant uptake allocation to each pool
            - NflxBNF2P, NflxBNF2L, NflxBNF2S: BNF allocation to each pool
            - NflxLP2L, NflxLP2S:  Litter production allocation
            - NflxLD2S, NflxLD2M: Litter decomposition to soil and mineral
            - NflxLUgrs2P, NflxLUgrs2L, NflxLUgrs2S: Land use emissions from each pool

            **Derived Fluxes:**

            - NflxNetMIN: Total net mineralization (nLD2M + nSR) [GtN/yr]
            - NflxNetPLSM: Net change in total terrestrial nitrogen [GtN/yr]

        Notes
        -----
        The net terrestrial nitrogen flux is calculated using numerical differentiation:

        .. math::

            nNetPLSM = \\frac{dN_{PLSM}}{dt} = \\frac{d(N_P + N_L + N_S + N_M)}{dt}

        This represents the net nitrogen balance of the terrestrial biosphere, with
        positive values indicating nitrogen accumulation and negative values indicating
        nitrogen loss.

        The net mineralization is calculated as:

        .. math::

            NetMIN = nLD2M + nSR

        where:

        - :math:`nLD2M = nLD \\times (1 - f_{nLD2S})` is litter decomposition to mineral
        - :math:`nSR` is soil respiration (all assumed to mineralize)

        This represents the total flux of nitrogen from organic to mineral forms,
        which is the internal source of plant-available nitrogen in the system.
        """
        state = self.n_state
        time = self.time_axis

        # Calculate turnover
        turnover = self.calc_n_turnover.calculate_nLPLDSR(
            NplsP=state["NplsP"],
            NplsL=state["NplsL"],
            NplsS=state["NplsS"],
            eff_dT_nLP=self.eff_dT_nLP_t(time),
            eff_dT_nLD=self.eff_dT_nLD_t(time),
            eff_dT_nSR=self.eff_dT_nSR_t(time),
            eff_N_nLP=self.eff_N_nLP_t(time),
            eff_N_nLD=self.eff_N_nLD_t(time),
            eff_N_nSR=self.eff_N_nSR_t(time),
            eff_LU_nLP=self.eff_LU_nLP_t(time),
            eff_LU_nLD=self.eff_LU_nLD_t(time),
            eff_LU_nSR=self.eff_LU_nSR_t(time),
        ) | self.calc_n_turnover.calculate_nLS(
            NplsM=state["NplsM"],
            eff_dT_nLSgas=self.eff_dT_nLSgas_t(time),
        )

        # Calculate all fluxes
        flux = self.calc_n_cycle.calculate_Nflx_all(
            NflxPU=self.NflxPU_t(time),
            NflxBNF=self.NflxBNF_t(time),
            NflxAD=self.NflxAD_t(time),
            NflxFT=self.NflxFT_t(time),
            NflxLUgrs=self.NflxLUgrs_t(time),
            NflxLUmin=self.NflxLUmin_t(time),
            NflxLP=turnover["NflxLP"],
            NflxLD=turnover["NflxLD"],
            NflxSR=turnover["NflxSR"],
            NflxLS=turnover["NflxLS"],
        )

        # Calculate combined pools
        extra_var = {"NplsPLS": state["NplsP"] + state["NplsL"] + state["NplsS"]}
        extra_var["NplsPLSM"] = extra_var["NplsPLS"] + state["NplsM"]
        extra_var["NflxNetPLSM"] = Q(np.gradient(extra_var["NplsPLSM"].m), "GtN/yr")

        return make_dataset_from_var_dict({**state, **flux, **extra_var}, time)


@define
class NitrogenCycleModel:
    """MAGICC's terrestrial nitrogen cycle model (nitrogen component of CNit).

    This model simulates the dynamics of four terrestrial nitrogen pools in the
    terrestrial biosphere by solving a system of ordinary differential equations (ODEs)
    using scipy.integrate.solve_ivp. The model tracks nitrogen transfers between pools
    through processes including:

    - External inputs:  Biological nitrogen fixation (BNF), atmospheric deposition (AD),
      and fertilizer application (FT)
    - Plant uptake (PU): Transfer from mineral to organic pools
    - Organic matter turnover (litter production, litter decomposition, soil respiration)
    - Mineralization: Release of plant-available nitrogen from organic matter
    - Mineral nitrogen losses: Gaseous emissions and leaching
    - Land use changes: Gross deforestation and direct mineral losses
    - Environmental modifiers (temperature, nitrogen availability, land use effects)

    The model uses annual timesteps and accounts for rapid nitrogen cycling through
    empirical allocation fractions that represent nitrogen cascading through multiple
    pools within a year.

    The nitrogen cycle is solved in two stages:
    (1) organic nitrogen pools (plant, litter, soil) are solved first, then
    (2) the mineral nitrogen pool is solved using net mineralization calculated from
    the organic pool dynamics.

    Time-continuous functions (suffixed with `_t`) are generated through linear
    interpolation of discrete time points and provide continuous values over the
    simulation period, allowing the ODE solver to evaluate fluxes at any time point.
    """

    calc_n_turnover: NitrogenTurnoverCalculator
    """Calculator for nitrogen turnover processes. 

    Instance of NitrogenTurnoverCalculator containing turnover time parameters
    (tau_NplsP, tau_NplsL, tau_NplsS, tau_NplsM) and methods for calculating
    organic matter breakdown and mineral N loss fluxes based on pool sizes and
    environmental effects.
    """

    calc_n_cycle: NitrogenCycleCalculator
    """Calculator for nitrogen cycle mass balance.

    Instance of NitrogenCycleCalculator containing allocation fractions
    (frc_BNF2P, frc_BNF2L, frc_PU2P, frc_PU2L, frc_nLP2L, frc_nLD2S,
    frc_nLSgas, frc_nLUgrs2P, frc_nLUgrs2L) and methods for calculating
    flux partitioning and pool dynamics. 
    """

    NplsP0: pint.Quantity = field(
        default=Q(10, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial plant nitrogen pool size [GtN].

    Nitrogen contained in living plant biomass at the start of the simulation
    (time0). Typical values range from 3-6 GtN for global simulations.
    """

    NplsL0: pint.Quantity = field(
        default=Q(1, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial litter nitrogen pool size [GtN]. 

    Nitrogen contained in dead plant material (litter) at the start of the simulation
    (time0). Typical values range from 0.5-2 GtN for global simulations.
    """

    NplsS0: pint.Quantity = field(
        default=Q(100, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial soil nitrogen pool size [GtN].

    Nitrogen contained in soil organic matter at the start of the simulation (time0).
    Typical values range from 95-140 GtN for global simulations.
    """

    NplsM0: pint.Quantity = field(
        default=Q(1, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial mineral nitrogen pool size [GtN].

    Plant-available inorganic nitrogen (primarily NH₄⁺ and NO₃⁻) at the start of
    the simulation (time0). Typical values range from 0.01-0.1 GtN for global
    simulations.  This pool has much faster turnover than organic pools.
    """

    time0: pint.Quantity = field(
        default=Q(1850, "yr"),
        validator=check_units("yr"),
    )
    """Initialization time [yr].

    Time point at which the initial pool sizes (NplsP0, NplsL0, NplsS0, NplsM0) apply.
    The simulation time_axis must start at or after this time.  Typically set to
    the pre-industrial baseline year (1850).
    """

    switch_Npls: list[int] = field(
        default=[1, 1, 1, 1],
    )
    """Nitrogen pool switches for [plant, litter, soil, mineral] [dimensionless].

    Binary switches (1=enabled, 0=disabled) to selectively enable or disable
    pool dynamics.  Useful for diagnostic purposes or simplified model configurations.
    Setting a switch to 0 freezes that pool at its initial value.

    - switch_Npls[0]: Plant pool switch
    - switch_Npls[1]: Litter pool switch
    - switch_Npls[2]: Soil pool switch
    - switch_Npls[3]: Mineral pool switch
    """


    def run(
            self,
            time_axis: pint.Quantity,
            NflxPU_t: Callable[[pint.Quantity], pint.Quantity],
            NflxBNF_t: Callable[[pint.Quantity], pint.Quantity],
            NflxAD_t: Callable[[pint.Quantity], pint.Quantity],
            NflxFT_t: Callable[[pint.Quantity], pint.Quantity],
            NflxLUgrs_t: Callable[[pint.Quantity], pint.Quantity],
            NflxLUmin_t: Callable[[pint.Quantity], pint.Quantity],
            eff_dT_nLP_t: Callable[[pint.Quantity], pint.Quantity],
            eff_dT_nLD_t: Callable[[pint.Quantity], pint.Quantity],
            eff_dT_nSR_t: Callable[[pint.Quantity], pint.Quantity],
            eff_dT_nLSgas_t: Callable[[pint.Quantity], pint.Quantity],
            eff_N_nLP_t: Callable[[pint.Quantity], pint.Quantity],
            eff_N_nLD_t: Callable[[pint.Quantity], pint.Quantity],
            eff_N_nSR_t: Callable[[pint.Quantity], pint.Quantity],
            eff_LU_nLP_t: Callable[[pint.Quantity], pint.Quantity],
            eff_LU_nLD_t: Callable[[pint.Quantity], pint.Quantity],
            eff_LU_nSR_t: Callable[[pint.Quantity], pint.Quantity],
    ) -> NitrogenCycleModelResult:
        """
            Run the nitrogen cycle model simulation.

            This method solves the nitrogen cycle in two stages using scipy.integrate.solve_ivp:
            (1) Solves organic nitrogen pools (plant, litter, soil) first, then
            (2) Solves mineral nitrogen pool using net mineralization from stage 1.

            The solver uses time-continuous interpolated functions for all forcing data
            and environmental effects, allowing flexible time step selection during integration.

            Parameters
            ----------
            time_axis
                Time points for simulation [yr].  Must start at or after time0.
                These are the time points at which solution is explicitly evaluated.
            NflxPU_t
                Time-interpolated plant nitrogen uptake flux function [GtN/yr].
                Maps time to uptake from mineral to organic pools.
            NflxBNF_t
                Time-interpolated biological nitrogen fixation flux function [GtN/yr].
                Maps time to new nitrogen entering from atmospheric N₂.
            NflxAD_t
                Time-interpolated atmospheric deposition flux function [GtN/yr].
                Maps time to external nitrogen input from wet and dry deposition.
            NflxFT_t
                Time-interpolated fertilizer application flux function [GtN/yr].
                Maps time to anthropogenic nitrogen input from agriculture.
            NflxLUgrs_t
                Time-interpolated gross land use change flux function [GtN/yr].
                Maps time to nitrogen lost from organic pools during land conversion.
            NflxLUmin_t
                Time-interpolated land use mineral loss flux function [GtN/yr].
                Maps time to direct mineral nitrogen losses during land conversion.
            eff_dT_nLP_t
                Time-interpolated temperature effect on litter production function [dimensionless].
                Maps time to temperature modifier for plant-to-litter turnover.
            eff_dT_nLD_t
                Time-interpolated temperature effect on litter decomposition function [dimensionless].
                Maps time to temperature modifier for litter-to-soil turnover.
            eff_dT_nSR_t
                Time-interpolated temperature effect on soil respiration function [dimensionless].
                Maps time to temperature modifier for soil organic nitrogen turnover.
            eff_dT_nLSgas_t
                Time-interpolated temperature effect on gaseous losses function [dimensionless].
                Maps time to temperature modifier for denitrification and volatilization.
            eff_N_nLP_t
                Time-interpolated nitrogen effect on litter production function [dimensionless].
                Maps time to nitrogen availability modifier for plant-to-litter turnover.
            eff_N_nLD_t
                Time-interpolated nitrogen effect on litter decomposition function [dimensionless].
                Maps time to nitrogen availability modifier for litter decomposition.
            eff_N_nSR_t
                Time-interpolated nitrogen effect on soil respiration function [dimensionless].
                Maps time to nitrogen availability modifier for soil nitrogen turnover.
            eff_LU_nLP_t
                Time-interpolated land use effect on litter production function [dimensionless].
                Maps time to land use modifier for plant-to-litter turnover.
            eff_LU_nLD_t
                Time-interpolated land use effect on litter decomposition function [dimensionless].
                Maps time to land use modifier for litter decomposition.
            eff_LU_nSR_t
                Time-interpolated land use effect on soil respiration function [dimensionless].
                Maps time to land use modifier for soil nitrogen turnover.

            Returns
            -------
            NitrogenCycleModelResult
                Simulation results containing:

                - Nitrogen pool time series (NplsP, NplsL, NplsS, NplsM)
                - References to calculator instances
                - All input forcing and effect functions
                - Methods to calculate derived fluxes and diagnostics

            Raises
            ------
            ValueError
                If time_axis starts before time0. The simulation cannot run backwards
                from the initialization time.
            SolveError
                If the numerical ODE solver fails to converge.  This can occur due to
                numerical instability, runaway feedbacks, or inappropriate tolerance settings.

            Notes
            -----
            **Two-Stage Solution:**

            The nitrogen cycle is solved in two stages because the mineral pool depends
            on net mineralization from organic pools:

            1. **Stage 1 - Organic pools (plant, litter, soil):**

            ..  math::

                \\frac{dN_P}{dt} = BNF \\times f_{BNF2P} + PU \\times f_{PU2P} - nLP -
                nLUgrs \\times f_{nLUgrs2P}

                \\frac{dN_L}{dt} = BNF \\times f_{BNF2L} + PU \\times f_{PU2L} +
                nLP \\times f_{nLP2L} - nLD - nLUgrs \\times f_{nLUgrs2L}

                \\frac{dN_S}{dt} = BNF \\times f_{BNF2S} + PU \\times f_{PU2S} +
                nLP \\times f_{nLP2S} + nLD \\times f_{nLD2S} - nSR - nLUgrs \\times f_{nLUgrs2S}

            2. **Stage 2 - Mineral pool:**

            After solving organic pools, net mineralization is calculated:

            .. math::

                NetMIN = nLD_{2M} + nSR

            where :math:`nLD_{2M} = nLD \\times (1 - f_{nLD2S})` is the mineralized
            fraction of litter decomposition.

            Then the mineral pool is solved:

            .. math::

                \\frac{dN_M}{dt} = AD + FT + NetMIN - PU - nLS - nLUmin

            where turnover fluxes (nLP, nLD, nSR, nLS) are calculated from pool sizes
            and environmental effects using first-order kinetics.

            **Numerical Methods:**

            Both stages use scipy.integrate.solve_ivp with:

            - Adaptive time stepping (solver selects appropriate substeps)
            - Absolute tolerance:  1e-6 GtN
            - Relative tolerance: 1e-3
            - Default solver method (typically RK45)
            - Stage 2 uses max_step=2 yr for stability due to fast mineral pool turnover

            **Pool Switches:**

            Pool dynamics can be selectively disabled using switch_Npls.  When a switch
            is set to 0, the corresponding pool's rate of change is forced to zero,
            effectively freezing that pool at its initial value throughout the simulation.
            If the mineral pool switch is 0, the mineral pool is set to zero throughout.

            **Time-Continuous Functions:**

            All `_t` arguments are callable functions created through linear interpolation
            of discrete time points.  This allows the ODE solver to evaluate forcing data
            and environmental effects at any time point during integration, not just at
            the discrete output times specified in time_axis.
            """

        def func_to_solve_NplsPLS(t: float, y: np.ndarray) -> list:
            """
            ODE system function for organic nitrogen pools (stage 1).

            Parameters
            ----------
            t
                Current time [yr] (unitless for solver)
            y
                Current state vector [NplsP, NplsL, NplsS] in GtN (unitless for solver)

            Returns
            -------
            list
                Time derivatives [dNplsP/dt, dNplsL/dt, dNplsS/dt] in GtN/yr (unitless)
            """
            t = Q(t, "yr")
            # Calculate turnover fluxes
            turnover = self.calc_n_turnover.calculate_nLPLDSR(
                NplsP=Q(y[0], "GtN"),
                NplsL=Q(y[1], "GtN"),
                NplsS=Q(y[2], "GtN"),
                eff_dT_nLP=eff_dT_nLP_t(t),
                eff_dT_nLD=eff_dT_nLD_t(t),
                eff_dT_nSR=eff_dT_nSR_t(t),
                eff_N_nLP=eff_N_nLP_t(t),
                eff_N_nLD=eff_N_nLD_t(t),
                eff_N_nSR=eff_N_nSR_t(t),
                eff_LU_nLP=eff_LU_nLP_t(t),
                eff_LU_nLD=eff_LU_nLD_t(t),
                eff_LU_nSR=eff_LU_nSR_t(t),
            )

            # Calculate rate of change
            dydt = self.calc_n_cycle.calculate_dNplsPLS_dt(
                NflxPU=NflxPU_t(t),
                NflxBNF=NflxBNF_t(t),
                NflxLUgrs=NflxLUgrs_t(t),
                NflxLP=turnover["NflxLP"],
                NflxLD=turnover["NflxLD"],
                NflxSR=turnover["NflxSR"],
            )
            dydt = [
                v.to("GtN/yr").m * switch
                for v, switch in zip(dydt, self.switch_Npls[:-1])
            ]
            if np.all(np.abs(dydt) < 1e-12):
                dydt = [0, 0, 0]
            return dydt

        def func_to_solve_NplsM(t: float, y: np.ndarray) -> list:
            """
            ODE system function for mineral nitrogen pool (stage 2).

            Parameters
            ----------
            t
                Current time [yr] (unitless for solver)
            y
                Current state vector [NplsM] in GtN (unitless for solver)

            Returns
            -------
            list
                Time derivative [dNplsM/dt] in GtN/yr (unitless)
            """
            t = Q(t, "yr")
            # Calculate turnover fluxes
            turnover = self.calc_n_turnover.calculate_nLS(
                NplsM=Q(y[0], "GtN"),
                eff_dT_nLSgas=eff_dT_nLSgas_t(t),
            )

            # Calculate rate of change
            dydt = self.calc_n_cycle.calculate_dNplsM_dt(
                NflxPU=NflxPU_t(t),
                NflxNetMIN=NflxNetMIN_t(t),
                NflxAD=NflxAD_t(t),
                NflxFT=NflxFT_t(t),
                NflxLUmin=NflxLUmin_t(t),
                NflxLS=turnover["NflxLS"],
            )
            dydt = [dydt[0].to("GtN/yr").m]
            if np.all(np.abs(dydt) < 1e-12):
                dydt = [0]
            return dydt

        if time_axis[0] < self.time0:
            raise ValueError(
                f"time_axis starts before time0: {time_axis[0]} < {self.time0}"
            )

        # Prepare solver arguments
        t_eval = time_axis.to("yr").m
        t_span = (self.time0.to("yr").m, t_eval[-1])

        raw = scipy.integrate.solve_ivp(
            func_to_solve_NplsPLS,
            t_span=t_span,
            t_eval=t_eval,
            y0=(
                self.NplsP0.to("GtN").m,
                self.NplsL0.to("GtN").m,
                self.NplsS0.to("GtN").m,
            ),
            atol=1e-6,
            rtol=1e-3,
            # max_step=2,
            # method='LSODA',
        )

        if not raw.success:
            raise SolveError(
                f"Model failed to solve (possible runaway feedback)\n{raw}"
            )
        res_NplsPLS = {
            "NplsP": Q(raw.y[0, :], "GtN"),
            "NplsL": Q(raw.y[1, :], "GtN"),
            "NplsS": Q(raw.y[2, :], "GtN"),
        }

        if self.switch_Npls[-1] != 0:
            NflxLD2M = (
                    res_NplsPLS["NplsL"]
                    / self.calc_n_turnover.tau_NplsL
                    * eff_dT_nLD_t(time_axis)
                    * eff_N_nLD_t(time_axis)
                    * eff_LU_nLD_t(time_axis)
                    * (1 - self.calc_n_cycle.frc_nLD2S)
            )
            NflxSR = (
                    res_NplsPLS["NplsS"]
                    / self.calc_n_turnover.tau_NplsS
                    * eff_dT_nSR_t(time_axis)
                    * eff_N_nSR_t(time_axis)
                    * eff_LU_nSR_t(time_axis)
            )
            NflxNetMIN = NflxLD2M + NflxSR
            NflxNetMIN_t = _get_interp(NflxNetMIN, time_axis)

            raw = scipy.integrate.solve_ivp(
                func_to_solve_NplsM,
                t_span=t_span,
                t_eval=t_eval,
                y0=(self.NplsM0.to("GtN").m,),
                atol=1e-6,
                rtol=1e-3,
                max_step=2,
                # method='LSODA',
            )

            if not raw.success:
                raise SolveError(
                    f"Model failed to solve (possible runaway feedback)\n{raw}"
                )
            res_NplsM = {"NplsM": Q(raw.y[0, :], "GtN")}
        else:
            res_NplsM = {
                "NplsM": Q(np.zeros_like(time_axis), "GtN"),
            }

        return NitrogenCycleModelResult(
            n_state=res_NplsPLS | res_NplsM,
            calc_n_turnover=self.calc_n_turnover,
            calc_n_cycle=self.calc_n_cycle,
            time_axis=time_axis,
            NflxPU_t=NflxPU_t,
            NflxBNF_t=NflxBNF_t,
            NflxAD_t=NflxAD_t,
            NflxFT_t=NflxFT_t,
            NflxLUgrs_t=NflxLUgrs_t,
            NflxLUmin_t=NflxLUmin_t,
            eff_dT_nLP_t=eff_dT_nLP_t,
            eff_dT_nLD_t=eff_dT_nLD_t,
            eff_dT_nSR_t=eff_dT_nSR_t,
            eff_dT_nLSgas_t=eff_dT_nLSgas_t,
            eff_N_nLP_t=eff_N_nLP_t,
            eff_N_nLD_t=eff_N_nLD_t,
            eff_N_nSR_t=eff_N_nSR_t,
            eff_LU_nLP_t=eff_LU_nLP_t,
            eff_LU_nLD_t=eff_LU_nLD_t,
            eff_LU_nSR_t=eff_LU_nSR_t,
        )