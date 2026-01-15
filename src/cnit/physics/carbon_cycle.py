"""
MAGICC Terrestrial Carbon Cycle Model Implementation

This module provides the implementation of the carbon cycle component of MAGICC's
coupled carbon-nitrogen cycle model (CNit). It simulates the dynamics of terrestrial
carbon pools and their responses to climate change, CO2 fertilization, nitrogen
availability, and land use change.

**Carbon Cycle Components:**

- :class:`CarbonNPPLPRCalculator`: Calculates net primary production (NPP) and litter
  production respiration (LPR) fluxes with environmental modifiers
- :class:`CarbonTurnoverCalculator`: Computes organic matter breakdown fluxes (litter
  production, litter decomposition, soil respiration) based on pool sizes and turnover times
- :class:`CarbonCycleCalculator`: Handles carbon mass balance and flux partitioning among
  plant, litter, and soil pools using empirical allocation fractions
- :class:`CarbonCycleModel`: Solves the carbon cycle ODEs using scipy.integrate.solve_ivp
- :class:`CarbonCycleModelResult`: Stores and processes carbon cycle simulation results

**Key Model Features:**

1. **Three-Pool Structure**: The model tracks carbon in plant, litter, and soil organic
   matter pools, representing the major terrestrial carbon reservoirs.

2. **Annual Timestep with Cascading**: The model operates on annual timesteps but accounts
   for rapid sub-annual cycling through empirical allocation fractions that represent carbon
   cascading through multiple pools within a year.

3. **Environmental Modifiers**: All process rates are modified by dimensionless effect factors
   representing:

   - Temperature effects on decomposition and respiration
   - CO2 fertilization effects on NPP
   - Nitrogen availability effects (carbon-nitrogen coupling)
   - Land use change effects on turnover rates

   (see :mod:`effect_factor_calculators`)

4. **First-Order Kinetics**: Carbon turnover from each pool follows first-order kinetics
   (flux = pool size / turnover time), modified by environmental effects.

5. **Flux Partitioning**: Input fluxes (NPP, land use) are partitioned among pools using
   empirical fractions that capture both direct allocation and rapid within-year cycling.

6. **Carbon-Nitrogen Coupling**:  Nitrogen availability affects carbon cycle processes through
   limitation effects on NPP and decomposition rates (calculated via nitrogen cycle effects).

**Model Structure:**

The carbon cycle includes three terrestrial pools:

- **Plant pool (CplsP)**: Carbon in living plant biomass (wood, leaves, roots, storage)
- **Litter pool (CplsL)**: Carbon in dead plant material undergoing decomposition
- **Soil pool (CplsS)**: Carbon in soil organic matter (humus, stable organic compounds)

Key processes:

- **Carbon fixation**: Net primary production (NPP) converting atmospheric CO₂ to organic carbon
- **Autotrophic respiration**: Litter production respiration (LPR) from fast turnover of plant tissues
- **Organic matter turnover**: Litter production (cLP), litter decomposition (cLD),
  soil respiration (cSR)
- **Heterotrophic respiration**: Total CO₂ release from decomposition (RH = LPR + cLD2A + cSR)
- **Land use emissions**: Carbon release from deforestation and land conversion (cLUgrs)
"""

from typing import Callable, Dict, List
from attrs import define, field
import numpy as np
import scipy
import xarray as xr
import pint

from ..utils.units import check_units, Q
from ..utils.maths import SolveError
from ..utils.data_utils import make_dataset_from_var_dict


@define
class CarbonNPPLPRCalculator:
    """
    Calculator for Net Primary Production and Litter Production Respiration carbon fluxes
    """

    NPP0: pint.Quantity = field(
        default=Q(50, "GtC/yr"),
        validator=check_units("GtC/yr"),
    )
    """
        Initial net primary production (NPP), base NPP without any effects [GtC/yr]
    """

    LPR0: pint.Quantity = field(
        default=Q(0, "GtC/yr"),
        validator=check_units("GtC/yr"),
    )
    """
        Initial Litter production respiration (LPR), base LPR without any effects [GtC/yr]
    """

    def calculate(
            self,
            CflxLUrgr: pint.Quantity,
            eff_CO2_NPP: pint.Quantity,
            eff_dT_NPP: pint.Quantity,
            eff_N_NPP: pint.Quantity,
            eff_LU_NPP: pint.Quantity,
            eff_dT_LPR: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate carbon fluxes for Net Primary Production and Litter Production Respiration

        Parameters
        ----------
        eff_CO2_NPP
            CO2 effect on NPP [dimensionless]
        eff_dT_NPP
            Temperature effect on NPP [dimensionless]
        eff_dT_LPR
            Temperature effect on LPR [dimensionless]
        eff_N_NPP
            Carbon-nitrogen coupling effect on NPP [dimensionless]
        eff_LU_NPP
            Land use effect on NPP [dimensionless]
        CflxLUrgr
            Carbon flux from land-use regrowth [GtC/yr]

        Returns
        -------
            Dictionary with keys:

            - CflxNPP: Net Primary Production carbon flux [GtC/yr]
            - CflxLPR: Litter Production Respiration carbon flux [GtC/yr]

        Notes
        -----
        The following formulas are used:

        .. math::

            NPP = (NPP_0 + LUrgr) \\times \\epsilon_{T(NPP)} \\times \\epsilon_{CO2}
            \\times \\epsilon_{N(NPP)} \\times \\epsilon_{LU(NPP)}

            LPR = LPR0 \\times \\epsilon_{T(LPR)} \\times \\epsilon_{CO2} \\times
            \\epsilon_{N(LPR)} \\times \\epsilon_{LU(NPP)}
        """
        return {
            "CflxNPP": (self.NPP0 + CflxLUrgr)
                       * eff_dT_NPP
                       * eff_CO2_NPP
                       * eff_N_NPP
                       * eff_LU_NPP,
            "CflxLPR": self.LPR0 * eff_dT_LPR * eff_CO2_NPP * eff_N_NPP * eff_LU_NPP,
        }


@define
class CarbonTurnoverCalculator:
    """Calculator for carbon turnover processes in terrestrial pools.

    This calculator determines carbon fluxes between different pools in the
    terrestrial carbon cycle based on pool sizes, turnover times, and
    environmental effects.  The turnover processes include:

    - Litter production (cLP): Transfer from plant to litter pool
    - Litter decomposition (cLD): Transfer from litter to soil pool
    - Soil respiration (cSR): Release from soil organic matter to atmosphere

    Each flux is controlled by pool size, baseline turnover time, and environmental
    effects including temperature, nitrogen availability, and land use change.
    """

    tau_CplsP: pint.Quantity = field(
        default=Q(10, "yr"), validator=check_units("yr")
    )
    """Plant carbon pool turnover time under reference conditions [yr]. 
    
    Determines the baseline rate at which carbon is transferred from plant
    biomass to litter through senescence and mortality.
    """

    tau_CplsL: pint.Quantity = field(
        default=Q(1, "yr"), validator=check_units("yr")
    )
    """Litter carbon pool turnover time under reference conditions [yr].
    
    Determines the baseline rate at which carbon is released from litter
    and transferred to soil organic matter through decomposition.
    """

    tau_CplsS: pint.Quantity = field(
        default=Q(100, "yr"), validator=check_units("yr")
    )
    """Soil carbon pool turnover time under reference conditions [yr].
    
    Determines the baseline rate at which carbon is released from soil
    organic matter through decomposition and respiration processes.
    """

    def calculate(
            self,
            CplsP: pint.Quantity,
            CplsL: pint.Quantity,
            CplsS: pint.Quantity,
            eff_dT_cLP: pint.Quantity,
            eff_dT_cLD: pint.Quantity,
            eff_dT_cSR: pint.Quantity,
            eff_N_cLP: pint.Quantity,
            eff_N_cLD: pint.Quantity,
            eff_N_cSR: pint.Quantity,
            eff_LU_cLP: pint.Quantity,
            eff_LU_cLD: pint.Quantity,
            eff_LU_cSR: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate carbon turnover fluxes from plant, litter, and soil carbon
        pools based on pool sizes, turnover times, and environmental effect factors.

        Parameters
        ----------
        CplsP
            Plant carbon pool size [GtC].  Total carbon contained in living
            plant biomass.
        CplsL
            Litter carbon pool size [GtC]. Total carbon contained in dead
            plant material (litter).
        CplsS
            Soil carbon pool size [GtC]. Total carbon contained in soil
            organic matter.
        eff_dT_cLP
            Temperature effect on litter production carbon turnover [dimensionless].
        eff_dT_cLD
            Temperature effect on litter decomposition carbon turnover [dimensionless].
        eff_dT_cSR
            Temperature effect on soil respiration carbon turnover [dimensionless].
        eff_N_cLP
            Carbon-nitrogen coupling effect on litter production [dimensionless].
        eff_N_cLD
            Carbon-nitrogen coupling effect on litter decomposition [dimensionless].
        eff_N_cSR
            Carbon-nitrogen coupling effect on soil respiration [dimensionless].
        eff_LU_cLP
            Land use effect on litter production carbon turnover [dimensionless].
        eff_LU_cLD
            Land use effect on litter decomposition carbon turnover [dimensionless].
        eff_LU_cSR
            Land use effect on soil respiration carbon turnover [dimensionless].

        Returns
        -------
            Dictionary with keys:

            - CflxLP:  Carbon flux from litter production [GtC/yr]
            - CflxLD: Carbon flux from litter decomposition [GtC/yr]
            - CflxSR: Carbon flux from soil respiration [GtC/yr]

        Notes
        -----
        The following formulas are used:

        .. math::

            cLP = \\frac{C_P}{\\tau_{C_P}} \\times \\epsilon_{T(cLP)} \\times
            \\epsilon_{N(cLP)} \\times \\epsilon_{LU(cLP)}

            cLD = \\frac{C_L}{\\tau_{C_L}} \\times \\epsilon_{T(cLD)} \\times
            \\epsilon_{N(cLD)} \\times \\epsilon_{LU(cLD)}

            cSR = \\frac{C_S}{\\tau_{C_S}} \\times \\epsilon_{T(cSR)} \\times
            \\epsilon_{N(cSR)} \\times \\epsilon_{LU(cSR)}

        where:

        - :math:`C_P`, :math:`C_L`, :math:`C_S` are the plant, litter, and soil
          carbon pool sizes
        - :math:`\\tau_{C_P}`, :math:`\\tau_{C_L}`, :math:`\\tau_{C_S}` are the
          respective turnover times
        - :math:`\\epsilon_{T}` represents temperature effects
        - :math:`\\epsilon_{N}` represents carbon-nitrogen coupling effects
        - :math:`\\epsilon_{LU}` represents land use effects

        Each flux is calculated using first-order kinetics where the turnover rate
        equals pool size divided by turnover time, modified by environmental effects.
        This formulation assumes that turnover processes follow exponential decay
        dynamics at baseline conditions, with multiplicative effects from environmental
        factors.
        """
        return {
            "CflxLP": CplsP / self.tau_CplsP * eff_dT_cLP * eff_N_cLP * eff_LU_cLP,
            "CflxLD": CplsL / self.tau_CplsL * eff_dT_cLD * eff_N_cLD * eff_LU_cLD,
            "CflxSR": CplsS / self.tau_CplsS * eff_dT_cSR * eff_N_cSR * eff_LU_cSR,
        }


@define
class CarbonCycleCalculator:
    """Calculator for terrestrial carbon cycle dynamics.

    This class handles the complete carbon cycle by tracking carbon transfers
    between pools and calculating pool dynamics. The carbon cycle includes three
    main terrestrial pools:

    - Plant carbon pool (CplsP): Carbon in living plant biomass
    - Litter carbon pool (CplsL): Carbon in dead plant material
    - Soil carbon pool (CplsS): Carbon in soil organic matter

    Key processes include:

    - Carbon fixation:  Net primary production (NPP) converting atmospheric CO₂
      to organic carbon
    - Autotrophic respiration: Litter production respiration (LPR) releasing CO₂
      from plant metabolism
    - Organic matter turnover: Litter production (cLP), litter decomposition (cLD),
      and soil respiration (cSR)
    - Heterotrophic respiration: Release of CO₂ from decomposition (RH = LPR + cLD2A + cSR)
    - Land use changes: Gross deforestation (cLUgrs) releasing carbon to atmosphere

    Each input flux is partitioned among pools using empirically-derived fractions
    that represent the distribution of carbon through different pathways in the
    terrestrial biosphere.  These fractions account for the cascade of carbon through
    multiple pools within the annual timestep, capturing both immediate allocation
    and subsequent rapid turnover processes.
    """

    frc_NPP2P: pint.Quantity = field(
        default=Q(0.55, "1"),
        validator=check_units("1"),
    )
    """Fraction of net primary production allocated to plant pool [dimensionless]. 
    
    Represents the portion of newly fixed carbon remaining in plant biomass at the
    end of the annual timestep, primarily in long-lived tissues (wood, structural
    tissues, storage organs) that do not turn over within the year.
    """

    frc_NPP2L: pint.Quantity = field(
        default=Q(0.35, "1"),
        validator=check_units("1"),
    )
    """Fraction of net primary production allocated to litter pool [dimensionless].
    
    Represents newly fixed carbon that enters the litter pool within the annual
    timestep through rapid plant-to-litter turnover (e.g., fine root turnover, leaf
    senescence). The remainder (1 - frc_NPP2P - frc_NPP2L) goes to the soil pool,
    representing carbon that cycles through both plant and litter pools and enters
    soil within the annual timestep.
    """

    frc_cLP2L: pint.Quantity = field(
        default=Q(0.7, "1"),
        validator=check_units("1"),
    )
    """Fraction of litter production carbon retained in litter pool [dimensionless]. 
    
    Represents carbon in plant materials that remain in the litter pool at the end
    of the annual timestep as recognizable plant structures.  The remainder (1 - frc_cLP2L)
    goes directly to soil, representing materials that decompose rapidly and enter
    the soil pool within the annual timestep, bypassing extended residence in the
    litter stage.
    """

    frc_cLD2S: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """Fraction of litter decomposition carbon transferred to soil pool [dimensionless]. 
    
    Represents carbon incorporated into stable soil organic matter during decomposition.
    The remainder (1 - frc_cLD2S) is respired to the atmosphere (cLD2A), representing
    carbon released as CO₂ during litter breakdown.
    """

    frc_cLUgrs2P: pint.Quantity = field(
        default=Q(0.6, "1"),
        validator=check_units("1"),
    )
    """Fraction of gross deforestation carbon loss from plant pool [dimensionless].
    
    Represents the distribution of carbon lost during land conversion, with this
    fraction coming from plant biomass removal or burning. 
    """

    frc_cLUgrs2L: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """Fraction of gross deforestation carbon loss from litter pool [dimensionless]. 
    
    Represents carbon lost from litter during land conversion. The remainder
    (1 - frc_cLUgrs2P - frc_cLUgrs2L) comes from the soil pool. 
    """

    def calculate_dCpls_dt(
            self,
            CflxNPP: pint.Quantity,
            CflxLPR: pint.Quantity,
            CflxLP: pint.Quantity,
            CflxLD: pint.Quantity,
            CflxSR: pint.Quantity,
            CflxLUgrs: pint.Quantity,
    ) -> List[pint.Quantity]:
        """
        Calculate rate of change for plant, litter, and soil carbon pools.

        This method computes the time derivatives of terrestrial carbon pools based
        on inputs (NPP), autotrophic respiration (LPR), internal transfers (cLP, cLD),
        heterotrophic respiration (cSR), and land use changes (cLUgrs). Each input
        flux is partitioned among pools according to the fractional allocation parameters,
        which account for carbon cascading through multiple pools within the annual timestep.

        Parameters
        ----------
        CflxNPP
            Net primary production carbon flux [GtC/yr].  Total carbon fixation by
            plants after autotrophic respiration.
        CflxLPR
            Litter production respiration carbon flux [GtC/yr].  Autotrophic respiration
            associated with plant tissue turnover.
        CflxLP
            Litter production flux [GtC/yr]. Carbon transfer from plant to litter
            through senescence and mortality.
        CflxLD
            Litter decomposition flux [GtC/yr]. Carbon release from litter through
            decomposition.
        CflxSR
            Soil respiration flux [GtC/yr]. Carbon release from soil organic matter
            through decomposition (also called soil turnover).
        CflxLUgrs
            Gross deforestation carbon flux [GtC/yr].  Carbon lost from system
            due to land conversion.

        Returns
        -------
            List containing:

            - dCplsP_dt: Rate of change of plant carbon pool [GtC/yr]
            - dCplsL_dt: Rate of change of litter carbon pool [GtC/yr]
            - dCplsS_dt: Rate of change of soil carbon pool [GtC/yr]

        Notes
        -----
        The pool dynamics follow mass balance equations:

        **Plant Pool:**

        .. math::

            \\frac{dC_P}{dt} = NPP \\times f_{NPP2P} - LPR - cLP - cLUgrs \\times f_{cLUgrs2P}

        **Litter Pool:**

        .. math::

            \\frac{dC_L}{dt} = NPP \\times f_{NPP2L} + cLP \\times f_{cLP2L} - cLD -
            cLUgrs \\times f_{cLUgrs2L}

        **Soil Pool:**

        .. math::

            \\frac{dC_S}{dt} = NPP \\times f_{NPP2S} + cLP \\times f_{cLP2S} +
            cLD \\times f_{cLD2S} - cSR - cLUgrs \\times f_{cLUgrs2S}

        where:

        - :math:`f_{NPP2S} = 1 - f_{NPP2P} - f_{NPP2L}` represents NPP carbon
          cascading through plant→litter→soil within the annual timestep
        - :math:`f_{cLP2S} = 1 - f_{cLP2L}` represents rapid litter decomposition
          directly entering soil within the annual timestep
        - :math:`f_{cLUgrs2S} = 1 - f_{cLUgrs2P} - f_{cLUgrs2L}` (remainder from soil)

        **Important Note on Annual Timestep:**

        All allocation fractions account for the cascade of carbon through multiple
        pools within the annual timestep. For example, :math:`f_{NPP2S}` does not
        represent direct allocation from atmosphere to soil, but rather carbon that
        is fixed by plants, then enters litter through turnover, and subsequently
        enters soil through decomposition—all within a single year.  This approach
        effectively captures rapid cycling processes without explicitly modeling
        sub-annual dynamics.
        """
        # Calculate NPP fractions to each pool
        CflxNPP2P = CflxNPP * self.frc_NPP2P
        CflxNPP2L = CflxNPP * self.frc_NPP2L
        CflxNPP2S = CflxNPP - CflxNPP2P - CflxNPP2L

        # Calculate land use emissions fractions to each pool
        CflxLUgrs2P = CflxLUgrs * self.frc_cLUgrs2P
        CflxLUgrs2L = CflxLUgrs * self.frc_cLUgrs2L
        CflxLUgrs2S = CflxLUgrs - CflxLUgrs2P - CflxLUgrs2L

        # Calculate litter production fractions
        CflxLP2L = CflxLP * self.frc_cLP2L
        CflxLP2S = CflxLP - CflxLP2L

        # Calculate litter decomposition fractions
        CflxLD2S = CflxLD * self.frc_cLD2S

        # Calculate pool size changes
        dCplsP_dt = CflxNPP2P - CflxLPR - CflxLP - CflxLUgrs2P
        dCplsL_dt = CflxNPP2L + CflxLP2L - CflxLD - CflxLUgrs2L
        dCplsS_dt = CflxNPP2S + CflxLP2S + CflxLD2S - CflxSR - CflxLUgrs2S

        return [dCplsP_dt, dCplsL_dt, dCplsS_dt]

    def calculate_Cflx_all(
            self,
            CflxNPP: pint.Quantity,
            CflxLPR: pint.Quantity,
            CflxLP: pint.Quantity,
            CflxLD: pint.Quantity,
            CflxSR: pint.Quantity,
            CflxLUgrs: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate all carbon fluxes and their partitioning in the system.

        This method computes the complete set of carbon fluxes including all
        partitioned flows between pools.  Useful for detailed budget analysis and
        diagnostic output.

        Parameters
        ----------
        CflxNPP
            Net primary production carbon flux [GtC/yr].
        CflxLPR
            Litter production respiration carbon flux [GtC/yr].
        CflxLP
            Litter production flux [GtC/yr].
        CflxLD
            Litter decomposition flux [GtC/yr].
        CflxSR
            Soil respiration flux [GtC/yr] (also called soil turnover).
        CflxLUgrs
            Gross deforestation carbon flux [GtC/yr].

        Returns
        -------
            Dictionary containing all carbon fluxes and their partitioned components:

            - Primary fluxes: CflxNPP, CflxLPR, CflxLP, CflxLD, CflxSR, CflxLUgrs
            - NPP partitioning: CflxNPP2P, CflxNPP2L, CflxNPP2S
            - LP partitioning: CflxLP2L, CflxLP2S
            - LD partitioning: CflxLD2S, CflxLD2A
            - LUgrs partitioning: CflxLUgrs2P, CflxLUgrs2L, CflxLUgrs2S
            - Derived flux: CflxRH (total heterotrophic respiration = LPR + cLD2A + cSR)

        Notes
        -----
        Heterotrophic respiration (RH) is calculated as:

        .. math::

            RH = LPR + cLD2A + cSR

        where:

        - :math:`LPR` is litter production respiration (plant autotrophic respiration)
        - :math:`cLD2A = cLD \\times (1 - f_{cLD2S})` is litter decomposition to atmosphere
        - :math:`cSR` is soil respiration (all assumed to be respired to atmosphere)

        This represents the total flux of carbon from terrestrial organic matter to
        the atmosphere through respiration processes, which is a key component of
        the terrestrial carbon budget and climate feedbacks.
        """
        # Calculate component fluxes
        # Calculate NPP fractions to each pool
        CflxNPP2P = CflxNPP * self.frc_NPP2P
        CflxNPP2L = CflxNPP * self.frc_NPP2L
        CflxNPP2S = CflxNPP - CflxNPP2P - CflxNPP2L

        # Calculate land use emissions fractions
        CflxLUgrs2P = CflxLUgrs * self.frc_cLUgrs2P
        CflxLUgrs2L = CflxLUgrs * self.frc_cLUgrs2L
        CflxLUgrs2S = CflxLUgrs - CflxLUgrs2P - CflxLUgrs2L

        CflxLP2L = CflxLP * self.frc_cLP2L
        CflxLP2S = CflxLP - CflxLP2L

        CflxLD2S = CflxLD * self.frc_cLD2S
        CflxLD2A = CflxLD - CflxLD2S

        # Calculate heterotrophic respiration
        CflxRH = CflxLPR + CflxLD2A + CflxSR

        # Merge all dictionaries and return
        return {
            "CflxNPP": CflxNPP,
            "CflxNPP2P": CflxNPP2P,
            "CflxNPP2L": CflxNPP2L,
            "CflxNPP2S": CflxNPP2S,
            "CflxLPR": CflxLPR,
            "CflxLP": CflxLP,
            "CflxLP2L": CflxLP2L,
            "CflxLP2S": CflxLP2S,
            "CflxLD": CflxLD,
            "CflxLD2S": CflxLD2S,
            "CflxLD2A": CflxLD2A,
            "CflxSR": CflxSR,
            "CflxRH": CflxRH,
            "CflxLUgrs": CflxLUgrs,
            "CflxLUgrs2P": CflxLUgrs2P,
            "CflxLUgrs2L": CflxLUgrs2L,
            "CflxLUgrs2S": CflxLUgrs2S,
        }


@define
class CarbonCycleModelResult:
    """Results from running the carbon cycle model.

    This class stores the complete output from a carbon cycle model simulation,
    including carbon pool states, fluxes, and all environmental effect modifiers.
    Results can be converted to an xarray Dataset for analysis and visualization.

    The carbon cycle model tracks three terrestrial carbon pools:

    - Plant carbon pool (CplsP): Carbon in living plant biomass
    - Litter carbon pool (CplsL): Carbon in dead plant material
    - Soil carbon pool (CplsS): Carbon in soil organic matter

    Key processes represented include:

    - Carbon fixation through net primary production (NPP)
    - Litter production respiration (LPR)
    - Organic matter turnover (cLP, cLD, cSR)
    - Land use change emissions (cLUgrs)
    - Environmental modifiers: temperature, nitrogen, and land use effects

    Note:
        All `_t` attributes are time-continuous functions created through linear
        interpolation of discrete time points. These functions map from time [yr]
        to their respective units, allowing for smooth temporal interpolation of
        forcing data and effect modifiers.
    """

    c_state: Dict[str, pint.Quantity]
    """Carbon pool states over time [GtC].  
    
    Dictionary containing time series of carbon pools: 
    - 'CplsP': Plant carbon pool
    - 'CplsL': Litter carbon pool
    - 'CplsS':  Soil carbon pool
    """

    calc_c_turnover: CarbonTurnoverCalculator
    """Calculator for carbon turnover processes. 
    
    Instance of CarbonTurnoverCalculator containing turnover time parameters
    and methods for calculating organic matter breakdown fluxes (cLP, cLD, cSR).
    """

    calc_c_cycle: CarbonCycleCalculator
    """Calculator for carbon cycle mass balance.
    
    Instance of CarbonCycleCalculator containing allocation fractions and
    methods for calculating flux partitioning and pool dynamics.
    """

    time_axis: pint.Quantity
    """Time points for simulation [yr].
    
    Array of discrete time points at which the model state is evaluated,
    typically representing annual timesteps.
    """

    # Time-interpolated flux functions
    CflxNPP_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated net primary production flux function [GtC/yr]. 
    
    Continuous function mapping time to NPP flux. Interpolated
    from discrete forcing data to allow smooth temporal variation.
    """

    CflxLPR_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated litter production respiration flux function [GtC/yr].
    
    Continuous function mapping time to LPR flux. Interpolated from discrete
    forcing data. 
    """

    CflxLUgrs_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated gross land use change flux function [GtC/yr].
    
    Continuous function mapping time to land use emissions. Represents carbon
    lost to atmosphere during land conversion (deforestation, agricultural expansion).
    Interpolated from discrete forcing data.
    """

    # Time-interpolated effect functions
    eff_dT_cLP_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated temperature effect on litter production function [dimensionless].
    
    Continuous function mapping time to temperature modifier for plant-to-litter
    carbon turnover. Interpolated from calculated temperature effects.
    """

    eff_dT_cLD_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated temperature effect on litter decomposition function [dimensionless].
    
    Continuous function mapping time to temperature modifier for litter-to-soil
    carbon turnover. Interpolated from calculated temperature effects.
    """

    eff_dT_cSR_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated temperature effect on soil respiration function [dimensionless].
    
    Continuous function mapping time to temperature modifier for soil organic
    matter decomposition. Interpolated from calculated temperature effects.
    """

    eff_N_cLP_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated nitrogen effect on litter production function [dimensionless].
    
    Continuous function mapping time to nitrogen modifier for plant-to-litter carbon 
    turnover. Represents carbon-nitrogen coupling effects on plant tissue turnover.
    Interpolated from calculated nitrogen effects.
    """

    eff_N_cLD_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated nitrogen effect on litter decomposition function [dimensionless].
    
    Continuous function mapping time to nitrogen modifier for litter decomposition. 
    Represents how nitrogen limitation affects decomposer activity. Interpolated from 
    calculated nitrogen effects.
    """

    eff_N_cSR_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated nitrogen effect on soil respiration function [dimensionless]. 
    
    Continuous function mapping time to nitrogen modifier for soil organic matter 
    decomposition. Represents carbon-nitrogen coupling in soil respiration. 
    Interpolated from calculated nitrogen effects.
    """

    eff_LU_cLP_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated land use effect on litter production function [dimensionless].
    
    Continuous function mapping time to land use modifier for plant-to-litter
    carbon turnover. Represents how land use change affects plant mortality rates.
    Interpolated from calculated land use effects.
    """

    eff_LU_cLD_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated land use effect on litter decomposition function [dimensionless]. 
    
    Continuous function mapping time to land use modifier for litter decomposition.
    Represents how land use change affects litter breakdown rates. Interpolated
    from calculated land use effects. 
    """

    eff_LU_cSR_t: Callable[[pint.Quantity], pint.Quantity]
    """Time-interpolated land use effect on soil respiration function [dimensionless].
    
    Continuous function mapping time to land use modifier for soil organic matter
    decomposition. Represents how land use change affects soil carbon dynamics.
    Interpolated from calculated land use effects.
    """

    def add_non_state_variables(self) -> xr.Dataset:
        """
        Calculate derived variables and create xarray Dataset with complete results.

        This method computes all diagnostic fluxes and derived quantities from the
        model state variables and returns a complete xarray Dataset containing:

        - State variables:  Carbon pools (CplsP, CplsL, CplsS)
        - Turnover fluxes:  Calculated from pools and environmental effects (cLP, cLD, cSR)
        - All partitioned fluxes: NPP, LPR, and land use change distributed among pools
        - Derived quantities: Total terrestrial carbon (CplsPLS), net terrestrial flux (CflxNetPLS)
        - Heterotrophic respiration:  Total CO₂ release from decomposition (CflxRH)

        Returns
        -------
        xr.Dataset
            Complete model results as xarray Dataset with time dimension and all
            variables labeled with descriptive names and units.  The Dataset includes:

            **State Variables:**

            - CplsP, CplsL, CplsS: Individual pool sizes [GtC]
            - CplsPLS: Total terrestrial carbon [GtC]

            **Primary Fluxes:**

            - CflxNPP: Net primary production [GtC/yr]
            - CflxLPR:  Litter production respiration [GtC/yr]
            - CflxLP: Litter production [GtC/yr]
            - CflxLD: Litter decomposition [GtC/yr]
            - CflxSR: Soil respiration [GtC/yr]
            - CflxLUgrs:  Gross land use emissions [GtC/yr]

            **Partitioned Fluxes:**

            - CflxNPP2P, CflxNPP2L, CflxNPP2S: NPP allocation to each pool
            - CflxLP2L, CflxLP2S:  Litter production allocation
            - CflxLD2S, CflxLD2A: Litter decomposition to soil and atmosphere
            - CflxLUgrs2P, CflxLUgrs2L, CflxLUgrs2S: Land use emissions from each pool

            **Derived Fluxes:**

            - CflxRH: Total heterotrophic respiration [GtC/yr]
            - CflxNetPLS: Net change in terrestrial carbon [GtC/yr]

        Notes
        -----
        The net terrestrial carbon flux is calculated using numerical differentiation:

        .. math::

            cNetPLS = \\frac{dC_{PLS}}{dt} = \\frac{d(C_P + C_L + C_S)}{dt}

        This represents the net carbon balance of the terrestrial biosphere, with
        positive values indicating carbon uptake (sink) and negative values indicating
        carbon release (source).

        The heterotrophic respiration is calculated as:

        .. math::

            RH = LPR + cLD2A + cSR

        This represents total CO₂ release from decomposition processes, a key
        component of the terrestrial carbon budget and climate feedbacks.
        """
        state = self.c_state
        time = self.time_axis

        # Calculate turnover
        turnover = self.calc_c_turnover.calculate(
            CplsP=state["CplsP"],
            CplsL=state["CplsL"],
            CplsS=state["CplsS"],
            eff_dT_cLP=self.eff_dT_cLP_t(time),
            eff_dT_cLD=self.eff_dT_cLD_t(time),
            eff_dT_cSR=self.eff_dT_cSR_t(time),
            eff_N_cLP=self.eff_N_cLP_t(time),
            eff_N_cLD=self.eff_N_cLD_t(time),
            eff_N_cSR=self.eff_N_cSR_t(time),
            eff_LU_cLP=self.eff_LU_cLP_t(time),
            eff_LU_cLD=self.eff_LU_cLD_t(time),
            eff_LU_cSR=self.eff_LU_cSR_t(time),
        )

        # Calculate all fluxes
        flux = self.calc_c_cycle.calculate_Cflx_all(
            CflxNPP=self.CflxNPP_t(time),
            CflxLPR=self.CflxLPR_t(time),
            CflxLP=turnover["CflxLP"],
            CflxLD=turnover["CflxLD"],
            CflxSR=turnover["CflxSR"],
            CflxLUgrs=self.CflxLUgrs_t(time),
        )

        extra_var = {"CplsPLS": state["CplsP"] + state["CplsL"] + state["CplsS"]}
        extra_var["CflxNetPLS"] = Q(np.gradient(extra_var["CplsPLS"].m), "GtC/yr")

        return make_dataset_from_var_dict({**state, **flux, **extra_var}, time)


@define
class CarbonCycleModel:
    """MAGICC's terrestrial carbon cycle model (carbon component of CNit).

    This model simulates the dynamics of three terrestrial carbon pools in the
    terrestrial biosphere by solving a system of ordinary differential equations (ODEs)
    using scipy.integrate.solve_ivp.  The model tracks carbon transfers between pools
    through processes including:

    - Carbon fixation via net primary production (NPP)
    - Litter production respiration (LPR)
    - Organic matter turnover (litter production, litter decomposition, soil respiration)
    - Land use change emissions
    - Environmental modifiers (temperature, nitrogen availability, land use effects)

    The model uses annual timesteps and accounts for rapid carbon cycling through
    empirical allocation fractions that represent carbon cascading through multiple
    pools within a year.

    Time-continuous functions (suffixed with `_t`) are generated through linear
    interpolation of discrete time points and provide continuous values over the
    simulation period, allowing the ODE solver to evaluate fluxes at any time point.
    """

    calc_c_turnover: CarbonTurnoverCalculator
    """Calculator for carbon turnover processes. 
    
    Instance of CarbonTurnoverCalculator containing turnover time parameters
    (tau_CplsP, tau_CplsL, tau_CplsS) and methods for calculating organic matter
    breakdown fluxes based on pool sizes and environmental effects.
    """

    calc_c_cycle: CarbonCycleCalculator
    """Calculator for carbon cycle mass balance. 
    
    Instance of CarbonCycleCalculator containing allocation fractions
    (frc_NPP2P, frc_NPP2L, frc_cLP2L, frc_cLD2S, frc_cLUgrs2P, frc_cLUgrs2L)
    and methods for calculating flux partitioning and pool dynamics.
    """

    CplsP0: pint.Quantity = field(
        default=Q(400, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial plant carbon pool size [GtC]. 
    
    Carbon contained in living plant biomass at the start of the simulation
    (time0). Typical values range from 400-550 GtC for global simulations.
    """

    CplsL0: pint.Quantity = field(
        default=Q(50, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial litter carbon pool size [GtC]. 
    
    Carbon contained in dead plant material (litter) at the start of the simulation
    (time0). Typical values range from 40-70 GtC for global simulations.
    """

    CplsS0: pint.Quantity = field(
        default=Q(1550, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial soil carbon pool size [GtC].
    
    Carbon contained in soil organic matter at the start of the simulation (time0).
    Typical values range from 1400-1700 GtC for global simulations.
    """

    time0: pint.Quantity = field(
        default=Q(1850, "yr"),
        validator=check_units("yr"),
    )
    """Initialization time [yr].
    
    Time point at which the initial pool sizes (CplsP0, CplsL0, CplsS0) apply.
    The simulation time_axis must start at or after this time. Typically set to
    the pre-industrial baseline year (1850).
    """

    switch_Cpls: list[int] = field(
        default=[1, 1, 1],
    )
    """Carbon pool switches for [plant, litter, soil] [dimensionless].
    
    Binary switches (1=enabled, 0=disabled) to selectively enable or disable
    pool dynamics. Useful for diagnostic purposes or simplified model configurations.
    Setting a switch to 0 freezes that pool at its initial value.
    
    - switch_Cpls[0]:  Plant pool switch
    - switch_Cpls[1]: Litter pool switch
    - switch_Cpls[2]: Soil pool switch
    """

    def run(
            self,
            time_axis: pint.Quantity,
            CflxNPP_t: Callable[[pint.Quantity], pint.Quantity],
            CflxLPR_t: Callable[[pint.Quantity], pint.Quantity],
            CflxLUgrs_t: Callable[[pint.Quantity], pint.Quantity],
            eff_dT_cLP_t: Callable[[pint.Quantity], pint.Quantity],
            eff_dT_cLD_t: Callable[[pint.Quantity], pint.Quantity],
            eff_dT_cSR_t: Callable[[pint.Quantity], pint.Quantity],
            eff_N_cLP_t: Callable[[pint.Quantity], pint.Quantity],
            eff_N_cLD_t: Callable[[pint.Quantity], pint.Quantity],
            eff_N_cSR_t: Callable[[pint.Quantity], pint.Quantity],
            eff_LU_cLP_t: Callable[[pint.Quantity], pint.Quantity],
            eff_LU_cLD_t: Callable[[pint.Quantity], pint.Quantity],
            eff_LU_cSR_t: Callable[[pint.Quantity], pint.Quantity],
    ) -> CarbonCycleModelResult:
        """
        Run the carbon cycle model simulation.

        This method solves the system of ODEs for carbon pool dynamics over the
        specified time period using scipy.integrate.solve_ivp.  The solver uses
        time-continuous interpolated functions for all forcing data and environmental
        effects, allowing flexible time step selection during integration.

        Parameters
        ----------
        time_axis
            Time points for simulation [yr].  Must start at or after time0.
            These are the time points at which solution is explicitly evaluated.
        CflxNPP_t
            Time-interpolated net primary production flux function [GtC/yr].
            Maps time to NPP, representing carbon fixation by plants.
        CflxLPR_t
            Time-interpolated litter production respiration flux function [GtC/yr].
            Maps time to LPR, representing autotrophic respiration from tissue turnover.
        CflxLUgrs_t
            Time-interpolated gross land use change flux function [GtC/yr].
            Maps time to land use emissions from deforestation and land conversion.
        eff_dT_cLP_t
            Time-interpolated temperature effect on litter production function [dimensionless].
            Maps time to temperature modifier for plant-to-litter turnover.
        eff_dT_cLD_t
            Time-interpolated temperature effect on litter decomposition function [dimensionless].
            Maps time to temperature modifier for litter-to-soil turnover.
        eff_dT_cSR_t
            Time-interpolated temperature effect on soil respiration function [dimensionless].
            Maps time to temperature modifier for soil organic matter decomposition.
        eff_N_cLP_t
            Time-interpolated nitrogen effect on litter production function [dimensionless].
            Maps time to nitrogen modifier for plant-to-litter turnover.
        eff_N_cLD_t
            Time-interpolated nitrogen effect on litter decomposition function [dimensionless].
            Maps time to nitrogen modifier for litter decomposition.
        eff_N_cSR_t
            Time-interpolated nitrogen effect on soil respiration function [dimensionless].
            Maps time to nitrogen modifier for soil respiration.
        eff_LU_cLP_t
            Time-interpolated land use effect on litter production function [dimensionless].
            Maps time to land use modifier for plant-to-litter turnover.
        eff_LU_cLD_t
            Time-interpolated land use effect on litter decomposition function [dimensionless].
            Maps time to land use modifier for litter decomposition.
        eff_LU_cSR_t
            Time-interpolated land use effect on soil respiration function [dimensionless].
            Maps time to land use modifier for soil respiration.

        Returns
        -------
        CarbonCycleModelResult
            Simulation results containing:

            - Carbon pool time series (CplsP, CplsL, CplsS)
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
        **ODE System:**

        The model solves the following system of ODEs:

        .. math::

            \\frac{dC_P}{dt} = NPP \\times f_{NPP2P} - LPR - cLP - cLUgrs \\times f_{cLUgrs2P}

            \\frac{dC_L}{dt} = NPP \\times f_{NPP2L} + cLP \\times f_{cLP2L} - cLD -
            cLUgrs \\times f_{cLUgrs2L}

            \\frac{dC_S}{dt} = NPP \\times f_{NPP2S} + cLP \\times f_{cLP2S} +
            cLD \\times f_{cLD2S} - cSR - cLUgrs \\times f_{cLUgrs2S}

        where turnover fluxes (cLP, cLD, cSR) are calculated from pool sizes and
        environmental effects using first-order kinetics.

        **Numerical Methods:**

        The solver uses scipy.integrate.solve_ivp with:

        - Adaptive time stepping (solver selects appropriate substeps)
        - Absolute tolerance:  1e-6 GtC
        - Relative tolerance: 1e-3
        - Default solver method (typically RK45)

        **Pool Switches:**

        Pool dynamics can be selectively disabled using switch_Cpls. When a switch
        is set to 0, the corresponding pool's rate of change is forced to zero,
        effectively freezing that pool at its initial value throughout the simulation.

        **Time-Continuous Functions:**

        All `_t` arguments are callable functions created through linear interpolation
        of discrete time points. This allows the ODE solver to evaluate forcing data
        and environmental effects at any time point during integration, not just at
        the discrete output times specified in time_axis.
        """

        def func_to_solve(t: float, y: np.ndarray) -> list:
            """
            ODE system function for scipy.integrate.solve_ivp.

            Parameters
            ----------
            t
                Current time [yr] (unitless for solver)
            y
                Current state vector [CplsP, CplsL, CplsS] in GtC (unitless for solver)

            Returns
            -------
            list
                Time derivatives [dCplsP/dt, dCplsL/dt, dCplsS/dt] in GtC/yr (unitless)
            """
            t = Q(t, "yr")

            # Current state with units
            CplsP = Q(y[0], "GtC")
            CplsL = Q(y[1], "GtC")
            CplsS = Q(y[2], "GtC")

            # Calculate turnover fluxes
            Cflx_turnover = self.calc_c_turnover.calculate(
                CplsP=CplsP,
                CplsL=CplsL,
                CplsS=CplsS,
                eff_dT_cLP=eff_dT_cLP_t(t),
                eff_dT_cLD=eff_dT_cLD_t(t),
                eff_dT_cSR=eff_dT_cSR_t(t),
                eff_N_cLP=eff_N_cLP_t(t),
                eff_N_cLD=eff_N_cLD_t(t),
                eff_N_cSR=eff_N_cSR_t(t),
                eff_LU_cLP=eff_LU_cLP_t(t),
                eff_LU_cLD=eff_LU_cLD_t(t),
                eff_LU_cSR=eff_LU_cSR_t(t),
            )

            # Calculate rate of change
            dydt = self.calc_c_cycle.calculate_dCpls_dt(
                CflxNPP=CflxNPP_t(t),
                CflxLPR=CflxLPR_t(t),
                CflxLUgrs=CflxLUgrs_t(t),
                CflxLP=Cflx_turnover["CflxLP"],
                CflxLD=Cflx_turnover["CflxLD"],
                CflxSR=Cflx_turnover["CflxSR"],
            )
            # Apply pool switches
            dydt = [v.to("GtC/yr").m * switch for v, switch in
                    zip(dydt, self.switch_Cpls)]
            if np.all(np.abs(dydt) < 1e-12):
                dydt = [0, 0, 0]
            return dydt

        if time_axis[0] < self.time0:
            raise ValueError(
                f"time_axis starts before time0: {time_axis[0]} < {self.time0}"
            )

        # Prepare solver arguments
        t_eval = time_axis.to("yr").m
        t_span = (self.time0.to("yr").m, t_eval[-1])

        # Solve the system
        raw = scipy.integrate.solve_ivp(
            func_to_solve,
            t_span=t_span,
            t_eval=t_eval,
            y0=(
                self.CplsP0.to("GtC").m,
                self.CplsL0.to("GtC").m,
                self.CplsS0.to("GtC").m,
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

        # Create result collection
        c_state = {
            "CplsP": Q(raw.y[0, :], "GtC"),
            "CplsL": Q(raw.y[1, :], "GtC"),
            "CplsS": Q(raw.y[2, :], "GtC"),
        }

        return CarbonCycleModelResult(
            c_state=c_state,
            calc_c_turnover=self.calc_c_turnover,
            calc_c_cycle=self.calc_c_cycle,
            time_axis=time_axis,
            CflxNPP_t=CflxNPP_t,
            CflxLPR_t=CflxLPR_t,
            CflxLUgrs_t=CflxLUgrs_t,
            eff_dT_cLP_t=eff_dT_cLP_t,
            eff_dT_cLD_t=eff_dT_cLD_t,
            eff_dT_cSR_t=eff_dT_cSR_t,
            eff_N_cLP_t=eff_N_cLP_t,
            eff_N_cLD_t=eff_N_cLD_t,
            eff_N_cSR_t=eff_N_cSR_t,
            eff_LU_cLP_t=eff_LU_cLP_t,
            eff_LU_cLD_t=eff_LU_cLD_t,
            eff_LU_cSR_t=eff_LU_cSR_t,
        )
