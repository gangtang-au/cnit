"""
MAGICC Terrestrial Carbon Cycle Model Implementation

This module provides the implementation of the carbon cycle component of MAGICC's
coupled carbon-nitrogen cycle model (CNit). It simulates the dynamics of terrestrial
carbon pools and their responses to climate change, CO₂ fertilization, nitrogen
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

2. **Allocation Fractions with Cascading**: Empirical allocation fractions represent
   carbon transfers through multiple pools within each model timestep. Annual timesteps
   are used in the default examples; evenly spaced subannual time axes are also supported.

3. **Environmental Modifiers**: All process rates are modified by dimensionless effect factors
   representing:

   - Temperature effects on decomposition and respiration
   - CO₂ fertilization effects on NPP
   - Nitrogen availability effects (carbon-nitrogen coupling)
   - Land use change effects on turnover rates

   (see :mod:`effect_factor_calculators`)

4. **First-Order Kinetics**: Carbon turnover from each pool follows first-order kinetics
   (flux = pool size / turnover time), modified by environmental effects.

5. **Flux Partitioning**: Input fluxes (NPP, land use) are partitioned among pools using
   empirical fractions that capture both direct allocation and rapid cycling over their
   model timestep.

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

from typing import Dict, List
from attrs import define, field
import numpy as np
import scipy
import xarray as xr
import pint

from ..utils.units import check_units, Q
from ..utils.maths import SolveError, _make_packed_numeric_interp
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
    """Base net primary production (NPP) rate under reference conditions, before
    environmental effects [GtC/yr].
    """

    LPR0: pint.Quantity = field(
        default=Q(0, "GtC/yr"),
        validator=check_units("GtC/yr"),
    )
    """Base litter production respiration (LPR) rate under reference conditions,
    before environmental effects [GtC/yr].
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
            CO₂ effect on NPP [dimensionless]
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

            NPP = (NPP_0 + LUrgr) \\times \\epsilon_{T(NPP)} \\times \\epsilon_{CO_2}
            \\times \\epsilon_{N(NPP)} \\times \\epsilon_{LU(NPP)}

            LPR = LPR0 \\times \\epsilon_{T(LPR)} \\times \\epsilon_{CO_2} \\times
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

            - CflxLP: Carbon flux from litter production [GtC/yr]
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

    def calculate_numeric(
        self,
        CplsP_GtC: float,
        CplsL_GtC: float,
        CplsS_GtC: float,
        tau_CplsP_yr: float,
        tau_CplsL_yr: float,
        tau_CplsS_yr: float,
        eff_dT_cLP: float,
        eff_dT_cLD: float,
        eff_dT_cSR: float,
        eff_N_cLP: float,
        eff_N_cLD: float,
        eff_N_cSR: float,
        eff_LU_cLP: float,
        eff_LU_cLD: float,
        eff_LU_cSR: float,
    ) -> Dict[str, float]:
        """Numeric counterpart of :meth:`calculate` for ODE integration.

        This applies the same turnover equations as :meth:`calculate`, using plain
        numeric values to avoid Pint overhead in the solver callback. Pool values
        represent GtC, turnover times represent years, and effects are dimensionless.
        Returned numeric fluxes represent GtC/yr.
        """
        return {
            "CflxLP": CplsP_GtC / tau_CplsP_yr * eff_dT_cLP * eff_N_cLP * eff_LU_cLP,
            "CflxLD": CplsL_GtC / tau_CplsL_yr * eff_dT_cLD * eff_N_cLD * eff_LU_cLD,
            "CflxSR": CplsS_GtC / tau_CplsS_yr * eff_dT_cSR * eff_N_cSR * eff_LU_cSR,
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

    - Carbon fixation: Net primary production (NPP) converting atmospheric CO₂
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
    multiple pools within the model timestep, capturing both immediate allocation
    and subsequent rapid turnover processes.
    """

    frc_NPP2P: pint.Quantity = field(
        default=Q(0.55, "1"),
        validator=check_units("1"),
    )
    """Fraction of net primary production allocated to plant pool [dimensionless].

    Represents the portion of newly fixed carbon remaining in plant biomass at the
    end of a model timestep, primarily in long-lived tissues (wood, structural
    tissues, storage organs) that do not turn over within the model timestep.
    """

    frc_NPP2L: pint.Quantity = field(
        default=Q(0.35, "1"),
        validator=check_units("1"),
    )
    """Fraction of net primary production allocated to litter pool [dimensionless].

    Represents newly fixed carbon that enters the litter pool within the model
    timestep through rapid plant-to-litter turnover (e.g., fine root turnover, leaf
    senescence). The remainder (1 - frc_NPP2P - frc_NPP2L) goes to the soil pool,
    representing carbon that cycles through both plant and litter pools and enters
    soil within the model timestep.
    """

    frc_cLP2L: pint.Quantity = field(
        default=Q(0.7, "1"),
        validator=check_units("1"),
    )
    """Fraction of litter production carbon retained in litter pool [dimensionless].

    Represents carbon in plant materials that remain in the litter pool at the end
    of a model timestep as recognizable plant structures.  The remainder (1 - frc_cLP2L)
    goes directly to soil, representing materials that decompose rapidly and enter
    the soil pool within the model timestep, bypassing extended residence in the
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
        which represent rapid carbon cascading through multiple pools within the model
        parameterization.

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
          cascading through plant→litter→soil within the model timestep
        - :math:`f_{cLP2S} = 1 - f_{cLP2L}` represents rapid litter decomposition
          directly entering soil within the model timestep
        - :math:`f_{cLUgrs2S} = 1 - f_{cLUgrs2P} - f_{cLUgrs2L}` (remainder from soil)

        **Allocation Fractions and Time Resolution:**

        All allocation fractions account for the rapid cascade of carbon through
        multiple pools over the model timestep. For example,
        :math:`f_{NPP2S}` does not represent direct allocation from atmosphere to
        soil, but rather carbon that is fixed by plants, then enters litter through
        turnover, and subsequently enters soil through decomposition over that model
        timestep.
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

    def calculate_dCpls_dt_numeric(
        self,
        CflxNPP_GtC_yr: float,
        CflxLPR_GtC_yr: float,
        CflxLP_GtC_yr: float,
        CflxLD_GtC_yr: float,
        CflxSR_GtC_yr: float,
        CflxLUgrs_GtC_yr: float,
        frc_NPP2P: float,
        frc_NPP2L: float,
        frc_cLP2L: float,
        frc_cLD2S: float,
        frc_cLUgrs2P: float,
        frc_cLUgrs2L: float,
    ) -> np.ndarray:
        """Numeric counterpart of :meth:`calculate_dCpls_dt` for ODE integration.

        This applies the same pool-balance equations as :meth:`calculate_dCpls_dt`,
        using plain numeric values to avoid Pint overhead in the solver callback.
        Fluxes represent GtC/yr and allocation fractions are dimensionless. The
        returned array contains numeric derivatives representing GtC/yr.
        """
        CflxNPP2P = CflxNPP_GtC_yr * frc_NPP2P
        CflxNPP2L = CflxNPP_GtC_yr * frc_NPP2L
        CflxNPP2S = CflxNPP_GtC_yr - CflxNPP2P - CflxNPP2L

        CflxLUgrs2P = CflxLUgrs_GtC_yr * frc_cLUgrs2P
        CflxLUgrs2L = CflxLUgrs_GtC_yr * frc_cLUgrs2L
        CflxLUgrs2S = CflxLUgrs_GtC_yr - CflxLUgrs2P - CflxLUgrs2L

        CflxLP2L = CflxLP_GtC_yr * frc_cLP2L
        CflxLP2S = CflxLP_GtC_yr - CflxLP2L
        CflxLD2S = CflxLD_GtC_yr * frc_cLD2S

        return np.array(
            [
                CflxNPP2P - CflxLPR_GtC_yr - CflxLP_GtC_yr - CflxLUgrs2P,
                CflxNPP2L + CflxLP2L - CflxLD_GtC_yr - CflxLUgrs2L,
                CflxNPP2S + CflxLP2S + CflxLD2S - CflxSR_GtC_yr - CflxLUgrs2S,
            ],
            dtype=float,
        )

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
    """Results returned by :meth:`CarbonCycleModel.run`.

    Stores the time-aligned carbon pool states and the flux and effect series used
    during the simulation. Use :meth:`add_non_state_variables` to calculate derived
    fluxes and diagnostics and assemble the complete xarray dataset. Attribute
    descriptions below give the contents and units of the stored values.
    """

    c_state: Dict[str, pint.Quantity]
    """Carbon pool states over time [GtC].

    Dictionary containing time series of carbon pools:

    - 'CplsP': Plant carbon pool
    - 'CplsL': Litter carbon pool
    - 'CplsS': Soil carbon pool
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

    Array of evenly spaced time points at which the model state is evaluated.
    Annual timesteps are used in the default examples; subannual time axes are also
    supported.
    """

    CflxNPP_s: pint.Quantity
    """Net primary production flux series [GtC/yr].

    Series values are aligned with ``time_axis``.
    """

    CflxLPR_s: pint.Quantity
    """Litter production respiration flux series [GtC/yr].

    Series values are aligned with ``time_axis``.
    """

    CflxLUgrs_s: pint.Quantity
    """Gross land use carbon emissions series [GtC/yr].

    Series values are aligned with ``time_axis``.
    """

    eff_dT_cLP_s: pint.Quantity
    """Temperature effect on litter production carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_dT_cLD_s: pint.Quantity
    """Temperature effect on litter decomposition carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_dT_cSR_s: pint.Quantity
    """Temperature effect on soil respiration carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_N_cLP_s: pint.Quantity
    """Nitrogen effect on litter production carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_N_cLD_s: pint.Quantity
    """Nitrogen effect on litter decomposition carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_N_cSR_s: pint.Quantity
    """Nitrogen effect on soil respiration carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_LU_cLP_s: pint.Quantity
    """Land use effect on litter production carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_LU_cLD_s: pint.Quantity
    """Land use effect on litter decomposition carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    eff_LU_cSR_s: pint.Quantity
    """Land use effect on soil respiration carbon turnover [dimensionless].

    Series values are aligned with ``time_axis``.
    """

    def add_non_state_variables(self) -> xr.Dataset:
        """
        Calculate derived variables and create xarray Dataset with complete results.

        This method computes all diagnostic fluxes and derived quantities from the
        model state variables and returns a complete xarray Dataset containing:

        - State variables: Carbon pools (CplsP, CplsL, CplsS)
        - Turnover fluxes: Calculated from pools and environmental effects (cLP, cLD, cSR)
        - All partitioned fluxes: NPP, LPR, and land use change distributed among pools
        - Derived quantities: Total terrestrial carbon (CplsPLS), net terrestrial flux (CflxNetPLS)
        - Heterotrophic respiration: Total CO₂ release from decomposition (CflxRH)

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
            - CflxLPR: Litter production respiration [GtC/yr]
            - CflxLP: Litter production [GtC/yr]
            - CflxLD: Litter decomposition [GtC/yr]
            - CflxSR: Soil respiration [GtC/yr]
            - CflxLUgrs: Gross land use emissions [GtC/yr]

            **Partitioned Fluxes:**

            - CflxNPP2P, CflxNPP2L, CflxNPP2S: NPP allocation to each pool
            - CflxLP2L, CflxLP2S: Litter production allocation
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
        positive values indicating carbon uptake (a sink) and negative values
        indicating carbon release (a source). **In terrestrial carbon cycle studies,
        this is equivalent to net biome production (NBP)**.

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
            eff_dT_cLP=self.eff_dT_cLP_s,
            eff_dT_cLD=self.eff_dT_cLD_s,
            eff_dT_cSR=self.eff_dT_cSR_s,
            eff_N_cLP=self.eff_N_cLP_s,
            eff_N_cLD=self.eff_N_cLD_s,
            eff_N_cSR=self.eff_N_cSR_s,
            eff_LU_cLP=self.eff_LU_cLP_s,
            eff_LU_cLD=self.eff_LU_cLD_s,
            eff_LU_cSR=self.eff_LU_cSR_s,
        )

        # Calculate all fluxes
        flux = self.calc_c_cycle.calculate_Cflx_all(
            CflxNPP=self.CflxNPP_s,
            CflxLPR=self.CflxLPR_s,
            CflxLP=turnover["CflxLP"],
            CflxLD=turnover["CflxLD"],
            CflxSR=turnover["CflxSR"],
            CflxLUgrs=self.CflxLUgrs_s,
        )

        extra_var = {"CplsPLS": state["CplsP"] + state["CplsL"] + state["CplsS"]}
        extra_var["CflxNetPLS"] = Q(np.gradient(extra_var["CplsPLS"].m), "GtC/yr")

        return make_dataset_from_var_dict({**state, **flux, **extra_var}, time)


@define
class CarbonCycleModel:
    """ODE model for CNit's three terrestrial carbon pools.

    The model represents carbon in plant biomass (CplsP), litter (CplsL), and soil
    organic matter (CplsS). Net primary production adds carbon to the system; litter
    production respiration and soil respiration release carbon, while litter
    production and decomposition transfer carbon among pools. Land-use emissions
    remove carbon from the pools according to the configured allocation fractions.
    Temperature, nitrogen availability, and land-use history modify the relevant
    process rates.

    The model solves pool dynamics from carbon flux and effect series supplied by the
    coupled model. Evenly spaced subannual time axes are supported; the default
    examples use annual timesteps, with allocation fractions applying within each
    timestep. See :meth:`run` for the required series and simulation behavior.
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
    (time0).
    """

    CplsL0: pint.Quantity = field(
        default=Q(50, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial litter carbon pool size [GtC].

    Carbon contained in dead plant material (litter) at the start of the simulation
    (time0).
    """

    CplsS0: pint.Quantity = field(
        default=Q(1550, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial soil carbon pool size [GtC].

    Carbon contained in soil organic matter at the start of the simulation (time0).
    """

    time0: pint.Quantity = field(
        default=Q(1850, "yr"),
        validator=check_units("yr"),
    )
    """Initialization time [yr].

    Time point at which the initial pool sizes apply. The simulation time_axis must
    start at or after this time. The default is the pre-industrial baseline year
    (1850).
    """

    switch_Cpls: list[int] = field(
        factory=lambda: [1, 1, 1],
    )
    """Carbon pool switches for [plant, litter, soil] [dimensionless].

    Use 1 to enable or 0 to freeze a pool at its initial value. Other numeric values
    scale that pool's derivative. These switches are useful for diagnostic purposes
    or simplified model configurations.

    - switch_Cpls[0]: Plant pool switch
    - switch_Cpls[1]: Litter pool switch
    - switch_Cpls[2]: Soil pool switch
    """

    def run(
        self,
        time_axis: pint.Quantity,
        CflxNPP_s: pint.Quantity,
        CflxLPR_s: pint.Quantity,
        CflxLUgrs_s: pint.Quantity,
        eff_dT_cLP_s: pint.Quantity,
        eff_dT_cLD_s: pint.Quantity,
        eff_dT_cSR_s: pint.Quantity,
        eff_N_cLP_s: pint.Quantity,
        eff_N_cLD_s: pint.Quantity,
        eff_N_cSR_s: pint.Quantity,
        eff_LU_cLP_s: pint.Quantity,
        eff_LU_cLD_s: pint.Quantity,
        eff_LU_cSR_s: pint.Quantity,
    ) -> CarbonCycleModelResult:
        """
        Run the carbon cycle model simulation.

        This method solves the system of ODEs for carbon pool dynamics over the
        specified time period using scipy.integrate.solve_ivp. Numeric interpolators
        are built directly from the unit-aware flux and effect series for the
        solver's adaptive time steps.

        Parameters
        ----------
        time_axis
            Time points for simulation [yr].  Must start at or after time0.
            These are the time points at which solution is explicitly evaluated.
        CflxNPP_s, CflxLPR_s, CflxLUgrs_s
            Carbon flux series [GtC/yr], aligned with ``time_axis``.
        eff_*_s
            Dimensionless effect series aligned with ``time_axis``.

        Returns
        -------
        CarbonCycleModelResult
            Simulation results containing:

            - Carbon pool time series (CplsP, CplsL, CplsS)
            - References to calculator instances
            - All flux and effect series
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
        - Absolute tolerance: 1e-6 GtC
        - Relative tolerance: 1e-3
        - Default solver method (typically RK45)

        **Pool Switches:**

        Pool dynamics can be selectively disabled using switch_Cpls. When a switch
        is set to 0, the corresponding pool's rate of change is forced to zero,
        effectively freezing that pool at its initial value throughout the simulation.

        The flux and effect series are converted to numeric interpolators once before
        integration, allowing adaptive solver evaluations between output times.
        """

        t_eval = time_axis.to("yr").m
        # Pack the fluxes and environmental effect factors into a unified interpolator
        flux_and_effect_t = _make_packed_numeric_interp(
            [
                CflxNPP_s,
                CflxLPR_s,
                CflxLUgrs_s,
                eff_dT_cLP_s,
                eff_dT_cLD_s,
                eff_dT_cSR_s,
                eff_N_cLP_s,
                eff_N_cLD_s,
                eff_N_cSR_s,
                eff_LU_cLP_s,
                eff_LU_cLD_s,
                eff_LU_cSR_s,
            ],
            t_eval,
            [
                "GtC/yr",
                "GtC/yr",
                "GtC/yr",
                "1",
                "1",
                "1",
                "1",
                "1",
                "1",
                "1",
                "1",
                "1",
            ],
        )

        calc_turnover = self.calc_c_turnover
        calc_c_cycle = self.calc_c_cycle

        tau_CplsP_yr = calc_turnover.tau_CplsP.to("yr").m
        tau_CplsL_yr = calc_turnover.tau_CplsL.to("yr").m
        tau_CplsS_yr = calc_turnover.tau_CplsS.to("yr").m
        frc_NPP2P = calc_c_cycle.frc_NPP2P.to("1").m
        frc_NPP2L = calc_c_cycle.frc_NPP2L.to("1").m
        frc_cLP2L = calc_c_cycle.frc_cLP2L.to("1").m
        frc_cLD2S = calc_c_cycle.frc_cLD2S.to("1").m
        frc_cLUgrs2P = calc_c_cycle.frc_cLUgrs2P.to("1").m
        frc_cLUgrs2L = calc_c_cycle.frc_cLUgrs2L.to("1").m

        pool_switches = np.asarray(self.switch_Cpls, dtype=float)

        def func_to_solve(t: float, y: np.ndarray) -> np.ndarray:
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
            # Interpolate all fluxes and effect multipliers at continuous time t
            fluxes_and_effects = flux_and_effect_t(t)

            # Unpack carbon fluxes [GtC/yr]
            CflxNPP_GtC_yr = fluxes_and_effects[0]
            CflxLPR_GtC_yr = fluxes_and_effects[1]
            CflxLUgrs_GtC_yr = fluxes_and_effects[2]

            # Unpack environmental effects [1]
            eff_dT_cLP = fluxes_and_effects[3]
            eff_dT_cLD = fluxes_and_effects[4]
            eff_dT_cSR = fluxes_and_effects[5]
            eff_N_cLP = fluxes_and_effects[6]
            eff_N_cLD = fluxes_and_effects[7]
            eff_N_cSR = fluxes_and_effects[8]
            eff_LU_cLP = fluxes_and_effects[9]
            eff_LU_cLD = fluxes_and_effects[10]
            eff_LU_cSR = fluxes_and_effects[11]

            turnover = calc_turnover.calculate_numeric(
                CplsP_GtC=y[0],
                CplsL_GtC=y[1],
                CplsS_GtC=y[2],
                tau_CplsP_yr=tau_CplsP_yr,
                tau_CplsL_yr=tau_CplsL_yr,
                tau_CplsS_yr=tau_CplsS_yr,
                eff_dT_cLP=eff_dT_cLP,
                eff_dT_cLD=eff_dT_cLD,
                eff_dT_cSR=eff_dT_cSR,
                eff_N_cLP=eff_N_cLP,
                eff_N_cLD=eff_N_cLD,
                eff_N_cSR=eff_N_cSR,
                eff_LU_cLP=eff_LU_cLP,
                eff_LU_cLD=eff_LU_cLD,
                eff_LU_cSR=eff_LU_cSR,
            )

            dydt = (
                calc_c_cycle.calculate_dCpls_dt_numeric(
                    CflxNPP_GtC_yr=CflxNPP_GtC_yr,
                    CflxLPR_GtC_yr=CflxLPR_GtC_yr,
                    CflxLP_GtC_yr=turnover["CflxLP"],
                    CflxLD_GtC_yr=turnover["CflxLD"],
                    CflxSR_GtC_yr=turnover["CflxSR"],
                    CflxLUgrs_GtC_yr=CflxLUgrs_GtC_yr,
                    frc_NPP2P=frc_NPP2P,
                    frc_NPP2L=frc_NPP2L,
                    frc_cLP2L=frc_cLP2L,
                    frc_cLD2S=frc_cLD2S,
                    frc_cLUgrs2P=frc_cLUgrs2P,
                    frc_cLUgrs2L=frc_cLUgrs2L,
                )
                * pool_switches
            )
            if np.all(np.abs(dydt) < 1e-12):
                dydt[:] = 0
            return dydt

        if time_axis[0] < self.time0:
            raise ValueError(
                f"time_axis starts before time0: {time_axis[0]} < {self.time0}"
            )

        # Prepare solver arguments
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
            CflxNPP_s=CflxNPP_s,
            CflxLPR_s=CflxLPR_s,
            CflxLUgrs_s=CflxLUgrs_s,
            eff_dT_cLP_s=eff_dT_cLP_s,
            eff_dT_cLD_s=eff_dT_cLD_s,
            eff_dT_cSR_s=eff_dT_cSR_s,
            eff_N_cLP_s=eff_N_cLP_s,
            eff_N_cLD_s=eff_N_cLD_s,
            eff_N_cSR_s=eff_N_cSR_s,
            eff_LU_cLP_s=eff_LU_cLP_s,
            eff_LU_cLD_s=eff_LU_cLD_s,
            eff_LU_cSR_s=eff_LU_cSR_s,
        )
