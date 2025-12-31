"""
Land Use Flux Calculator for the CNit Model.

This module defines the :class:`LanduseCalculator` class, which processes land use
change forcing data and translates it into the specific fluxes required by the CNit
model.  The calculator handles both carbon and nitrogen land use fluxes with explicit
separation of deforestation, afforestation (or reforestation), and regrowth dynamics.

**Key Features:**

1. **Explicit Separation of Deforestation and Afforestation:**
   Land use gross emissions (LUgrs) are separated into deforestation (LUdfst, always
   positive) and afforestation/reforestation (LUafst, always positive). This separation
   is critical for scenarios involving carbon dioxide removal (CDR) and net-zero targets
   where afforestation plays a significant role.

2. **Regrowth Following Deforestation:**
   Natural recovery after deforestation is represented through partial regrowth (LUrgr),
   calculated as a fraction of deforestation distributed over a regrowth timescale.

3. **Decay Following Afforestation:**
   Carbon sequestered through afforestation may not remain permanently due to vegetation
   turnover. The calculator includes an afforestation decay flux (LUafst_decay) that
   represents the gradual return of a fraction of sequestered carbon.

4. **Flexible Input Methods:**
   Supports two calculation methods:

   - Method 0: Uses net emissions (LUnet) as input (e.g., AFOLU emissions)
   - Method 1: Uses gross emissions (LUgrs) as input (e.g., harvest/deforestation data)

**Land Use Flux Relationships:**

The net land use emission is related to gross emission and regrowth by:

..  math::

    LUnet(t) = LUgrs(t) - LUrgr(t)

Gross emissions are further separated into deforestation and afforestation:

..  math::

    LUnet(t) = LUdfst(t) - LUrgr(t) - LUafst(t)

where:

- :math:`LUdfst(t) = \\max(LUgrs(t), 0)` (deforestation, positive)
- :math:`LUafst(t) = \\max(-LUgrs(t), 0)` (afforestation, positive)

**Regrowth Dynamics:**

For a deforestation event at time t, regrowth occurs as a constant flux over the
regrowth period:

..  math::

    LUrgr(T) = \\frac{\\varphi \\times LUdfst(t)}{\\tau_{rgr}},
    \\quad \\forall T \\in [t+1, t+\\tau_{rgr}]

where :math:`\\varphi` is the fraction that can regrow and :math:`\\tau_{rgr}` is
the regrowth timescale.

Total regrowth at time t is the sum of contributions from all past deforestation
events still undergoing recovery:

.. math::

    LUrgr(t) = \\begin{cases}
        \\sum_{T = \\max(0, t - \\tau_{rgr})}^{t-1} \\frac{\\varphi \\times LUdfst(T)}{\\tau_{rgr}}
        & \\text{if } t \\geq 1 \\\\
        0 & \\text{if } t < 1
    \\end{cases}

**Afforestation Decay:**

Similarly, for afforestation at time t, a fraction of sequestered carbon decays over
time:

.. math::

    LUafst\_decay(T) = \\frac{\\phi \\times LUafst(t)}{\\tau_{afst\_decay}},
    \\quad \\forall T \\in [t+1, t+\\tau_{afst\_decay}]

Total decay at time t sums all contributions from past afforestation events:

.. math::

    LUafst\_decay(t) = \\begin{cases}
        \\sum_{T = \\max(0, t - \\tau_{afst\_decay})}^{t-1} \\frac{\\phi \\times LUafst(T)}{\\tau_{afst\_decay}}
        & \\text{if } t \\geq 1 \\\\
        0 & \\text{if } t < 1
    \\end{cases}

where :math:`\\phi` is the afforestation decay fraction (default :math:`\\phi = 1 - \\varphi`) and
:math:`\\tau_{afst\_decay}` is the decay timescale (default: :math:`\\tau_{afst\_decay} = \\tau_{rgr}`).

**Calculation Algorithms:**

The calculator uses two different algorithms depending on the input method:

**Method 0: Net Emission as Input**

When using net emissions as input (suitable for AFOLU emissions data), the algorithm
iteratively calculates gross flux by adding regrowth, then propagates regrowth and
decay effects forward in time:

1. Initialize first timestep with zero regrowth
2. For each timestep:

   a. Calculate gross flux:  :math:`LUgrs(t) = LUnet(t) + LUrgr(t)`
   b. If :math:`LUgrs(t) > 0`: Deforestation event

      - Set :math:`LUdfst(t) = LUgrs(t)`
      - Add regrowth contribution to future timesteps:
        :math:`LUrgr(T) += \\frac{\\varphi \\times LUdfst(t)}{\\tau_{rgr}}`
        for :math:`T \\in [t+1, t+\\tau_{rgr}]`

   c. If :math:`LUgrs(t) < 0`: Afforestation event

      - Set :math:`LUafst(t) = -LUgrs(t)`
      - Add decay contribution to future timesteps:
        :math:`LUafst\_decay(T) += \\frac{\\phi \\times LUafst(t)}{\\tau_{afst\_decay}}`
        for :math:`T \\in [t+1, t+\\tau_{afst\_decay}]`

3. Regrowth and decay contributions accumulate for all affected future timesteps

This method ensures that when net emissions are negative (indicating net sequestration
beyond regrowth), the resulting negative gross flux is properly interpreted as
afforestation rather than negative deforestation.

**Method 1: Gross Emission as Input**

When using gross emissions directly (suitable for harvest/deforestation data), the
algorithm uses the provided gross flux and calculates regrowth and decay by propagating
effects forward in time:

1. For each timestep with gross flux:

   a. If :math:`LUgrs(t) > 0`: Deforestation event

      - Set :math:`LUdfst(t) = LUgrs(t)`
      - Add regrowth contribution to future timesteps:
        :math:`LUrgr(T) += \\frac{\\varphi \\times LUdfst(t)}{\\tau_{rgr}}`
        for :math:`T \\in [t+1, t+\\tau_{rgr}]`

   b. If :math:`LUgrs(t) < 0`: Afforestation event

      - Set :math:`LUafst(t) = -LUgrs(t)`
      - Add decay contribution to future timesteps:
        :math:`LUafst\_decay(T) += \\frac{\\phi \\times LUafst(t)}{\\tau_{afst\_decay}}`
        for :math:`T \\in [t+1, t+\\tau_{afst\_decay}]`

2. Calculate net flux:  :math:`LUnet(t) = LUgrs(t) - LUrgr(t)`

This method is more straightforward when gross emissions are directly observed from
data sources like forest harvest statistics.

**Important Notes:**

- :math:`LUafst\_decay` is not part of the net land use flux equations.     It provides
  a reference for afforestation-related decay independent of environmental changes,
  analogous to :   math:`LUrgr`.
- :math:`LUrgr` is explicitly part of NPP and is further modified by CO₂
  fertilization and climate change.
- Decay of afforested carbon is implicitly included in vegetation turnover and is
  therefore already affected by environmental conditions.
- :math:`LUafst\_decay` supports parameterization of NPP effects associated with
  afforestation and deforestation.
- Both algorithms use in-place array operations for efficiency when processing long
  time series.
- All arrays are modified in-place to minimize memory allocation.
- The algorithms properly handle edge cases including zero fluxes and boundary effects
  at the start and end of time series.

**Typical Usage:**

..    code-block:: python

    # Initialize calculator
    lu_calc = LanduseCalculator(
        frc_rgr=Q(0.9, "1"),           # varphi = 0.9 (90% can regrow)
        tau_rgr=Q(100, "yr"),          # 100-year regrowth period
        method_cLU=Q(0, "1"),          # Use net carbon emissions
        method_nLU=Q(1, "1"),          # Use gross nitrogen emissions
    )

    # Calculate carbon fluxes
    carbon_fluxes = lu_calc.calculate_CflxLU_series(
        CemsLUnet_s=net_emissions,
        CemsLUgrs_s=gross_emissions,
    )

    # Access results
    deforestation = carbon_fluxes["CflxLUdfst"]
    afforestation = carbon_fluxes["CflxLUafst"]
    regrowth = carbon_fluxes["CflxLUrgr"]
    cumulative_regrowth = carbon_fluxes["cumsum_CflxLUrgr"]

See Also
--------
:py:meth:`cnit.physics.carbon_nitrogen_cycle.CarbonNitrogenCycleModel`: Main model using these land use fluxes
:py:meth:`cnit.physics.effects.EffectLanduseCalculator`: Environmental effects from land use changes
"""

from typing import Dict
from attrs import define, field
import numpy as np
import pint

from ..utils.units import check_units, Q


@define
class LanduseCalculator:
    """Calculator for land use change fluxes with explicit deforestation/afforestation separation.

    This calculator processes land use change forcing data (net or gross emissions) and
    generates detailed flux time series including deforestation, afforestation, regrowth,
    and decay components. The separation of deforestation and afforestation is essential
    for scenarios involving carbon dioxide removal (CDR) and net-zero targets.

    The calculator operates on annual time series and produces both instantaneous fluxes
    and cumulative quantities useful for tracking land use history and calculating
    environmental effects.
    """

    frc_rgr: pint.Quantity = field(
        default=Q(0.9, "1"),
        validator=check_units("1"),
    )
    """Fraction of deforestation that can regrow [dimensionless].  
    
    Represents the portion of deforested area that undergoes natural recovery.
    The remaining fraction (1 - frc_rgr) represents permanent land conversion with no
    natural regrowth.
    """

    tau_rgr: pint. Quantity = field(
        default=Q(100, "yr"),
        validator=check_units("yr"),
    )
    """Regrowth timescale [yr].
    
    Time period over which regrowth occurs following deforestation. 
    Regrowth is distributed as a constant flux over this period.
    """

    method_cLU:  pint.Quantity = field(
        default=Q(0, "1"),
        validator=check_units("1"),
    )
    """Carbon land use flux calculation method [dimensionless].
    
    Selects the input data type and calculation approach: 
    
    - 0: Use net emissions (CemsLUnet) as input, suitable for AFOLU emissions data
      where regrowth must be calculated from net emissions
    - 1: Use gross emissions (CemsLUgrs) as input, suitable for harvest/deforestation
      data where gross emissions are directly observed
    """

    method_nLU: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """Nitrogen land use flux calculation method [dimensionless].
    
    Selects the input data type and calculation approach:
    
    - 0: Use net emissions (NemsLUnet) as input (less commonly used)
    - 1: Use gross emissions (NemsLUgrs) as input, suitable for gross deforestation
      data (typical for nitrogen)
    """

    def calculate_CflxLU_series(
            self,
            CemsLUnet_s: pint.Quantity,
            CemsLUgrs_s: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate carbon flux time series from land use change forcing.

        This method processes land use change emissions and generates comprehensive
        flux time series including gross emissions, regrowth, deforestation,
        afforestation, and decay components. The method automatically separates
        deforestation (positive gross flux) from afforestation (negative gross flux).

        Parameters
        ----------
        CemsLUnet_s
            Net land use carbon emission time series [GtC/yr].  Represents the
            difference between gross emissions and regrowth.  Used when method_cLU=0.
        CemsLUgrs_s
            Gross land use carbon emission time series [GtC/yr].  Can be positive
            (deforestation) or negative (afforestation). Used when method_cLU=1.

        Returns
        -------
        Dict[str, pint.Quantity]
            Dictionary containing carbon flux time series:

            **Instantaneous Fluxes [GtC/yr]:**

            - CflxLUgrs: Gross land use flux (can be positive or negative)
            - CflxLUrgr: Regrowth flux following deforestation (always non-negative)
            - CflxLUdfst: Deforestation flux (always non-negative)
            - CflxLUafst: Afforestation flux (always non-negative)
            - CflxLUafst_decay: Decay flux following afforestation (always non-negative)

            **Cumulative Quantities [GtC]:**

            - cumsum_CflxLUrgr:   Cumulative regrowth carbon
            - cumsum_CflxLUdfst: Cumulative deforestation carbon
            - cumsum_CflxLUafst: Cumulative afforestation carbon
            - cumsum_CflxLUafst_decay: Cumulative decay from afforestation

        Notes
        -----
        **Flux Relationships:**

        Net emissions relate to gross emissions and regrowth:

        ..  math::

            LUnet(t) = LUgrs(t) - LUrgr(t)

        Gross emissions are separated into components:

        .. math::

            LUnet(t) = LUdfst(t) - LUrgr(t) - LUafst(t)

        where deforestation and afforestation are always positive:

        .. math::

            LUdfst(t) = \\max(LUgrs(t), 0)

            LUafst(t) = \\max(-LUgrs(t), 0)

        **Regrowth Calculation:**

        Regrowth following deforestation at time t:

        .. math::

            LUrgr(T) = \\frac{\\varphi \\times LUdfst(t)}{\\tau_{rgr}},
            \\quad T \\in [t+1, t+\\tau_{rgr}]

        Total regrowth is sum of all active regrowth contributions:

        .. math::

            LUrgr(t) = \\sum_{T=\\max(0, t-\\tau_{rgr})}^{t-1}
            \\frac{\\varphi \\times LUdfst(T)}{\\tau_{rgr}}

        **Afforestation Decay:**

        Decay following afforestation at time t:

        .. math::

            LUafst\_decay(T) = \\frac{\\phi \\times LUafst(t)}{\\tau_{afst\_decay}},
            \\quad T \\in [t+1, t+\\tau_{afst\_decay}]

        where by default :math:`\\phi = 1 - \\varphi` and :math:`\\tau_{afst\_decay} =
        \\tau_{rgr}`.

        **Method Selection:**

        - Method 0: Calculates LUgrs from LUnet by adding calculated regrowth.
          Suitable when net emissions are the primary data source (e.g., AFOLU).

        - Method 1: Uses LUgrs directly and calculates regrowth from it.
          Suitable when gross emissions are observed directly (e.g., harvest data).

        **Cumulative Quantities:**

        Cumulative values track the total historical land use changes and are used
        to calculate land use effects on ecosystem productivity (see
        :class:`EffectLanduseCalculator`).

        See Also
        --------
        calculate_NflxLU_series: Analogous nitrogen flux calculation
        :py:meth:`cnit.physics.effects.EffectLanduseCalculator`: Land use effects calculation using these fluxes
        :py:meth:`cnit.physics.carbon_cycle.CarbonCycleModel.run`: How these fluxes are used in carbon cycle
        """

        return self._calculate_flux_series(
            net_flux=CemsLUnet_s,
            gross_flux=CemsLUgrs_s,
            method_calc=self.method_cLU,
            unit="GtC/yr",
            prefix="Cflx",
        )

    def calculate_NflxLU_series(
            self,
            NemsLUnet_s: pint.Quantity,
            NemsLUgrs_s: pint.Quantity,
            NemsLUmin_s: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate nitrogen flux time series from land use change forcing.

        This method processes land use change emissions and generates comprehensive
        nitrogen flux time series.    In addition to organic nitrogen fluxes (analogous
        to carbon), it includes direct mineral nitrogen emissions (e.g., NH₃, NOₓ,
        N₂O from agricultural land use change).

        Parameters
        ----------
        NemsLUnet_s
            Net land use nitrogen emission time series [GtN/yr].  Represents emissions
            from organic nitrogen pools.    Used when method_nLU=0.
        NemsLUgrs_s
            Gross land use nitrogen emission time series [GtN/yr].  Can be positive
            (deforestation) or negative (afforestation). Used when method_nLU=1.
        NemsLUmin_s
            Direct land use nitrogen emission from mineral pool time series [GtN/yr].
            Represents gaseous losses (NH₃, NOₓ, N₂O) during land conversion,
            primarily from agricultural activities.

        Returns
        -------
        Dictionary containing nitrogen flux time series:

            **Instantaneous Fluxes [GtN/yr]:**

            - NflxLUgrs:  Gross land use flux from organic pools
            - NflxLUrgr:  Regrowth flux following deforestation
            - NflxLUdfst:    Deforestation flux (from organic pools)
            - NflxLUafst:  Afforestation flux (to organic pools)
            - NflxLUafst_decay:    Decay flux following afforestation
            - NflxLUmin:   Direct emission from mineral pool

            **Cumulative Quantities [GtN]:**

            - cumsum_NflxLUrgr: Cumulative regrowth nitrogen
            - cumsum_NflxLUdfst: Cumulative deforestation nitrogen
            - cumsum_NflxLUafst:    Cumulative afforestation nitrogen
            - cumsum_NflxLUafst_decay:   Cumulative decay from afforestation

        Notes
        -----
        **Nitrogen-Specific Considerations:**

        1. **Organic vs Mineral Emissions:**

           - LUgrs, LUrgr, LUdfst, LUafst: Affect organic nitrogen pools (plant, litter, soil)
           - LUmin: Directly affects mineral nitrogen pool, bypassing organic pools

        2. **Mineral Nitrogen Emissions:**
           Direct mineral emissions (LUmin) typically include:

           - NH₃ volatilization from fertilizer application
           - NOₓ emissions from agricultural activities
           - N₂O emissions from soil disturbance

           These emissions do not go through organic matter turnover and are subtracted
           directly from the mineral nitrogen pool.

        3. **Stoichiometric Coupling:**
           Organic nitrogen land use fluxes (LUgrs) are coupled to carbon fluxes through
           ecosystem C:N ratios, but mineral emissions (LUmin) are independent.

        **Typical Nitrogen Method:**

        Method 1 is typically used for nitrogen because gross deforestation data
        combined with C:N ratios provides more direct estimates of organic nitrogen
        losses than net emission calculations.

        See Also
        --------
        calculate_CflxLU_series: Analogous carbon flux calculation
        :py:meth:`cnit.physics.nitrogen_cycle.NitrogenCycleModel.run`: How these fluxes are used in nitrogen cycle
        """
        return self._calculate_flux_series(
            net_flux=NemsLUnet_s,
            gross_flux=NemsLUgrs_s,
            method_calc=self.method_nLU,
            unit="GtN/yr",
            prefix="Nflx",
        ) | {"NflxLUmin": NemsLUmin_s}

    def _calculate_flux_series(
            self,
            net_flux: pint.Quantity,
            gross_flux: pint.Quantity,
            method_calc: pint.Quantity,
            unit: str,
            prefix: str,
    ) -> Dict[str, pint.Quantity]:
        """
        Core flux calculation engine used by both carbon and nitrogen methods.

        This internal method handles the common flux calculation logic, including
        regrowth propagation, afforestation decay, and cumulative tracking.   It
        operates on dimensionless arrays and adds units at the end.

        Parameters
        ----------
        net_flux
            Net emission time series.
        gross_flux
            Gross emission time series.
        method_calc
            Calculation method (0 or 1).
        unit
            Unit string for flux quantities (e.g., "GtC/yr" or "GtN/yr").
        prefix
            Prefix for output dictionary keys ("Cflx" or "Nflx").

        Returns
        -------
        Dict[str, pint.Quantity]
            Dictionary with named flux series (see public methods for details).

        Raises
        ------
        ValueError
            If method_calc is not 0 or 1.

        Notes
        -----
        This method implements the core algorithms described in the public methods.
        It uses in-place array operations for efficiency when processing long time
        series.

        The method handles edge cases:
        - Zero fluxes: Returns arrays of zeros
        - Boundary effects: Properly handles regrowth/decay at start and end of series
        - Numerical precision: Uses float64 for accumulation calculations

        **Parameter Relationships:**

        - frc_rgr: Regrowth fraction (varphi)
        - frc_afst_decay:  Afforestation decay fraction (phi = 1 - varphi)
        - tau_rgr_int: Regrowth timescale (tau_rgr)
        - tau_afst_decay_int: Afforestation decay timescale (tau_afst_decay = tau_rgr)
        """


        frc_rgr = self.frc_rgr.m
        tau_rgr_int = int(self.tau_rgr.m)
        frc_afst_decay = 1 - frc_rgr
        tau_afst_decay_int = tau_rgr_int

        length_input = len(net_flux.m)

        # Handle zero flux case
        if np.allclose(net_flux, 0) and np.allclose(gross_flux, 0):
            gross_flux_a = regrowth_flux_a = deforest_flux_a = afforest_flux_a = (
                afforest_decay_flux_a
            ) = np.zeros(length_input, dtype=float)
        else:
            net_flux_a = net_flux.astype(float).m
            gross_flux_a = gross_flux.astype(float).m
            regrowth_flux_a = np.zeros(length_input, dtype=float)
            deforest_flux_a = np.zeros(length_input, dtype=float)
            afforest_flux_a = np.zeros(length_input, dtype=float)
            afforest_decay_flux_a = np.zeros(length_input, dtype=float)

            if method_calc.m == 0:
                # Method 0: Use net flux as input (e.g., AFOLU emissions)
                self._calculate_from_net_flux(
                    net_flux_a,
                    gross_flux_a,
                    regrowth_flux_a,
                    deforest_flux_a,
                    afforest_flux_a,
                    afforest_decay_flux_a,
                    tau_rgr_int,
                    frc_rgr,
                    tau_afst_decay_int,
                    frc_afst_decay,
                )
            elif method_calc.m == 1:
                # Method 1: Use gross flux as input (e.g., fCharvest emissions)
                # When using gross flux directly, we need to copy the input values first
                self._calculate_from_gross_flux(
                    net_flux_a,
                    gross_flux_a,
                    regrowth_flux_a,
                    deforest_flux_a,
                    afforest_flux_a,
                    afforest_decay_flux_a,
                    tau_rgr_int,
                    frc_rgr,
                    tau_afst_decay_int,
                    frc_afst_decay,
                )
            else:
                raise ValueError(
                    f"Invalid method_calc value: {method_calc.m}. Must be 0 or 1."
                )

        # Calculate cumulative quantities
        non_time_unit = unit.replace("/yr", "")

        return {
            f"{prefix}LUgrs": Q(gross_flux_a, unit),
            f"{prefix}LUrgr": Q(regrowth_flux_a, unit),
            f"{prefix}LUdfst": Q(deforest_flux_a, unit),
            f"{prefix}LUafst": Q(afforest_flux_a, unit),
            f"{prefix}LUafst_decay": Q(afforest_decay_flux_a, unit),
            f"cumsum_{prefix}LUrgr": Q(np.cumsum(regrowth_flux_a), non_time_unit),
            f"cumsum_{prefix}LUdfst": Q(np.cumsum(deforest_flux_a), non_time_unit),
            f"cumsum_{prefix}LUafst": Q(np.cumsum(afforest_flux_a), non_time_unit),
            f"cumsum_{prefix}LUafst_decay": Q(
                np.cumsum(afforest_decay_flux_a), non_time_unit
            ),
        }

    def _calculate_from_net_flux(
            self,
            net_flux_values: np.ndarray,
            gross_flux_values: np.ndarray,
            regrowth_flux_values: np.ndarray,
            deforest_flux_values: np.ndarray,
            afforest_flux_values: np.ndarray,
            afforest_decay_flux_values: np.ndarray,
            tau_rgr_int: int,
            frc_rgr: float,
            tau_afst_decay_int: int,
            frc_afst_decay: float,
    ) -> None:
        """
        Calculate flux values when using net flux as input.
        Updates the gross_flux_values array and regrowth_flux_values array in-place.

        Parameters
        ----------
        net_flux_values : array-like
            Array of net flux values.
        gross_flux_values : array-like
            Array to store gross flux values (modified in-place).
        regrowth_flux_values : array-like
            Array to store regrowth flux values (modified in-place).
        length : int
            Length of the time series.
        tau_rgr_int : int
            Integer regrowth time constant.
        frc_rgr : float
            Fraction that can regrow.
        """

        length_input = len(net_flux_values)
        # Initialize first time step regrowth_flux_value
        regrowth_flux_values[0] = 0

        # Process the rest of the time steps
        for i in range(0, length_input - 1):
            # Current gross flux equals net flux plus regrowth
            gross_flux_values[i] = net_flux_values[i] + regrowth_flux_values[i]

            # Calculate regrowth contribution for future time steps
            # Regrowth occurs over the next tau_rgr_int steps only when the gross
            # flux is positive
            if gross_flux_values[i] > 0:
                deforest_flux_values[i] = gross_flux_values[i]
                regrowth_per_step = gross_flux_values[i] * frc_rgr / tau_rgr_int
                regrowth_flux_values[
                i + 1: min(i + 1 + tau_rgr_int, length_input)
                ] += regrowth_per_step
            else:
                afforest_flux_values[i] = -gross_flux_values[i]
                afforest_decay_per_step = (
                        afforest_flux_values[i] * frc_afst_decay / tau_afst_decay_int
                )
                afforest_decay_flux_values[
                i + 1: min(i + 1 + tau_afst_decay_int, length_input)
                ] += afforest_decay_per_step

        # Final step
        gross_flux_values[-1] = net_flux_values[-1] + regrowth_flux_values[-1]
        if gross_flux_values[-1] > 0:
            deforest_flux_values[-1] = gross_flux_values[-1]
        else:
            afforest_flux_values[-1] = -gross_flux_values[-1]

    def _calculate_from_gross_flux(
            self,
            net_flux_values: np.ndarray,
            gross_flux_values: np.ndarray,
            regrowth_flux_values: np.ndarray,
            deforest_flux_values: np.ndarray,
            afforest_flux_values: np.ndarray,
            afforest_decay_flux_values: np.ndarray,
            tau_rgr_int: int,
            frc_rgr: float,
            tau_afst_decay_int: int,
            frc_afst_decay: float,
    ) -> None:
        """
        Calculate flux values when using gross flux as input.
        Updates the regrowth_flux_values array in-place.

        Parameters
        ----------
        gross_flux_values : array-like
            Array of gross flux values.
        regrowth_flux_values : array-like
            Array to store regrowth flux values (modified in-place).
        length : int
            Length of the time series.
        tau_rgr_int : int
            Integer regrowth time constant.
        frc_rgr : float
            Fraction that can regrow.
        """
        length_input = len(gross_flux_values)
        # Calculate regrowth for each time step
        for i in range(length_input):
            if gross_flux_values[i] > 0:
                deforest_flux_values[i] = gross_flux_values[i]
                regrowth_per_step = deforest_flux_values[i] * frc_rgr / tau_rgr_int
                regrowth_flux_values[
                i + 1: min(i + 1 + tau_rgr_int, length_input)
                ] += regrowth_per_step
            else:
                afforest_flux_values[i] = -gross_flux_values[i]
                afforest_decay_per_step = (
                        afforest_flux_values[i] * frc_afst_decay / tau_afst_decay_int
                )
                afforest_decay_flux_values[
                i + 1: min(i + 1 + tau_afst_decay_int, length_input)
                ] += afforest_decay_per_step

        # Calculate net flux as gross flux minus regrowth
        net_flux_values[:] = gross_flux_values - regrowth_flux_values
