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

For an interval with positive gross land use flux rate :math:`F_i` and duration
:math:`\\Delta t`, the deforestation amount is :math:`D_i = F_i \\Delta t`.
Regrowth from that amount starts at the following interval and occurs at a constant
rate over the configured regrowth timescale:

..  math::

    LUrgr(T) = \\frac{\\varphi \\times D_i}{\\tau_{rgr}},
    \\quad T \\in [t_{i+1}, t_{i+1}+\\tau_{rgr})

where :math:`\\varphi` is the fraction that can regrow.

In each interval, total regrowth is the sum of contributions from all past deforestation
events whose regrowth periods are active:

..  math::

    LUrgr(t) = \\sum_i \\frac{\\varphi D_i}{\\tau_{rgr}}
    \\mathbb{1}_{[t_{i+1}, t_{i+1}+\\tau_{rgr})}(t)

**Afforestation Decay:**

Similarly, a fraction of each interval's afforestation amount decays at a constant
rate over the afforestation-decay timescale, beginning at the next interval boundary.

..  math::

    LUafst\_decay(T) = \\frac{\\phi A_i}{\\tau_{afst\_decay}},
    \\quad T \\in [t_{i+1}, t_{i+1}+\\tau_{afst\_decay})

where :math:`\\phi` is the afforestation decay fraction (default
:math:`\\phi = 1 - \\varphi`). In the current calculator,
:math:`\\tau_{afst\_decay} = \\tau_{rgr}`.

**Implementation note:** The calculation rounds each configured duration to the
nearest whole number of model intervals, with a minimum of one, and spreads the full
regrowth or afforestation-decay amount evenly across those intervals. Using whole
intervals keeps the regrowth and afforestation-decay amounts straightforward to track
in the model's timestep-based budgets. For example, with annual timesteps, a configured
timescale of 100.5 years would end halfway through the 101st year. Annual output has
one flux value for that whole year, so it cannot clearly show that regrowth or
afforestation decay acts for only half of it. Using a whole number of annual
intervals keeps the reported rate aligned with the timestep. The implemented rates
can therefore differ slightly from the rates in the equations above.

**Calculation Algorithms:**

The calculator uses two different algorithms depending on the input method:

**Method 0: Net Emission as Input**

When using net emissions as input (suitable for AFOLU emissions data), the algorithm
iteratively calculates gross flux by adding regrowth, then propagates regrowth and
afforestation-decay effects forward in time:

1. Initialize first timestep with zero regrowth
2. For each timestep:

   a. Calculate gross flux: :math:`LUgrs(t) = LUnet(t) + LUrgr(t)`
   b. If :math:`LUgrs(t) > 0`: Deforestation event

      - Set :math:`LUdfst(t) = LUgrs(t)`
      - Schedule regrowth over ``tau_rgr`` years, rounding the duration to the nearest
        whole number of model intervals (at least one).

   c. If :math:`LUgrs(t) < 0`: Afforestation event

      - Set :math:`LUafst(t) = -LUgrs(t)`
      - Schedule afforestation decay over its rounded duration, based on the
        afforestation amount in the current interval.

3. Record each contribution's rate change at its start and end intervals. During the
   time loop, accumulate these changes to obtain the active regrowth and
   afforestation-decay rates. This gives the same rates as adding each contribution to
   every affected future interval, while avoiding repeated updates across those
   intervals.

This method ensures that when net emissions are negative (indicating net sequestration
beyond regrowth), the resulting negative gross flux is properly interpreted as
afforestation rather than negative deforestation.

**Method 1: Gross Emission as Input**

When using gross emissions directly (suitable for harvest/deforestation data), the
algorithm uses the provided gross flux and calculates regrowth and afforestation decay
by propagating effects forward in time:

1. For each timestep with gross flux:

   a. If :math:`LUgrs(t) > 0`: Deforestation event

      - Set :math:`LUdfst(t) = LUgrs(t)`
      - Schedule regrowth over ``tau_rgr`` years, rounding the duration to the nearest
        whole number of model intervals (at least one).

   b. If :math:`LUgrs(t) < 0`: Afforestation event

      - Set :math:`LUafst(t) = -LUgrs(t)`
      - Schedule afforestation decay over its rounded duration, based on the
        afforestation amount in the current interval.

2. Calculate net flux: :math:`LUnet(t) = LUgrs(t) - LUrgr(t)`

The model uses range-difference scheduling for each contribution. For example, a
contribution of 6 units of rate in intervals 1, 2, and 3 is recorded as ``+6`` at
interval 1 and ``-6`` at interval 4. The changes ``[0, 6, 0, 0, -6]`` have the running
sum ``[0, 6, 6, 6, 0]``, which gives the rate in each interval. Overlapping
contributions add their start and stop changes to the same array. This takes constant
time per event and one pass over the time series, rather than repeatedly adding rates
to every affected future interval. ``_add_future_pulse`` remains as a direct,
interval-by-interval reference implementation of the same calculation.

This method is more straightforward when gross emissions are directly observed from
data sources like forest harvest statistics.

**Important Notes:**

- :math:`LUafst\_decay` is not part of the net land use flux equations. It provides
  a reference for afforestation-related decay independent of environmental changes,
  analogous to :math:`LUrgr`.
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

from typing import Dict, Optional
from attrs import define, field
import numpy as np
import pint

from ..utils.units import check_units, Q


@define
class LanduseCalculator:
    """Calculator for land use change fluxes with explicit deforestation/afforestation separation.

    This calculator processes land use change forcing data (net or gross emissions) and
    generates detailed flux time series including deforestation, afforestation,
    regrowth, and afforestation-decay components. The separation of deforestation and
    afforestation is essential for scenarios involving carbon dioxide removal (CDR)
    and net-zero targets.

    The calculator operates on evenly spaced time series aligned with ``time_axis``.
    Each input rate is held constant over its following interval. Regrowth and
    afforestation decay are scheduled over ``tau_rgr`` years, with the duration
    rounded to the nearest whole number of model intervals (at least one). If no time
    axis is supplied, samples are treated as annual. Experiment configs validate input
    units, series lengths, and time-axis regularity before model execution. Direct calls
    to this calculator assume those inputs are valid and aligned. Direct callers must also
    provide ``0 <= frc_rgr <= 1`` and a positive ``tau_rgr``.
    """

    frc_rgr: pint.Quantity = field(
        default=Q(0.9, "1"),
        validator=check_units("1"),
    )
    """Fraction of deforestation that can regrow [dimensionless].

    Represents the portion of deforested area that can regrow naturally.
    The remaining fraction (1 - frc_rgr) represents permanent land conversion with no
    natural regrowth.
    """

    tau_rgr: pint.Quantity = field(
        default=Q(100, "yr"),
        validator=check_units("yr"),
    )
    """Regrowth timescale [yr].

    Time period over which regrowth occurs following deforestation.
    Regrowth is distributed as a constant flux over this period.
    """

    method_cLU: pint.Quantity = field(
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
        time_axis: Optional[pint.Quantity] = None,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate carbon flux time series from land use change forcing.

        This method processes land use change emissions and generates comprehensive
        flux time series including gross emissions, regrowth, deforestation,
        afforestation, and afforestation-decay components. The method automatically separates
        deforestation (positive gross flux) from afforestation (negative gross flux).

        Parameters
        ----------
        CemsLUnet_s
            Net land use carbon emission time series [GtC/yr].  Represents the
            difference between gross emissions and regrowth.  Used when method_cLU=0.
        CemsLUgrs_s
            Gross land use carbon emission time series [GtC/yr].  Can be positive
            (deforestation) or negative (afforestation). Used when method_cLU=1.
        time_axis
            Time coordinates aligned with the input series [yr]. They must be
            one-dimensional, strictly increasing, evenly spaced, and have the same
            length as the flux series. Each coordinate marks the start of the
            corresponding flux interval. Thus, the last value is also applied over
            one interval after the final coordinate. These conditions are checked by
            the experiment config when used through the model. Direct callers must
            supply an axis satisfying these requirements. If omitted, each value
            represents a one-year interval.

        Returns
        -------
        Dict[str, pint.Quantity]
            Dictionary containing carbon flux time series:

            **Instantaneous Fluxes [GtC/yr]:**

            - CflxLUgrs: Gross land use flux (can be positive or negative)
            - CflxLUrgr: Regrowth flux following deforestation (always non-negative)
            - CflxLUdfst: Deforestation flux (always non-negative)
            - CflxLUafst: Afforestation flux (always non-negative)
            - CflxLUafst_decay: Afforestation decay flux (always non-negative)

            **Cumulative Quantities [GtC]:**

            - cumsum_CflxLUrgr: Cumulative regrowth carbon
            - cumsum_CflxLUdfst: Cumulative deforestation carbon
            - cumsum_CflxLUafst: Cumulative afforestation carbon
            - cumsum_CflxLUafst_decay: Cumulative afforestation decay

        Notes
        -----
        **Flux Relationships:**

        Net emissions relate to gross emissions and regrowth:

        .. math::

            LUnet(t) = LUgrs(t) - LUrgr(t)

        Gross emissions are separated into components:

        .. math::

            LUnet(t) = LUdfst(t) - LUrgr(t) - LUafst(t)

        where deforestation and afforestation are always positive:

        .. math::

            LUdfst(t) = \\max(LUgrs(t), 0)

            LUafst(t) = \\max(-LUgrs(t), 0)

        **Regrowth Calculation:**

        Each positive gross flux rate is multiplied by the interval duration to get
        the amount of deforestation during that interval. The fraction ``frc_rgr``
        of this amount is returned at a constant rate over the configured
        ``tau_rgr`` timescale, starting at the beginning of the next interval. The
        duration is rounded to the nearest whole number of model intervals using
        ``max(1, round(tau_rgr / time_step))``; the full regrowth amount is spread
        over those intervals. Regrowth beyond the end of the supplied series is
        outside the returned result.

        **Afforestation Decay:**

        For each negative gross flux rate, the absolute interval amount is multiplied
        by ``1 - frc_rgr``. This afforestation-decay amount is spread uniformly over
        the same rounded ``tau_rgr`` duration as regrowth, starting at the next
        interval boundary.

        where :math:`\\phi = 1 - \\varphi` by default. The current calculator
        sets :math:`\\tau_{afst\_decay} = \\tau_{rgr}`.

        **Method Selection:**

        - Method 0: Calculates LUgrs from LUnet by adding calculated regrowth.
          Suitable when net emissions are the primary data source (e.g., AFOLU).

        - Method 1: Uses LUgrs directly and calculates regrowth from it.
          Suitable when gross emissions are observed directly (e.g., harvest data).

        **Cumulative Quantities:**

        Cumulative values track the total historical land use changes and are used
        to calculate land use effects on ecosystem productivity (see
        :class:`EffectLanduseCalculator`).

        **Time Resolution:**

        The interval width is the difference between adjacent time coordinates. For
        example, coordinates spaced by ``0.5 yr`` represent half-year flux intervals,
        and each input rate contributes an amount equal to its rate multiplied by
        ``0.5 yr``. A ``1.4 yr`` regrowth timescale is represented by three intervals
        (``1.5 yr``), with the full regrowth amount spread across them. The final
        coordinate starts the final represented interval.

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
            time_axis=time_axis,
        )

    def calculate_NflxLU_series(
        self,
        NemsLUnet_s: pint.Quantity,
        NemsLUgrs_s: pint.Quantity,
        NemsLUmin_s: pint.Quantity,
        time_axis: Optional[pint.Quantity] = None,
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
        time_axis
            Time coordinates aligned with the input series [yr]. They must be
            one-dimensional, strictly increasing, evenly spaced, and have the same
            length as the flux series. Each coordinate marks the start of the
            corresponding flux interval. Thus, the last value is also applied over
            one interval after the final coordinate. These conditions are checked by
            the experiment config when used through the model. Direct callers must
            supply an axis satisfying these requirements. If omitted, each value
            represents a one-year interval.

        Returns
        -------
        Dictionary containing nitrogen flux time series:

            **Instantaneous Fluxes [GtN/yr]:**

            - NflxLUgrs: Gross land use flux from organic pools
            - NflxLUrgr: Regrowth flux following deforestation
            - NflxLUdfst: Deforestation flux (from organic pools)
            - NflxLUafst: Afforestation flux (to organic pools)
            - NflxLUafst_decay: Afforestation decay flux
            - NflxLUmin: Direct emission from mineral pool

            **Cumulative Quantities [GtN]:**

            - cumsum_NflxLUrgr: Cumulative regrowth nitrogen
            - cumsum_NflxLUdfst: Cumulative deforestation nitrogen
            - cumsum_NflxLUafst: Cumulative afforestation nitrogen
            - cumsum_NflxLUafst_decay: Cumulative afforestation decay

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

        ``NemsLUmin_s`` is a direct mineral-pool loss flux. It is converted to GtN/yr
        and returned as ``NflxLUmin``.

        **Typical Nitrogen Method:**

        Method 1 is typically used for nitrogen because gross deforestation data
        combined with C:N ratios provides more direct estimates of organic nitrogen
        losses than net emission calculations.

        See Also
        --------
        calculate_CflxLU_series: Analogous carbon flux calculation
        :py:meth:`cnit.physics.nitrogen_cycle.NitrogenCycleModel.run`: How these fluxes are used in nitrogen cycle
        """
        fluxes = self._calculate_flux_series(
            net_flux=NemsLUnet_s,
            gross_flux=NemsLUgrs_s,
            method_calc=self.method_nLU,
            unit="GtN/yr",
            prefix="Nflx",
            time_axis=time_axis,
        )

        # This direct mineral-pool loss needs unit conversion, not land-use processing.
        return fluxes | {"NflxLUmin": NemsLUmin_s.to("GtN/yr")}

    def _calculate_flux_series(
        self,
        net_flux: pint.Quantity,
        gross_flux: pint.Quantity,
        method_calc: pint.Quantity,
        unit: str,
        prefix: str,
        time_axis: Optional[pint.Quantity],
    ) -> Dict[str, pint.Quantity]:
        """
        Core flux calculation engine used by both carbon and nitrogen methods.

        This method is shared by the carbon and organic-nitrogen calculations. It
        converts the unit-aware flux inputs to numeric arrays once, performs the
        interval calculations using NumPy, and attaches units to the returned
        quantities. Rates are treated as constant over equal-width intervals.

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
        time_axis
            Time coordinates aligned with each flux value [yr]. Their spacing sets
            the duration represented by each rate. If ``None``, values are treated
            as annual. Otherwise, callers must provide an axis matching the flux
            length that is finite, strictly increasing, and evenly spaced. The model
            experiment config checks these conditions; this internal method assumes
            they are satisfied. Each coordinate starts an interval, including the
            final coordinate.

        Returns
        -------
        Dict[str, pint.Quantity]
            Dictionary with named flux series (see public methods for details).

        ValueError
            If a nonzero input uses a method other than 0 or 1.

        Notes
        -----
        Method 0 derives gross flux sequentially from net flux. At each interval,
        current gross flux is net flux plus regrowth already scheduled for that
        interval. Positive gross flux schedules future regrowth; negative gross flux
        schedules future afforestation decay.

        Method 1 takes gross flux as given and calculates net flux as gross flux minus
        regrowth. For either method, a pulse is applied only within the simulated
        horizon; any part extending beyond the final interval is not included in the
        returned fluxes or cumulative quantities.

        Cumulative outputs are discrete integrals of their corresponding rates:
        each rate is multiplied by the interval duration before summation. When both
        input series are exactly zero, all component outputs are zero.

        **Parameter Relationships:**

        - frc_rgr: Regrowth fraction (varphi)
        - frc_afst_decay: Afforestation decay fraction (phi = 1 - varphi)
        - tau_rgr: Regrowth and afforestation decay timescale in years
        """

        net_flux_a = self._to_numeric_series(net_flux, unit)
        gross_flux_a = self._to_numeric_series(gross_flux, unit)

        if time_axis is None:
            time_step = 1.0
        else:
            time_values = np.asarray(time_axis.to("yr").m, dtype=float)
            time_step = float(np.diff(time_values)[0])

        # Physical ranges for these parameters are checked by the model config.
        # This calculation path only converts them to the units used below.
        frc_rgr = float(self.frc_rgr.to("1").m)
        tau_rgr = float(self.tau_rgr.to("yr").m)
        frc_afst_decay = 1 - frc_rgr

        length_input = net_flux_a.size

        # Handle zero flux case
        if not np.any(net_flux_a) and not np.any(gross_flux_a):
            gross_flux_a = regrowth_flux_a = deforest_flux_a = afforest_flux_a = (
                afforest_decay_flux_a
            ) = np.zeros(length_input, dtype=float)
        else:
            regrowth_flux_a = np.zeros(length_input, dtype=float)
            deforest_flux_a = np.zeros(length_input, dtype=float)
            afforest_flux_a = np.zeros(length_input, dtype=float)
            afforest_decay_flux_a = np.zeros(length_input, dtype=float)

            method_value = method_calc.to("1").m
            if method_value == 0:
                # Method 0: Use net flux as input (e.g., AFOLU emissions)
                self._calculate_from_net_flux(
                    net_flux_a,
                    gross_flux_a,
                    regrowth_flux_a,
                    deforest_flux_a,
                    afforest_flux_a,
                    afforest_decay_flux_a,
                    time_step,
                    frc_rgr,
                    tau_rgr,
                    frc_afst_decay,
                )
            elif method_value == 1:
                # Method 1: Use gross flux as input (e.g., fCharvest emissions)
                # When using gross flux directly, we need to copy the input values first
                self._calculate_from_gross_flux(
                    net_flux_a,
                    gross_flux_a,
                    regrowth_flux_a,
                    deforest_flux_a,
                    afforest_flux_a,
                    afforest_decay_flux_a,
                    time_step,
                    frc_rgr,
                    tau_rgr,
                    frc_afst_decay,
                )
            else:
                raise ValueError(
                    f"Invalid method_calc value: {method_value}. Must be 0 or 1."
                )

        # Integrate each piecewise-constant flux rate over its interval.
        amount_unit = unit.replace("/yr", "")

        return {
            f"{prefix}LUgrs": Q(gross_flux_a, unit),
            f"{prefix}LUrgr": Q(regrowth_flux_a, unit),
            f"{prefix}LUdfst": Q(deforest_flux_a, unit),
            f"{prefix}LUafst": Q(afforest_flux_a, unit),
            f"{prefix}LUafst_decay": Q(afforest_decay_flux_a, unit),
            f"cumsum_{prefix}LUrgr": Q(
                np.cumsum(regrowth_flux_a * time_step), amount_unit
            ),
            f"cumsum_{prefix}LUdfst": Q(
                np.cumsum(deforest_flux_a * time_step), amount_unit
            ),
            f"cumsum_{prefix}LUafst": Q(
                np.cumsum(afforest_flux_a * time_step), amount_unit
            ),
            f"cumsum_{prefix}LUafst_decay": Q(
                np.cumsum(afforest_decay_flux_a * time_step), amount_unit
            ),
        }

    @staticmethod
    def _to_numeric_series(
        values: pint.Quantity,
        unit: str,
    ) -> np.ndarray:
        """Convert a prevalidated flux series to a copied numeric array.

        Experiment configs validate forcing units and series lengths before model
        execution. This helper assumes those preconditions and only converts the
        values to the calculation unit. The copy ensures in-place calculations
        cannot alter caller inputs.

        Parameters
        ----------
        values
            Pint quantity containing one flux value per time interval.
        unit
            Unit to convert values to before dropping their units, such as
            ``"GtC/yr"`` or ``"GtN/yr"``.
        Returns
        -------
        numpy.ndarray
            Copy of the converted flux magnitudes, with units removed.
        """
        return np.asarray(values.to(unit).m, dtype=float).copy()

    def _calculate_from_net_flux(
        self,
        net_flux_values: np.ndarray,
        gross_flux_values: np.ndarray,
        regrowth_flux_values: np.ndarray,
        deforest_flux_values: np.ndarray,
        afforest_flux_values: np.ndarray,
        afforest_decay_flux_values: np.ndarray,
        time_step: float,
        frc_rgr: float,
        tau_rgr: float,
        frc_afst_decay: float,
    ) -> None:
        """Derive gross flux from net input and schedule regrowth and afforestation decay.

        At each interval, ``gross = net + regrowth`` already scheduled for that
        interval. Positive gross flux is deforestation: its interval amount is
        ``gross * time_step``, and the ``frc_rgr`` fraction is spread uniformly over
        ``tau_rgr`` years, rounded to the nearest whole number of model intervals
        (at least one). Negative gross flux is afforestation: its magnitude determines
        the amount whose ``frc_afst_decay`` fraction decays over the same duration.
        The current interval's event starts affecting the following interval.

        All flux arrays are numeric and modified in place. A pulse is truncated at
        the end of the simulation array.
        """
        length = net_flux_values.size
        regrowth_rate_deltas = np.zeros(length + 1, dtype=float)
        decay_rate_deltas = np.zeros(length + 1, dtype=float)
        regrowth_rate = 0.0
        decay_rate = 0.0
        # Each land-use legacy flux (regrowth or afforestation decay) stays active
        # for this many whole model intervals after its triggering event.
        pulse_interval_count = max(1, round(tau_rgr / time_step))
        pulse_duration = pulse_interval_count * time_step

        for i, net_flux in enumerate(net_flux_values):
            # Apply changes scheduled for this interval before reading its fluxes.
            regrowth_rate += regrowth_rate_deltas[i]
            decay_rate += decay_rate_deltas[i]
            regrowth_flux_values[i] = regrowth_rate
            afforest_decay_flux_values[i] = decay_rate

            gross_flux = net_flux + regrowth_flux_values[i]
            gross_flux_values[i] = gross_flux

            if gross_flux > 0:
                deforest_flux_values[i] = gross_flux
                deforestation_amount = gross_flux * time_step
                self._schedule_future_pulse(
                    regrowth_rate_deltas,
                    start_index=i + 1,
                    amount=deforestation_amount * frc_rgr,
                    interval_count=pulse_interval_count,
                    pulse_duration=pulse_duration,
                )
            else:
                afforest_flux_values[i] = -gross_flux
                afforestation_amount = afforest_flux_values[i] * time_step
                self._schedule_future_pulse(
                    decay_rate_deltas,
                    start_index=i + 1,
                    amount=afforestation_amount * frc_afst_decay,
                    interval_count=pulse_interval_count,
                    pulse_duration=pulse_duration,
                )

    def _calculate_from_gross_flux(
        self,
        net_flux_values: np.ndarray,
        gross_flux_values: np.ndarray,
        regrowth_flux_values: np.ndarray,
        deforest_flux_values: np.ndarray,
        afforest_flux_values: np.ndarray,
        afforest_decay_flux_values: np.ndarray,
        time_step: float,
        frc_rgr: float,
        tau_rgr: float,
        frc_afst_decay: float,
    ) -> None:
        """Use gross input to derive net flux and schedule regrowth and afforestation decay.

        Positive gross flux is deforestation; negative gross flux is afforestation.
        Their interval amounts are multiplied by ``time_step`` before applying the
        regrowth or afforestation-decay fraction. Each resulting amount is spread over
        ``tau_rgr`` years, rounded to the nearest whole number of model intervals (at
        least one), starting in the following interval. After all events are scheduled,
        ``net_flux_values`` is updated in place to ``gross_flux_values`` minus
        ``regrowth_flux_values``.

        Inputs and outputs are numeric arrays in compatible rate units. The caller's
        original quantities are not modified because the shared calculation method
        works on copies.
        """
        length = gross_flux_values.size
        regrowth_rate_deltas = np.zeros(length + 1, dtype=float)
        decay_rate_deltas = np.zeros(length + 1, dtype=float)
        regrowth_rate = 0.0
        decay_rate = 0.0
        # Each land-use legacy flux (regrowth or afforestation decay) stays active
        # for this many whole model intervals after its triggering event.
        pulse_interval_count = max(1, round(tau_rgr / time_step))
        pulse_duration = pulse_interval_count * time_step

        for i, gross_flux in enumerate(gross_flux_values):
            # Apply changes scheduled for this interval before recording its fluxes.
            regrowth_rate += regrowth_rate_deltas[i]
            decay_rate += decay_rate_deltas[i]
            regrowth_flux_values[i] = regrowth_rate
            afforest_decay_flux_values[i] = decay_rate

            if gross_flux > 0:
                deforest_flux_values[i] = gross_flux
                deforestation_amount = gross_flux * time_step
                self._schedule_future_pulse(
                    regrowth_rate_deltas,
                    start_index=i + 1,
                    amount=deforestation_amount * frc_rgr,
                    interval_count=pulse_interval_count,
                    pulse_duration=pulse_duration,
                )
            else:
                afforest_flux_values[i] = -gross_flux
                afforestation_amount = afforest_flux_values[i] * time_step
                self._schedule_future_pulse(
                    decay_rate_deltas,
                    start_index=i + 1,
                    amount=afforestation_amount * frc_afst_decay,
                    interval_count=pulse_interval_count,
                    pulse_duration=pulse_duration,
                )

        net_flux_values[:] = gross_flux_values - regrowth_flux_values

    @staticmethod
    def _schedule_future_pulse(
        rate_deltas: np.ndarray,
        start_index: int,
        amount: float,
        interval_count: int,
        pulse_duration: float,
    ) -> None:
        """Schedule a land-use flux pulse with two range-difference updates.

        Here, a *pulse* is the finite-duration flux contribution from one land-use
        event: regrowth after deforestation or afforestation decay. It is applied at
        a constant rate over its scheduled intervals, not as an instantaneous flux.

        Rather than adding a pulse to every affected future interval, this records
        where its rate starts and stops. The running sum of ``rate_deltas`` gives the
        total active rate at each interval. Scheduling one event is O(1), so processing
        ``n`` input intervals takes O(n) time regardless of the pulse duration.

        For example, suppose a pulse contributes 6 units of rate in intervals 1, 2,
        and 3. Record ``+6`` at the start of interval 1 and ``-6`` at the start of
        interval 4. The delta array is ``[0, 6, 0, 0, -6]``; its running sum is
        ``[0, 6, 6, 6, 0]``, the rate to use in each interval. An overlapping pulse
        adds its own start and stop deltas, so the same running sum also combines
        overlapping contributions.

        ``interval_count`` is the configured duration rounded to a whole number of
        intervals, with a minimum of one. It and ``pulse_duration`` are calculated
        once before the time loop because they are the same for every event. If a
        pulse extends past the simulation horizon, its start is recorded but its stop
        is outside the simulated intervals, so its rate remains active through the
        final interval. The simpler :meth:`_add_future_pulse` reference implementation
        adds the same rate directly to each affected interval. The rate is based on
        the full rounded duration even when only part of the pulse is simulated.

        Parameters
        ----------
        rate_deltas
            Length ``n + 1`` array of rate changes at interval boundaries. The extra
            element represents the boundary immediately after the final interval.
        start_index
            Index of the first interval affected by the pulse.
        amount
            Total quantity represented by the pulse.
        interval_count
            Rounded number of model intervals in the pulse.
        pulse_duration
            Rounded pulse duration in years, equal to ``interval_count * time_step``.
        """
        final_boundary = len(rate_deltas) - 1
        if amount == 0 or start_index >= final_boundary:
            return

        rate = amount / pulse_duration
        stop_index = start_index + interval_count

        rate_deltas[start_index] += rate
        if stop_index <= final_boundary:
            rate_deltas[stop_index] -= rate

    @staticmethod
    def _add_future_pulse(
        flux_values: np.ndarray,
        start_index: int,
        amount: float,
        duration: float,
        time_step: float,
    ) -> None:
        """Straightforward interval-by-interval equivalent of ``_schedule_future_pulse``.

        Here, a *pulse* is the finite-duration regrowth or afforestation-decay
        contribution triggered by one land-use event. The contribution is applied at
        a constant rate over multiple intervals, not as an instantaneous flux.

        This method produces the same interval rates as ``_schedule_future_pulse``.
        It is easier to follow because it adds the rate directly to each affected
        interval. Production calculations use ``_schedule_future_pulse``, which
        represents the same interval rates with range differences and a running sum.

        The pulse begins at ``start_index`` and has total amount ``amount`` spread
        uniformly over ``duration`` years. The duration is rounded to the nearest
        whole number of intervals, with a minimum of one interval. A NumPy slice adds
        the resulting constant rate to each affected interval.

        The output array stores interval-average rates. The function modifies it in
        place and truncates the pulse if the available array ends before the pulse
        does. ``time_step`` and ``duration`` must be positive.

        Parameters
        ----------
        flux_values
            Numeric output array of interval-average rates to update.
        start_index
            Index of the first interval affected by the pulse.
        amount
            Total quantity represented by the pulse.
        duration
            Configured duration of the pulse in years.
        time_step
            Width of each model interval in years.
        """
        if amount == 0 or start_index >= len(flux_values):
            return

        interval_count = max(1, round(duration / time_step))
        rounded_duration = interval_count * time_step
        rate = amount / rounded_duration
        full_stop = min(start_index + interval_count, len(flux_values))
        flux_values[start_index:full_stop] += rate
