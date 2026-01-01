"""
Effect Factor Calculators for the CNit Model.

This module defines the calculators responsible for computing the
dimensionless "effect factors" (e.g., from temperature change, CO2, land use)
that modify the rates of core carbon and nitrogen cycle processes.

The effect factors represent multiplicative modifiers applied to baseline process
rates to account for environmental changes and feedbacks. Each calculator class
handles a specific category of effects:

- :class:`EffectLanduseCalculator`: Land use change effects (deforestation, afforestation, regrowth)
- :class:`EffectCO2Calculator`: Atmospheric CO2 fertilization effects
- :class:`EffectTemperatureCalculator`: Temperature-dependent process modifications
- :class:`EffectCarbonNitrogenCouplingCalculator`: Carbon-nitrogen cycle feedback effects

All effect factors are dimensionless multipliers (typically near 1.0) that scale
process rates up or down based on environmental conditions.
"""

from typing import Tuple
from attrs import define, field
import numpy as np
import pint

from ..utils.units import check_units, Q

one = Q(1, "1")


@define
class EffectLanduseCalculator:
    """
    Calculator for land use effects on carbon-nitrogen cycle processes.

    This calculator computes land-use–induced modification factors for carbon and
    nitrogen cycle processes, accounting for deforestation, afforestation, and regrowth
    dynamics. It tracks cumulative land use changes and computes their effects on both
    Net Primary Production (NPP) and organic matter turnover rates.

    The calculator implements two categories of land use effects:

    1. **Effects on NPP**: Adjustments to primary production based on cumulative land
       use changes, blending instantaneous and equilibrium perspectives to capture
       both immediate impacts and long-term recovery dynamics
    2. **Effects on Turnovers**: Short-term immediate effects on litter production,
       decomposition, and soil respiration based on current regrowth flux

    Notes
    -----
    **Box Model Context and the Need for Land Use Effects**

    In a box model with constant turnover times, pool sizes naturally return to
    equilibrium after any perturbation. Without land use effects, the system would
    fully recover to its original state after deforestation or afforestation. The
    land use effects implemented here prevent this unrealistic full return by:

    - Adjusting NPP to reflect permanent changes in productive land area
    - Modifying turnover rates during recovery periods to capture transient dynamics

    **Dual Effect Framework**

    The calculator distinguishes between two complementary mechanisms:

    - **NPP adjustment** (``calculate_eff_LU_NPP``): Modifies primary production
      to reflect net changes in productive land area, preventing full equilibrium
      return. This operates over longer timescales and controls the ultimate
      equilibrium state.

    - **Turnover adjustment** (``calculate_eff_LU_cLPLDSR``, ``calculate_eff_LU_nLPLDSR``):
      Capture immediate, short-term perturbations to organic matter cycling during
      recovery periods. These effects decay as regrowth completes, avoiding
      permanent equilibrium shifts.

    **Regrowth Flux as Proxy for Turnover Effects**

    The turnover effects use regrowth flux as a proxy for land use perturbation
    state because:

    1. Regrowth reflects the recovering state of all past perturbations
    2. It captures lagging effects that extend into the future
    3. It embodies the principle "where there is regrowth, land use
       perturbation still takes effect"
    4. It naturally decays to zero as recovery completes, ensuring no permanent
       change to equilibrium turnover rates

    The regrowth flux provides a single, time-varying indicator of ongoing recovery
    that can be used as a proxy to parameterize how disturbance alters cycling rates
    during the transient recovery period.
    """

    sns_cLP2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"),
        validator=check_units("yr/GtC"),
    )
    """Sensitivity of carbon litter production to regrowth flux [yr/GtC]."""

    sns_cLD2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"),
        validator=check_units("yr/GtC"),
    )
    """Sensitivity of carbon litter decomposition to regrowth flux [yr/GtC]."""

    sns_cSR2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"),
        validator=check_units("yr/GtC"),
    )
    """Sensitivity of carbon soil respiration to regrowth flux [yr/GtC]"""

    sns_nLP2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"),
        validator=check_units("yr/GtC"),
    )
    """Sensitivity of nitrogen litter production to regrowth flux [yr/GtC]"""

    sns_nLD2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"),
        validator=check_units("yr/GtC"),
    )
    """Sensitivity of nitrogen litter decomposition to regrowth flux [yr/GtC]"""

    sns_nSR2LUrgr: pint.Quantity = field(
        default=Q(0, "yr/GtC"),
        validator=check_units("yr/GtC"),
    )
    """Sensitivity of nitrogen soil respiration to regrowth flux [yr/GtC]"""

    def calculate_eff_LU_NPP(
            self,
            cumsum_CflxLUrgr: pint.Quantity,
            cumsum_CflxLUdfst: pint.Quantity,
            cumsum_CflxLUafst: pint.Quantity,
            cumsum_CflxLUafst_decay: pint.Quantity,
            CplsPLS0: pint.Quantity,
            frc_rgr: pint.Quantity,
            frc_afst_decay: pint.Quantity,
    ) -> pint.Quantity:
        """
        Calculate the land use effect on Net Primary Production (NPP).

        This method computes a multiplicative effect factor that accounts for both NPP
        reductions from deforestation and NPP enhancements from afforestation. The
        calculation addresses a fundamental limitation of box models:  with constant
        turnover times, carbon pools always return to their equilibrium state following
        any perturbation.  Land use effects on NPP prevent this unrealistic behavior by
        scaling NPP to reflect persistent changes in productive land area.

        The method uses a refined parameterization that distinguishes between equilibrium
        (long-term potential) and instantaneous (current realized) effects, blending them
        to capture transient dynamics during ecosystem recovery and degradation.

        Parameters
        ----------
        cumsum_CflxLUrgr
            Cumulative regrowth carbon flux [GtC].  Total carbon accumulated through
            natural recovery following past deforestation events.
        cumsum_CflxLUdfst
            Cumulative deforestation carbon flux [GtC]. Total carbon lost through all
            deforestation events.
        cumsum_CflxLUafst
            Cumulative afforestation carbon flux [GtC]. Total carbon gained through
            afforestation and reforestation.
        cumsum_CflxLUafst_decay
            Cumulative afforestation decay carbon flux [GtC]. Total carbon lost from
            afforested areas due to vegetation turnover.
        CplsPLS0
            Initial equilibrium total terrestrial carbon pool size (Plant + Litter + Soil)
            [GtC]. Represents the baseline productive land carbon before land use changes.
        frc_rgr
            Regrowth fraction (varphi) [dimensionless]. Fraction of deforested area that
            can undergo natural recovery (typical values 0.8-0.95).
        frc_afst_decay
            Afforestation decay fraction (phi) [dimensionless]. Fraction of afforested
            carbon that decays due to turnover (phi = 1 - varphi).

        Returns
        -------
        Land use effect factor for NPP [dimensionless].  Multiplicative factor applied
            to baseline NPP:

            - Values < 1: Net NPP reduction (deforestation dominates)
            - Values > 1: Net NPP enhancement (afforestation dominates)
            - Value = 1: No net land use effect

        Notes
        -----
        **The Box Model Equilibrium Return**

        In a box model framework with constant turnover times, carbon pools always return
        to their initial equilibrium state following any perturbation, whether deforestation
        or afforestation. This occurs because:

        ..  math::

            \\frac{dC}{dt} = NPP - \\frac{C}{\\tau}

        At equilibrium (dC/dt = 0): :math:`C_{eq} = NPP \\times \\tau`

        Without land use effects on NPP, any perturbation to C eventually returns to the
        original equilibrium regardless of land use history. This is unrealistic because
        deforestation should reduce productive land area (lower equilibrium), while
        afforestation should increase it.

        Land use effects on NPP solve this by adjusting the equilibrium target to reflect
        net changes in productive land area, ensuring that long-term pool sizes reflect
        the persistent impact of land use change.

        **Deforestation Effect: Balancing Two Perspectives**

        Two complementary approaches characterize the deforestation effect:

        **1. Equilibrium Perspective** (:math:`\\epsilon_{LUdfst}^{eqm}`)

        Represents the ultimate long-term NPP reduction assuming all potential regrowth
        is complete:

        .. math::

            \\epsilon_{LUdfst}^{eqm} = \\frac{(1-\\varphi) \\sum_{0}^{t} LUdfst}{C_{Land0}}

        where :math:`\\varphi` is the regrowth fraction (`frc_rgr`) and :math:`C_{Land0}`
        is the initial equilibrium land carbon (`CplsPLS0`).

        **Physical interpretation**: Only the non-regrowing fraction : math:`(1-\\varphi)`
        represents permanent loss of productive land.  If :math:`\\varphi = 0. 9`, then 90%
        of deforested area eventually recovers, so only 10% contributes to long-term NPP
        reduction.

        **Limitation**: This assumes NPP immediately adjusts to the new equilibrium level
        after deforestation, potentially underestimating the effect during the regrowth
        period.  The new equilibrium corresponds to the end of the regrowth period
        (:math:`\\tau_{rgr}`), so at time t during active regrowth, :math:`\\epsilon_{LUdfst}^{eqm}`
        underestimates the realized impact because it doesn't account for lagged regrowth.

        **2. Instantaneous Perspective** (:math:`\\epsilon_{LUdfst}^{inst}`)

        Represents the realized net impact at the current time:

        ..  math::

            \\epsilon_{LUdfst}^{inst} = \\frac{\\sum_{0}^{t}(LUdfst - LUrgr)}{C_{Land0}}

        **Physical interpretation**: The actual net carbon loss to date. This captures the
        current imbalance between cumulative deforestation and cumulative regrowth.

        **Temporal behavior**:

        - Early in regrowth period: :math:`\\epsilon_{LUdfst}^{inst}` is larger than
          :math:`\\epsilon_{LUdfst}^{eqm}` because regrowth hasn't caught up
        - Long-term (after :math:`\\tau_{rgr}`): Both converge to the same value

        **Limitations of Using Either Perspective Alone**

        **Why** :math:`\\epsilon_{LUdfst}^{inst}` **alone overestimates the effect:**

        1. **Simplified regrowth dynamics**: The model assumes constant regrowth flux over
           :math:`\\tau_{rgr}`, which doesn't capture rapid initial regrowth in real
           ecosystems. This underestimates cumulative :math:`LUrgr`,
           exaggerating the net loss :math:`(LUdfst - LUrgr)`.

        2. **Missing high NPP of young vegetation**: Young regrowing forests have high NPP
           potential relative to their carbon stock, which isn't captured in the simple
           cumulative carbon balance.  The instantaneous effect treats all carbon equally,
           missing the fact that a small amount of regrowth carbon represents substantial
           NPP recovery.

        **Why** :math:`\\epsilon_{LUdfst}^{eqm}` **alone underestimates the effect:**

        - Assumes instantaneous adjustment to the post-regrowth equilibrium, ignoring the
          actual lag in recovery.  During active regrowth, the realized NPP reduction is
          larger than the long-term equilibrium would suggest.

        **Combined Approach:  Weighted Blending**

        The effective deforestation effect blends both perspectives with an adaptive
        weighting scheme:

        .. math::

            \\epsilon_{LUdfst} = \\frac{1}{2}\\left[(1+0.8\\varphi) \\times \\epsilon_{LUdfst}^{eqm}
            + (1-0.8\\varphi) \\times \\epsilon_{LUdfst}^{inst}\\right]

        **Rationale for the 0.8φ tuning factor:**

        The factor : math:`0.8\\varphi` prevents problematic behavior at extreme regrowth
        fractions:

        - **Problem at high φ** (e.g., :math:`\\varphi = 1`):

          * :math:`\\epsilon_{LUdfst}^{eqm} = 0` (all area recovers, no long-term effect)
          * :math:`\\epsilon_{LUdfst}^{inst}` remains large during regrowth
          * Simple average would give :math:`\\epsilon_{LUdfst} = 0.5 \\times \\epsilon_{LUdfst}^{inst}`,
            overemphasizing the exaggerated instantaneous response

        - **Solution with 0.8φ weighting**:

          * At :math:`\\varphi = 1`: weights are (0.9, 0.1) rather than (0.5, 0.5)
          * Minimum 10% contribution from instantaneous effect: :math:`(1 - 0.8 \\times 1) / 2 = 0.1`
          * Maximum 90% from equilibrium:  :math:`(1 + 0.8 \\times 1) / 2 = 0.9`
          * This preserves physically meaningful balance across the full regrowth range
            (0 ≤ φ ≤ 1)

        - **Behavior across regrowth fractions**:

          * Low φ (more permanent deforestation): Both perspectives matter similarly
          * High φ (more regrowth): Equilibrium perspective dominates, appropriately
            reducing reliance on the exaggerated instantaneous effect

        **Afforestation Effect: Similar Blending Without Tuning**

        The equilibrium and instantaneous afforestation effects mirror the deforestation
        formulation:

        .. math::

            \\epsilon_{LUafst}^{eqm} = \\frac{(1-\\phi) \\sum_{0}^{t} LUafst}{C_{Land,0}}

        .. math::

            \\epsilon_{LUafst}^{inst} = \\frac{\\sum_{0}^{t}(LUafst - LUafst\\_decay)}{C_{Land,0}}

        where :math:`\\phi` is the afforestation decay fraction (`frc_afst_decay`).

        The effective afforestation effect uses a simple unweighted average:

        .. math::

            \\epsilon_{LUafst} = \\frac{1}{2}\\left[\\epsilon_{LUafst}^{eqm} + \\epsilon_{LUafst}^{inst}\\right]

        **Why no φ adjustment for afforestation:**

        The decay of carbon in mature forests is a gradual process, unlike the potentially
        rapid initial regrowth after deforestation. The instantaneous effect for afforestation
        doesn't suffer from the same exaggeration issues:

        - Afforestation decay is slow and approximately linear (constant flux assumption is
          more realistic)
        - There's no equivalent to the "high NPP of young vegetation" bias working in the
          opposite direction
        - Therefore, equilibrium and instantaneous perspectives are equally valid and don't
          require adaptive weighting

        **Total Land Use Effect on NPP**

        The final multiplicative effect combines deforestation (NPP reduction) and
        afforestation (NPP amplification):

        .. math::

            \\epsilon_{LU(NPP)} = 1 - \\epsilon_{LUdfst} + \\epsilon_{LUafst}

        **Physical interpretation:**

        - :math:`-\\epsilon_{LUdfst}`: Reduces NPP proportional to net deforestation impact
        - :math:`+\\epsilon_{LUafst}`: Increases NPP proportional to net afforestation impact
        - Baseline of 1 represents no land use change

        **Example scenarios:**

        1. **Heavy deforestation, no afforestation**:  :math:`\\epsilon_{LUdfst} = 0.2`,
           :math:`\\epsilon_{LUafst} = 0` → :math:`\\epsilon_{LU(NPP)} = 0.8` (20% NPP reduction)

        2. **Heavy afforestation, no deforestation**: :math:`\\epsilon_{LUdfst} = 0`,
           :math:`\\epsilon_{LUafst} = 0.15` → :math:`\\epsilon_{LU(NPP)} = 1.15` (15% NPP increase)

        3. **Mixed**:  :math:`\\epsilon_{LUdfst} = 0.1`, :math:`\\epsilon_{LUafst} = 0.12`
           → :math:`\\epsilon_{LU(NPP)} = 1.02` (2% net NPP increase)

        See Also
        --------
        :py:meth:`cnit.physics.landuse.LanduseCalculator.calculate_CflxLU_series`: Calculates cumulative land use fluxes
        :py:meth:`cnit.physics.carbon_cycle.CarbonNPPLPRCalculator.calculate`: Uses this effect to modify NPP
        """

        # Calculate non-regrowing deforestation carbon for equilibrium perspective
        # This is the permanent carbon loss: only (1 - frc_rgr) fraction won't regrow
        cumsum_CflxLUdfst_nonrgr_eqm = (1 - frc_rgr) * cumsum_CflxLUdfst

        # Calculate non-regrowing deforestation carbon for instantaneous perspective
        # This is the actual net carbon loss to date: deforestation minus regrowth
        cumsum_CflxLUdfst_nonrgr_instant = cumsum_CflxLUdfst - cumsum_CflxLUrgr

        # Apply weighted blending with 0.8*frc_rgr tuning factor
        # This balances equilibrium and instantaneous perspectives, preventing
        # over-reliance on either at extreme regrowth fractions
        weight_frc_rgr = 0.8
        weight_rgr_eqm = weight_frc_rgr * frc_rgr
        eff_dfst_effective = (
                0.5
                * (
                        cumsum_CflxLUdfst_nonrgr_eqm * (1 + weight_rgr_eqm)
                        + cumsum_CflxLUdfst_nonrgr_instant * (1 - weight_rgr_eqm)
                )
                / CplsPLS0
        )

        # Calculate persistent afforestation carbon for equilibrium perspective
        # Only (1 - frc_afst_decay) fraction will persist long-term
        cumsum_CflxLUafst_persist_eqm = (1 - frc_afst_decay) * cumsum_CflxLUafst

        # Calculate persistent afforestation carbon for instantaneous perspective
        # Actual net carbon gain to date: afforestation minus decay
        cumsum_CflxLUafst_persist_instant = cumsum_CflxLUafst - cumsum_CflxLUafst_decay

        # Simple average for afforestation (no phi adjustment needed)
        # Mature forest decay is gradual, so no exaggeration issues
        eff_afst_effective = (
                0.5
                * (cumsum_CflxLUafst_persist_eqm + cumsum_CflxLUafst_persist_instant)
                / CplsPLS0
        )

        # Total effect: 1 (baseline) - reduction (dfst) + amplification (afst)
        return 1 - eff_dfst_effective + eff_afst_effective

    def calculate_eff_LU_cLPLDSR(
            self,
            CflxLUrgr: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity, pint.Quantity]:
        """
        Calculate land use effect factors on carbon cycle turnover processes.

        Computes exponential response factors for carbon litter production (cLP),
        litter decomposition (cLD), and soil respiration (cSR). These factors
        represent immediate, short-term effects of ongoing land use perturbations
        on organic matter cycling rates during recovery periods.

        Parameters
        ----------
        CflxLUrgr : pint.Quantity
            Current regrowth carbon flux [GtC/yr]. Serves as a proxy for the
            state of land use perturbation and ongoing recovery.

        Returns
        -------
        eff_LU_cLP : pint.Quantity
            Land use effect on carbon litter production [dimensionless].
        eff_LU_cLD : pint.Quantity
            Land use effect on carbon litter decomposition [dimensionless].
        eff_LU_cSR : pint.Quantity
            Land use effect on carbon soil respiration [dimensionless].

        Notes
        -----
        **Why Turnover Effects Are Necessary**

        The NPP adjustment (``calculate_eff_LU_NPP``) prevents full equilibrium
        return in the long term, but it primarily controls the ultimate equilibrium
        state. These direct turnover effects are needed to capture the shorter-term
        immediate effects of land use perturbation during active recovery periods.

        Without turnover effects, transient dynamics during recovery would be
        missing from the model, even though the long-term equilibrium is correctly
        adjusted.

        **Mathematical Formulation**

        Each effect factor follows an exponential response function:

        .. math::

            \\epsilon_{LU(i)} = \\exp(s_{i2LUrgr} \\times LUrgr)

        where :math:`s_{i2LUrgr}` is the sensitivity parameter for process
        :math:`i` (litter production, litter decomposition, or soil respiration)
        and :math:`LUrgr` is the current regrowth flux.

        **Why Regrowth Flux Is a Good Proxy**

        Using regrowth flux to parameterize turnover effects is appropriate because:

        1. **Reflects recovering state**: Regrowth magnitude indicates how much
           past perturbation is still influencing the system
        2. **Captures lagging effects**: Active regrowth means perturbation effects
           extend into the future (ongoing recovery)
        3. **Embodies "where regrowth occurs, perturbation persists"**: The spatial
           and temporal extent of regrowth maps to perturbation influence
        4. **Natural decay**: As regrowth completes, :math:`LUrgr \\to 0` and
           effects vanish, preventing permanent equilibrium shifts
        5. **Single time-varying indicator**: Provides a convenient scalar measure
           of system-wide recovery status

        **Interpretation of Sensitivities**

        - **Positive** :math:`s_{i2LUrgr}`: Process rate increases with regrowth
          (e.g., high litter production in rapidly accumulating young forests)
        - **Negative** :math:`s_{i2LUrgr}`: Process rate decreases with regrowth
          (e.g., reduced soil respiration as system recovers from disturbance pulse)
        - **Zero** :math:`s_{i2LUrgr}`: Process rate unaffected by land use state

        **No Long-Term Equilibrium Change**

        Because regrowth is not forever (recovery eventually completes), these
        turnover effects naturally decay to 1.0 as :math:`LUrgr \\to 0`. This means:

        - Short-term perturbations to cycling rates are captured
        - Long-term equilibrium turnover rates remain unchanged
        - Only the NPP effect permanently alters the equilibrium state

        """
        # Calculate effects using exponential response functions

        return (
            np.exp(self.sns_cLP2LUrgr * CflxLUrgr),
            np.exp(self.sns_cLD2LUrgr * CflxLUrgr),
            np.exp(self.sns_cSR2LUrgr * CflxLUrgr),
        )

    def calculate_eff_LU_nLPLDSR(
            self,
            CflxLUrgr: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity, pint.Quantity]:
        """
        Calculate land use effect factors on nitrogen turnover processes.

        Computes exponential response factors for nitrogen litter production (nLP),
        litter decomposition (nLD), and soil respiration (nSR). These factors modify
        nitrogen turnover process rates based on the current regrowth flux.

        The effect factors follow exponential response functions analogous to the
        carbon cycle effects (see :meth:`calculate_eff_LU_cLPLDSR`).

        Parameters
        ----------
        CflxLUrgr : pint.Quantity
            Current regrowth carbon flux [GtC/yr]. Used as a proxy for ecosystem
            recovery state affecting nitrogen cycling.

        Returns
        -------
        eff_LU_nLP : pint.Quantity
            Land use effect on nitrogen litter production [dimensionless].
        eff_LU_nLD : pint.Quantity
            Land use effect on nitrogen litter decomposition [dimensionless].
        eff_LU_nSR : pint.Quantity
            Land use effect on nitrogen soil respiration (mineralization) [dimensionless].

        Notes
        -----
        While parameterized using carbon regrowth flux, these factors affect nitrogen
        cycling rates. This reflects the coupled nature of C-N cycling in ecosystems,
        where recovery of carbon stocks is accompanied by changes in nitrogen turnover.

        Typical patterns include enhanced nitrogen mineralization during early regrowth
        (positive sensitivity for nSR) as accumulated organic matter decomposes.
        """
        # Calculate effects using exponential response functions
        return (
            np.exp(self.sns_nLP2LUrgr * CflxLUrgr),
            np.exp(self.sns_nLD2LUrgr * CflxLUrgr),
            np.exp(self.sns_nSR2LUrgr * CflxLUrgr),
        )


@define
class EffectCO2Calculator:
    """
    Calculator of CO2 fertilization effects on Net Primary Production (NPP).

    Implements multiple formulations of the CO2 fertilization effect, including
    logarithmic, rectangular hyperbolic (Gifford), and sigmoid forms. The calculator
    can blend between these formulations based on the ``method_CO2_NPP`` parameter
    to represent different assumptions about CO2 response saturation behavior.

    The CO2 effect represents how elevated atmospheric CO2 concentrations enhance
    photosynthetic rates and water use efficiency in vegetation, leading to increased
    NPP. Different formulations capture different assumptions about saturation behavior
    at high CO2 levels and the shape of the CO2 response curve.

    Notes
    -----
    **Physical Basis of CO2 Fertilization**

    Elevated CO2 enhances plant productivity through two main mechanisms:

    1. **Direct photosynthetic enhancement**: Higher CO2 increases the rate of
       carbon fixation in the Calvin cycle (substrate availability)
    2. **Improved water use efficiency**: Plants can partially close stomata while
       maintaining CO2 uptake, reducing water loss

    The strength of CO2 fertilization varies by:

    - Plant functional type (C3 vs. C4 photosynthesis)
    - Nutrient availability (progressive nitrogen limitation)
    - Water availability
    - Temperature regime

    **Choice of Formulation**

    Different formulations represent different hypotheses about CO2 response:

    - **Logarithmic**: No saturation, effect continues growing indefinitely.
      Simple and commonly used, but may overestimate effects at very high CO2.

    - **Rectangular hyperbolic**: Moderate saturation following Michaelis-Menten
      kinetics. Based on photosynthesis theory (Gifford 1980).

    - **Sigmoid**: Strong saturation with explicit maximum effect. Most conservative
      at high CO2, captures potential acclimation or nutrient limitation feedbacks.

    **Calibration Strategy**

    The rectangular hyperbolic sensitivity :math:`s^{rect}_{CO_2}` is automatically
    derived to match the logarithmic form at two reference concentrations (340 and
    680 ppm), ensuring smooth transitions when blending formulations.
    """

    sns_CO2_log: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """
    Sensitivity of NPP to atmospheric |CO2| concentration changes (logarithmic form) [dimensionless].
    """

    sns_CO2_sig: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """
    Sensitivity of NPP to atmospheric |CO2| concentration changes (sigmoid form) [ppm].
    """

    eff_CO2_sig_max: pint.Quantity = field(
        default=Q(2, "1"),
        validator=check_units("1"),
    )
    """
    Maximum CO2 effect factor for the sigmoid form [dimensionless].
    """

    CO2ref: pint.Quantity = field(
        default=Q(284.875, "ppm"),
        validator=check_units("ppm"),
    )
    """
    Atmospheric |CO2| concentration at pre-industrial level [ppm].
    """

    CO2b: pint.Quantity = field(
        default=Q(31, "ppm"),
        validator=check_units("ppm"),
    )
    """
    Atmospheric |CO2| concentration when NPP = 0 (for rectangular hyperbolic effect, 
    default 31 ppm) [ppm].
    """

    method_CO2_NPP: pint.Quantity = field(
        default=Q(0, "1"),
        validator=check_units("1"),
    )
    """
    |CO2| effect factor calculation method (value from 0 to 2), 
    otherwise eff_CO2 is set equal to 1.
    """

    def calculate_eff_CO2_NPP(
            self,
            CO2: pint.Quantity,
    ) -> pint.Quantity:
        """
        Calculate the CO2 fertilization effect factor for Net Primary Production (NPP).

        Computes a multiplicative factor representing how elevated CO2 enhances NPP
        relative to the reference concentration. The calculation can use pure or
        blended formulations depending on the ``method_CO2_NPP`` parameter.

        Parameters
        ----------
        CO2 : pint.Quantity
            Current atmospheric CO2 concentration [ppm].

        Returns
        -------
        pint.Quantity
            CO2 effect factor on NPP [dimensionless].

            - Values > 1: CO2 fertilization (enhanced NPP)
            - Values = 1: No effect (CO2 at reference level)
            - Values < 1: Suppression (rare, only at very low CO2 below reference)

        Notes
        -----
        **Three CO2 Effect Formulations**

        1. **Logarithmic Form** (simple, commonly used):

           .. math::

              \\epsilon^{log}_{CO_2} = 1 + s^{log}_{CO_2} \\times \\ln(CO_2/CO_{2,ref})

           where :math:`CO_2` is the current atmospheric CO2 concentration and
           :math:`CO_{2,ref}` is the reference concentration (typically pre-industrial).

           **Properties**:

           - Effect factor = 1 at :math:`CO_2 = CO_{2,ref}`
           - Logarithmic relationship implies diminishing but never-ceasing returns
           - No explicit saturation at high CO2
           - Commonly used in Earth system models (e.g., OSCAR, simple IAMs)

           **Interpretation of** :math:`s^{log}_{CO_2}` **(β parameter)**:

           - :math:`s^{log}_{CO_2} = 0.4`: Doubling CO2 increases NPP by ~28%
             (:math:`1 + 0.4 \\times \\ln(2) \\approx 1.277`)
           - :math:`s^{log}_{CO_2} = 0.6`: Doubling CO2 increases NPP by ~42%
           - Literature estimates typically range from 0.2 to 0.6

        2. **Rectangular Hyperbolic Form** (Gifford formulation, includes saturation):

           .. math::
              \\epsilon^{rect}_{CO_2} = \\frac{1/(CO_{2,ref}-CO_{2,b}) + s^{rect}_{CO_2}}{
              1/(CO_2-CO_{2,b}) + s^{rect}_{CO_2}}

           where :math:`CO_{2,b}` is the CO2 compensation point (NPP = 0, typically ~31 ppm).

           **Derivation of** :math:`s^{rect}_{CO_2}`:

           To ensure smooth blending with the logarithmic form, :math:`s^{rect}_{CO_2}`
           is derived by matching the effect ratio at two reference concentrations
           (340 and 680 ppm):

           .. math::
              r = \\frac{\\epsilon^{log}_{CO_2}(680)}{\\epsilon^{log}_{CO_2}(340)} =
              \\frac{\\epsilon^{rect}_{CO_2}(680)}{\\epsilon^{rect}_{CO_2}(340)}

           .. math::
              r = \\frac{1 + s^{log}_{CO_2} \\times \\ln(680/CO_{2,ref})}{
              1 + s^{log}_{CO_2} \\times \\ln(340/CO_{2,ref})}

           .. math::
              s^{rect}_{CO_2} = \\frac{(680-CO_{2,b}) - r(340-CO_{2,b})}{
              (r-1)(680-CO_{2,b})(340-CO_{2,b})}

           **Properties**:

           - Based on Michaelis-Menten enzyme kinetics
           - Exhibits saturation at high CO2 levels
           - Approaches asymptote as :math:`CO_2 \\to \\infty`
           - More physiologically grounded than logarithmic form
           - Used in some process-based terrestrial biosphere models

        3. **Sigmoid Form** (flexible, allows specified maximum effect):

           .. math::
              \\epsilon^{sig}_{CO_2} = \\frac{\\epsilon^{sig}_{CO_2,max}}{
              1 + (\\epsilon^{sig}_{CO_2,max}-1) \\times
              \\exp(-s^{sig}_{CO_2}(CO_2/CO_{2,ref} - 1))}

           **Properties**:

           - Transitions smoothly from 1 at :math:`CO_{2,ref}` to
             :math:`\\epsilon^{sig}_{CO_2,max}` at high CO2
           - Steepness controlled by :math:`s^{sig}_{CO_2}`
           - Allows explicit specification of maximum possible enhancement
           - Can represent nutrient limitation or acclimation feedbacks
           - Most conservative at very high CO2 concentrations

        **Blending Between Formulations**

        The final effect factor is a linear blend controlled by ``method_CO2_NPP``:

        .. math::

           \\epsilon_{CO_2} = \\begin{cases}
              (1-m) \\times \\epsilon^{log}_{CO_2} + m \\times \\epsilon^{rect}_{CO_2}
              & 0 \\le m \\le 1 \\\\
              (2-m) \\times \\epsilon^{rect}_{CO_2} + (m-1) \\times \\epsilon^{sig}_{CO_2}
              & 1 < m \\le 2 \\\\
              1 & \\text{otherwise}
           \\end{cases}

        where :math:`m = \\mathrm{method\\_CO2\\_NPP}`.

        **Blending rationale**:

        - Allows exploration of uncertainty in CO2 response shape
        - Smooth transitions between formulations (no discontinuities)
        - Can be calibrated against observations or expert judgment
        - Value outside [0,2] disables CO2 fertilization entirely
        """
        # calculate the logarithmic form CO2 effect factors
        eff_CO2_log = 1 + self.sns_CO2_log * np.log(CO2 / self.CO2ref)

        # calculate the rectangular hyperbolic form (GIFFORD) CO2 effect factors
        CO2_680 = Q(680, "ppm")
        CO2_340 = Q(340, "ppm")
        r = (1 + self.sns_CO2_log * np.log(CO2_680 / self.CO2ref)) / (
                1 + self.sns_CO2_log * np.log(CO2_340 / self.CO2ref)
        )
        s_CO2_rect = ((CO2_680 - self.CO2b) - r * (CO2_340 - self.CO2b)) / (
                (r - 1) * (CO2_680 - self.CO2b) * (CO2_340 - self.CO2b)
        )
        eff_CO2_rect = (1 / (self.CO2ref - self.CO2b) + s_CO2_rect) / (
                1 / (CO2 - self.CO2b) + s_CO2_rect
        )

        # calculate the sigmoid form (ALEX) CO2 effect factors
        eff_CO2_sig = self.eff_CO2_sig_max / (
                1
                + (self.eff_CO2_sig_max - 1)
                * np.exp(-self.sns_CO2_sig * (CO2 / self.CO2ref - 1))
        )

        # calculate the effective effect factor (linear combination of different
        # forms, default value = 1)
        if 0 <= self.method_CO2_NPP <= 1:
            eff_CO2_NPP = (
                                  1 - self.method_CO2_NPP
                          ) * eff_CO2_log + self.method_CO2_NPP * eff_CO2_rect
        elif 1 < self.method_CO2_NPP <= 2:
            eff_CO2_NPP = (2 - self.method_CO2_NPP) * eff_CO2_rect + (
                    self.method_CO2_NPP - 1
            ) * eff_CO2_sig
        else:
            eff_CO2_NPP = 1

        return eff_CO2_NPP


@define
class EffectTemperatureCalculator:
    """
    Calculator of temperature effects for ecological processes.

    Computes exponential and sigmoid temperature response factors for carbon and
    nitrogen cycle processes. Temperature effects generally accelerate biological
    rates (photosynthesis, respiration, decomposition, nitrogen transformations)
    following temperature-dependent kinetics described by Arrhenius-type or Q10
    relationships.

    All temperature effects are computed relative to a reference state (typically
    pre-industrial), using temperature change (dT) rather than absolute temperature.
    This approach simplifies parameter interpretation and maintains consistency
    with transient climate change simulations in Earth system models.

    Notes
    -----
    **Exponential Temperature Response (Q10 Framework)**

    Most temperature effects use exponential formulations:

    .. math::
       \\epsilon_T = \\exp(s_T \\times \\Delta T)

    where :math:`s_T` is the sensitivity parameter and :math:`\\Delta T` is
    temperature change from reference. This corresponds to the Q10 relationship:

    .. math::
       Q_{10} = \\exp(10 \\times s_T)

    Common Q10 values and corresponding sensitivities:

    - Q10 = 1.5 → s = 0.0405 K⁻¹ (weak temperature dependence)
    - Q10 = 2.0 → s = 0.0693 K⁻¹ (moderate, commonly used for soil respiration)
    - Q10 = 2.5 → s = 0.0916 K⁻¹ (strong temperature dependence)
    - Q10 = 3.0 → s = 0.1099 K⁻¹ (very strong, e.g., denitrification)

    **NPP Temperature Response: Monotonic vs. Peaked**

    The NPP temperature effect can represent two different hypotheses:

    1. **Monotonic increase** (``method_dT_NPP = 0``): Simple exponential response,
       suitable for temperature-limited ecosystems (e.g., boreal, arctic)

    2. **Peaked response** (``method_dT_NPP = 1``): Sigmoid with asymptote at 2×,
       representing thermal optima with reduced productivity at high temperatures
       (heat stress, drought coupling). Suitable for temperate to tropical systems.
    """

    sns_NPP2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of NPP to temperature changes, exponential form [K\ :sup:`-1`]
    """

    sns_NPP2dT_sig: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of NPP to temperature changes, sigmoid form [K\ :sup:`-1`]
    """

    method_dT_NPP: pint.Quantity = field(
        default=Q(0, "1"),
        validator=check_units("1"),
    )
    """
        Blending parameter between exponential (0) and sigmoid (1) NPP responses
        [dimensionless, 0-1] [dimensionless]
    """

    sns_LPR2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of litter production respiration to temperature changes [K\ :sup:`-1`]
    """

    sns_cLP2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of carbon litter production to temperature changes [K\ :sup:`-1`]
    """

    sns_cLD2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of carbon litter decomposition to temperature changes [K\ :sup:`-1`]
    """

    sns_cSR2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of carbon soil respiration (soil organic matter decomposition) to 
        temperature changes [K\ :sup:`-1`]
    """

    sns_PU2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of nitrogen plant uptake to temperature changes [K\ :sup:`-1`]
    """

    sns_BNF2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of biological nitrogen fixation to temperature changes [K\ :sup:`-1`]
    """

    sns_nLP2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of nitrogen litter production to temperature changes [K\ :sup:`-1`]
    """

    sns_nLD2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of nitrogen litter decomposition to temperature changes [K\ :sup:`-1`]
    """

    sns_nSR2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of nitrogen soil respiration (soil organic matter decomposition) to
        temperature changes [K\ :sup:`-1`]
    """

    sns_nLSgas2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """
        Sensitivity of mineral nitrogen gaseous loss to temperature changes [K\
        :sup:`-1`]
    """

    def calculate_eff_dT_NPPLPR(
            self,
            dT: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity]:
        """
        Calculate temperature effects for NPP and litter production respiration (LPR).

        Computes multiplicative effect factors representing how temperature changes
        modify NPP and autotrophic respiration. The NPP effect uses a blended
        exponential-sigmoid formulation to optionally capture both temperature
        limitation and high-temperature stress.

        Parameters
        ----------
        dT
            Change in temperature from pre-industrial [K]

        Returns
        -------
        eff_dT_NPP : pint.Quantity
            Temperature effect on Net Primary Production [dimensionless].
            Values > 1 indicate warming enhances NPP; values < 1 indicate suppression.
        eff_dT_LPR : pint.Quantity
            Temperature effect on litter production respiration [dimensionless].
            Typically > 1 for warming as autotrophic respiration increases.

        Notes
        -----
        **NPP Temperature Effect (Blended Exponential-Sigmoid)**

            .. math::

                \\epsilon_{T(NPP)} = (1 - m) \\times e^{s_{NPP2dT}^{exp} \\times
                \\Delta T} +
                m \\times \\frac{2}{1 + e^{-s_{NPP2dT}^{sig} \\times \\Delta T}}

            where:

            - :math:`m = \\mathrm{method\\_dT\\_NPP}` (blending parameter, 0-1)
            - :math:`s_{NPP2dT}^{exp} = \\mathrm{sns\\_NPP2dT}` (exponential sensitivity)
            - :math:`s_{NPP2dT}^{sig} = \\mathrm{sns\\_NPP2dT\\_sig}` (sigmoid sensitivity)
            - :math:`\\Delta T = dT` (temperature change)

        **Interpretation by** :math:`m` **value**:

        - :math:`m = 0`: Pure exponential, monotonic increase with temperature.
          Suitable for temperature-limited ecosystems (boreal, arctic) where
          warming consistently enhances productivity.

        - :math:`m = 1`: Pure sigmoid, approaches asymptote of 2× at high temperatures.
          Represents systems with thermal optima where extreme warming causes
          heat stress, drought, or other limitations. Suitable for temperate to
          tropical ecosystems.

        - :math:`0 < m < 1`: Weighted blend, allows intermediate behavior.

        **LPR Temperature Effect (Simple Exponential)**

        .. math::

            \\epsilon_{T(LPR)} = e^{s_{LPR2dT} \\times \\Delta T}

        where :math:`s_{LPR2dT} = \\mathrm{sns\\_LPR2dT}`.

        Autotrophic respiration (growth and maintenance respiration associated
        with litter production) typically follows simple exponential temperature
        dependence with Q10 ~ 2.0-2.5 (:math:`s_{LPR}` ~ 0.07-0.09 K⁻¹).
        """
        # Calculate NPP temperature effect using a combination of exponential and sigmoid models
        exponential_component = np.exp(dT * self.sns_NPP2dT)
        sigmoid_component = 2 / (1 + np.exp(-self.sns_NPP2dT_sig * dT))
        eff_dT_NPP = (
                             1 - self.method_dT_NPP
                     ) * exponential_component + self.method_dT_NPP * sigmoid_component

        return eff_dT_NPP, np.exp(dT * self.sns_LPR2dT)

    def calculate_eff_dT_cLPLDSR(
            self,
            dT: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity, pint.Quantity]:
        """
        Calculate temperature effects for carbon cycle turnover processes.

        Computes exponential temperature response factors for carbon litter
        production (LP), litter decomposition (LD), and soil respiration (SR).
        These factors modify organic matter turnover rates in response to
        temperature changes.

        Parameters
        ----------
        dT : pint.Quantity
            Change in temperature from reference [K].

        Returns
        -------
        eff_dT_cLP : pint.Quantity
            Temperature effect on carbon litter production [dimensionless].
        eff_dT_cLD : pint.Quantity
            Temperature effect on carbon litter decomposition [dimensionless].
        eff_dT_cSR : pint.Quantity
            Temperature effect on carbon soil respiration (SOM decomposition)
            [dimensionless].

        Notes
        -----
        **Mathematical Formulation**

        Each effect follows a simple exponential (Q10-type) response:

        .. math::
           \\epsilon_{T(i)} = \\exp(s_{i2dT} \\times \\Delta T)

        **Soil Respiration: Critical for Climate-Carbon Feedback**

        The soil respiration temperature sensitivity (``sns_cSR2dT``) is a key
        parameter determining climate-carbon feedback strength:

        .. math::
           \\epsilon_{T(cSR)} = \\exp(s_{cSR2dT} \\times \\Delta T)

        **Typical Q10 values for soil respiration**:

        - Conservative: Q10 = 1.5 (:math:`s_{cSR2dT}` = 0.041 K⁻¹)
        - Moderate: Q10 = 2.0 (:math:`s_{cSR2dT}` = 0.069 K⁻¹)
        - High: Q10 = 2.5 (:math:`s_{cSR2dT}` = 0.092 K⁻¹)
        - Very high: Q10 = 3.0 (:math:`s_{cSR2dT}` = 0.110 K⁻¹)

        **Feedback loop**: Warming → Enhanced soil respiration → CO2 release →
        More warming. Magnitude depends critically on Q10.

        **Example**: For 3 K warming:

        - Q10 = 1.5 → 16% increase in soil respiration
        - Q10 = 2.0 → 23% increase
        - Q10 = 3.0 → 37% increase

        This difference propagates to multi-decadal cumulative carbon losses.

        **Litter Production and Decomposition**

        - **Litter production** (``sns_cLP2dT``): May be positive (stress-induced
          turnover) or near zero (controlled by other factors like NPP and
          allocation). Often less temperature-sensitive than decomposition.

        - **Litter decomposition** (``sns_cLD2dT``): Typically positive with
          Q10 ~ 2.0-2.5, similar to soil respiration. Microbial decomposer
          activity accelerates with warmth.

        **Limitations and Caveats**

        1. **Moisture coupling**: Temperature effects on decomposition are often
           modulated by moisture availability (not explicitly represented). In
           dry conditions, warming may reduce decomposition.

        2. **Substrate quality**: Temperature sensitivity may vary with litter
           and soil organic matter chemistry (labile vs. recalcitrant).

        3. **Acclimation**: Long-term studies suggest microbial communities may
           acclimate, reducing effective Q10 over time. Fixed Q10 may overestimate
           long-term feedback.

        4. **Depth dependence**: Deep soil carbon may have different temperature
           sensitivity than surface layers.
        """
        # Use direct exponential calculation for each effect
        return (
            np.exp(dT * self.sns_cLP2dT),
            np.exp(dT * self.sns_cLD2dT),
            np.exp(dT * self.sns_cSR2dT),
        )

    def calculate_eff_dT_PUBNF(
            self,
            dT: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity]:
        """
        Calculate temperature effects for nitrogen acquisition processes.

        Computes exponential temperature response factors for plant nitrogen
        uptake (PU) and biological nitrogen fixation (BNF). These factors modify
        nitrogen acquisition rates affecting ecosystem nitrogen supply.

        Parameters
        ----------
        dT : pint.Quantity
            Change in temperature from reference [K].

        Returns
        -------
        eff_dT_PU : pint.Quantity
            Temperature effect on plant nitrogen uptake [dimensionless].
        eff_dT_BNF : pint.Quantity
            Temperature effect on biological nitrogen fixation [dimensionless].

        Notes
        -----
        **Mathematical Formulation**

        Both effects use simple exponential responses:

        .. math::
           \\epsilon_{T(PU)} = \\exp(s_{PU2dT} \\times \\Delta T)

        .. math::
           \\epsilon_{T(BNF)} = \\exp(s_{BNF2dT} \\times \\Delta T)

        **Plant Nitrogen Uptake**

        Temperature affects N uptake through multiple mechanisms:

        1. **Root metabolic activity**: Active transport across membranes is
           temperature-dependent (enzyme kinetics)
        2. **Root growth**: Warmer soils promote root proliferation
        3. **Soil solution mobility**: Diffusion and mass flow increase with
           temperature
        4. **Mycorrhizal activity**: Symbiotic fungi are temperature-sensitive

        **Biological Nitrogen Fixation**

        Symbiotic N fixation (e.g., legumes + rhizobia, actinorhizal plants) is
        strongly temperature-dependent:

        1. **Nitrogenase enzyme kinetics**: Optimal temperature typically 20-30°C
        2. **Bacterial growth rates**: Exponential temperature dependence
        3. **Plant carbon allocation**: Warmer temperatures may increase carbon
           supply to symbionts
        4. **Nodule formation**: Temperature affects infection and nodulation

        **Interaction with Carbon-Nitrogen Coupling**

        These temperature effects interact with C-N coupling effects:

        - Enhanced N uptake/fixation (temperature) may partially alleviate N
          limitation of CO2-fertilized NPP
        - However, if warming also increases NPP (temperature effect), N demand
          increases, potentially intensifying N limitation
        - Net effect on N limitation depends on the balance between enhanced
          supply (uptake, fixation, mineralization) and enhanced demand (NPP)
        """
        return np.exp(dT * self.sns_PU2dT), np.exp(dT * self.sns_BNF2dT)

    def calculate_eff_dT_nLPLDSRLS(
            self,
            dT: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity, pint.Quantity, pint.Quantity]:
        """
        Calculate temperature effects for nitrogen cycle turnover processes.

        Computes exponential temperature response factors for nitrogen litter
        production (LP), litter decomposition (LD), soil respiration/mineralization
        (SR), and gaseous losses (LSgas). These factors modify nitrogen cycling
        rates affecting ecosystem nitrogen availability.

        Parameters
        ----------
        dT : pint.Quantity
            Change in temperature from reference [K].

        Returns
        -------
        eff_dT_nLP : pint.Quantity
            Temperature effect on nitrogen litter production [dimensionless].
        eff_dT_nLD : pint.Quantity
            Temperature effect on nitrogen litter decomposition [dimensionless].
        eff_dT_nSR : pint.Quantity
            Temperature effect on nitrogen mineralization (soil respiration)
            [dimensionless].
        eff_dT_nLSgas : pint.Quantity
            Temperature effect on gaseous nitrogen losses (denitrification,
            volatilization) [dimensionless].

        Notes
        -----
        **Mathematical Formulation**

        All effects use exponential responses:

        .. math::
           \\epsilon_{T(i)} = \\exp(s_{i2dT} \\times \\Delta T)

        **Nitrogen Mineralization (nSR): Supply of Plant-Available N**

        Nitrogen mineralization (conversion of organic N to mineral NH4+ and NO3-)
        is temperature-sensitive, typically with Q10 ~ 2.0-2.5:

        .. math::
           \\epsilon_{T(nSR)} = \\exp(s_{nSR2dT} \\times \\Delta T)

        This is a **critical process** controlling N availability:

        - Warming → Enhanced mineralization → More plant-available N
        - May alleviate progressive nitrogen limitation under CO2 fertilization
        - Often mirrors carbon soil respiration temperature sensitivity

        Typical sensitivities: :math:`s_{nSR2dT}` ~ 0.05-0.09 K⁻¹ (Q10 ~ 1.6-2.5).

        **Gaseous Nitrogen Losses (nLSgas): N Retention vs. Loss**

        Gaseous N losses include:

        1. **Denitrification** (NO3- → N2O → N2): Anaerobic bacterial process,
           very temperature-sensitive (Q10 often > 2.5)
        2. **Volatilization** (NH3 loss): Less temperature-sensitive but increases
           with warmth
        3. **Nitrification-denitrification** (coupled processes): Both steps
           temperature-dependent

        Denitrification is particularly important:

        - Occurs in waterlogged/anoxic microsites
        - Strong temperature dependence: :math:`s_{nLSgas2dT}` ~ 0.08-0.12 K⁻¹
          (Q10 ~ 2.2-3.3)
        - Produces N2O (potent greenhouse gas) as intermediate

        **Trade-off**: Warming increases both N mineralization (supply) and gaseous
        losses. Net effect on N availability depends on which process is more
        sensitive and on soil moisture/oxygen status.
        """
        return (
            np.exp(dT * self.sns_nLP2dT),
            np.exp(dT * self.sns_nLD2dT),
            np.exp(dT * self.sns_nSR2dT),
            np.exp(dT * self.sns_nLSgas2dT),
        )


@define
class EffectCarbonNitrogenCouplingCalculator:
    """
    Calculator of carbon-nitrogen coupling effects for ecological processes.

    Computes bidirectional feedback effects between carbon and nitrogen cycles,
    representing the fundamental co-limitation of ecosystems by carbon (energy)
    and nitrogen (nutrients). This calculator implements two complementary types
    of coupling:

    1. **Nitrogen limitation of carbon processes** (N → C effects): How nitrogen
       availability constrains NPP and carbon turnover
    2. **Carbon enhancement of nitrogen processes** (C → N effects): How carbon
       availability (via NPP) enhances plant uptake and fixation

    These coupling effects are critical for predicting ecosystem responses to
    CO2 fertilization and climate change, particularly the phenomenon of
    "progressive nitrogen limitation" where CO2-stimulated growth is increasingly
    constrained by nitrogen availability.

    **Switch Parameter**

    The ``switch_N`` parameter toggles nitrogen limitation effects:

    - ``switch_N = 0``: No N limitation (effects = 1.0)
    - ``switch_N ≠ 0``: N limitation active
    """

    sns_NPP2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """
        Sensitivity of net primary production to nitrogen plant uptake deficit [GtN/yr]
    """

    sns_cLP2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """
        Sensitivity of litter production (carbon) to nitrogen plant uptake deficit [GtN/yr]
    """

    sns_cLD2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """
        Sensitivity of litter decomposition (carbon) to nitrogen plant uptake deficit [GtN/yr]
    """

    sns_cSR2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """
        Sensitivity of soil respiration (carbon) to nitrogen plant uptake deficit [GtN/yr]
    """

    sns_PU2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """Sensitivity of plant uptake to N deficit.

    Controls how strongly plant uptake responds to nitrogen deficiency [yr/GtN].
    Higher values indicate stronger uptake response to N limitation.
    """

    sns_BNF2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """Sensitivity of biological nitrogen fixation to N deficit.

    Controls how strongly fixation responds to nitrogen deficiency [yr/GtN].
    Higher values indicate stronger fixation response to N limitation.
    """

    sns_nLP2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """
        Sensitivity of litter production (nitrogen) to nitrogen plant uptake deficit [GtN/yr]
    """

    sns_nLD2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """
        Sensitivity of litter decomposition (nitrogen) to nitrogen plant uptake deficit [GtN/yr]
        """

    sns_nSR2PUdef: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """
        Sensitivity of sol respiration (nitrogen) to nitrogen plant uptake deficit [GtN/yr]
    """

    eff_C_PU_max: pint.Quantity = field(
        default=Q(2, "1"),
        validator=check_units("1"),
    )
    """Maximum effect of NPP on plant uptake.

    Upper limit of the NPP-dependent multiplier on plant uptake rate [dimensionless].
    Values > 1 allow enhanced uptake under high NPP conditions.
    """

    sns_PU2NPPrd: pint.Quantity = field(
        default=Q(2.5, "1"),
        validator=check_units("1"),
    )
    """Sensitivity of plant uptake to NPP.

    Controls how strongly plant uptake responds to changes in NPP [dimensionless].
    Higher values indicate stronger coupling between NPP and uptake.
    """

    eff_C_BNF_max: pint.Quantity = field(
        default=Q(2, "1"),
        validator=check_units("1"),
    )
    """Maximum effect of NPP on biological nitrogen fixation.

    Upper limit of the NPP-dependent multiplier on fixation rate [dimensionless].
    Values > 1 allow enhanced fixation under high NPP conditions.
    """

    sns_BNF2NPPrd: pint.Quantity = field(
        default=Q(2.5, "1"),
        validator=check_units("1"),
    )
    """Sensitivity of biological nitrogen fixation to NPP.

    Controls how strongly fixation responds to changes in NPP [dimensionless].
    Higher values indicate stronger coupling between NPP and fixation.
    """

    # Reference Values
    NPP0: pint.Quantity = field(
        default=Q(50, "GtC/yr"),
        validator=check_units("GtC/yr"),
    )
    """Reference net primary production.

    Baseline NPP used to calculate relative changes in production [GtC/yr].
    Used to normalize NPP effects on various processes.
    """

    def calculate_eff_N_NPP(
            self,
            NflxPUdef: pint.Quantity,
            switch_N: pint.Quantity,
    ) -> pint.Quantity:
        """
        Calculate the carbon-nitrogen coupling effect for Net Primary Production (NPP).

        Parameters
        ----------
        NflxPUdef : pint.Quantity
            Nitrogen plant uptake deficit [GtN/yr].
        switch_N : pint.Quantity
            Switch to enable/disable nitrogen limitation effects [0 or 1].

        Returns
        -------
        pint.Quantity
            Carbon-nitrogen coupling effect on NPP.

        Notes
        -----
        The carbon-nitrogen coupling effect on NPP is calculated as:

        .. math::

            \\epsilon_{N(NPP)} = \\begin{cases}
              e^{-s_{NPP2PUdef} \\times PUdef} & \\text{if } \\mathrm{
              switch\\_N} \\neq 0 \\\\
              1 & \\text{if } \\mathrm{switch\\_N} = 0
           \\end{cases}
        """
        # Apply nitrogen limitation effect on NPP only when nitrogen switch is on
        return np.where(switch_N != 0, np.exp(-self.sns_NPP2PUdef * NflxPUdef), one)

    def calculate_eff_N_cLPLDSR(
            self,
            NflxPUdef: pint.Quantity,
            switch_N: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity, pint.Quantity]:
        """
        Calculate carbon-nitrogen coupling effects for carbon cycle processes i: LP,
        LD, and SR.

        Parameters
        ----------
        NflxPUdef : pint.Quantity
            Nitrogen plant uptake deficit [GtN/yr].
        switch_N : pint.Quantity
            Switch to enable/disable nitrogen limitation effects [0 or 1].

        Returns
        -------
        eff_N_cLP : pint.Quantity
            N effect on carbon litter production [dimensionless].
        eff_N_cLD : pint.Quantity
            N effect on carbon litter decomposition [dimensionless].
        eff_N_cSR : pint.Quantity
            N effect on carbon soil respiration [dimensionless].

        Notes
        -----
        Each effect is calculated as:

        .. math::
           \\epsilon_{N(i)} = \\begin{cases}
              e^{s_{i2PUdef} \\times PUdef} & \\text{if } \\mathrm{switch\\_N}
              \\neq 0 \\\\
              1 & \\text{if } \\mathrm{switch\\_N} = 0
           \\end{cases}
        """
        # Apply nitrogen effects on carbon processes only when nitrogen switch is on
        return (
            np.where(switch_N != 0, np.exp(self.sns_cLP2PUdef * NflxPUdef), one),
            np.where(switch_N != 0, np.exp(self.sns_cLD2PUdef * NflxPUdef), one),
            np.where(switch_N != 0, np.exp(self.sns_cSR2PUdef * NflxPUdef), one),
        )

    def calculate_eff_C_PUBNF(
            self,
            CflxNPP: pint.Quantity,
            CflxNPP0: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity]:
        """
        Calculate carbon availability effects on nitrogen acquisition processes.

        Parameters
        ----------
        CflxNPP : pint.Quantity
            Current Net Primary Production flux [GtC/yr].
        CflxNPP0 : pint.Quantity
            Reference Net Primary Production flux [GtC/yr].

        Returns
        -------
        eff_C_PU : pint.Quantity
            Carbon enhancement of plant N uptake [dimensionless].
        eff_C_BNF : pint.Quantity
            Carbon enhancement of biological N fixation [dimensionless].

        Notes
        -----
        Both effects use sigmoid formulations based on relative NPP change:

        .. math::
           NPPrd = \\Delta NPPrel = \\frac{NPP}{NPP_0} - 1

        .. math::
           \\epsilon_{C(PU)} = \\frac{\\epsilon_{C(PU)_{max}}}{1 + (\\epsilon_{C(
           PU)_{max}} - 1)
           \\times e^{-s_{PU2NPPrd} \\times NPPrd}}

        .. math::
           \\epsilon_{C(BNF)} = \\frac{\\epsilon_{C(BNF)_{max}}}{1 + (\\epsilon_{C(
           BNF)_{max}} - 1)
           \\times e^{-s_{BNF2NPPrd} \\times NPPrd}}

        The sigmoid form ensures the effect equals 1 at :math:`NPP = NPP_0` and
        approaches :math:`\\epsilon_{max}` at high NPP, representing saturation
        of the carbon enhancement effect.
        """
        CflxNPPrd = CflxNPP / CflxNPP0 - 1

        return (
                self.eff_C_PU_max
                / (1 + (self.eff_C_PU_max - 1) * np.exp(-self.sns_PU2NPPrd * CflxNPPrd))
        ), (
                self.eff_C_BNF_max
                / (1 + (self.eff_C_BNF_max - 1) * np.exp(
            -self.sns_BNF2NPPrd * CflxNPPrd))
        )

    def calculate_eff_N_PUBNF(
            self,
            NflxPUdef: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity]:
        """
        Calculate nitrogen deficit effects on nitrogen acquisition processes.

        Parameters
        ----------
        NflxPUdef : pint.Quantity
            Nitrogen plant uptake deficit [GtN/yr].

        Returns
        -------
        eff_N_PU : pint.Quantity
            N deficit effect on plant N uptake [dimensionless].
        eff_N_BNF : pint.Quantity
            N deficit effect on biological N fixation [dimensionless].

        Notes
        -----
        Both effects use exponential responses with positive exponents:

        .. math::
           \\epsilon_{N(PU)} = \\exp(s_{PU2PUdef} \\times PUdef)

        .. math::
           \\epsilon_{N(BNF)} = \\exp(s_{BNF2PUdef} \\times PUdef)
        """
        return (
            np.exp(self.sns_PU2PUdef * NflxPUdef),
            np.exp(self.sns_BNF2PUdef * NflxPUdef),
        )

    def calculate_eff_N_nLPLDSR(
            self,
            NflxPUdef: pint.Quantity,
    ) -> Tuple[pint.Quantity, pint.Quantity, pint.Quantity]:
        """
        Calculate carbon-nitrogen coupling effects for nitrogen cycle processes: LP, LD, and SR

        Parameters
        ----------
        NflxPUdef
            Nitrogen plant uptake deficit [GtN/yr]

        Returns
        -------
        eff_N_nLP : pint.Quantity
            N deficit effect on nitrogen litter production [dimensionless].
        eff_N_nLD : pint.Quantity
            N deficit effect on nitrogen litter decomposition [dimensionless].
        eff_N_nSR : pint.Quantity
            N deficit effect on nitrogen mineralization [dimensionless].

        Notes
        -----
        All effects use exponential responses with positive exponents:

        .. math::
           \\epsilon_{N(i)} = \\exp(s_{i2PUdef} \\times PUdef)
        """
        return (
            np.exp(self.sns_nLP2PUdef * NflxPUdef),
            np.exp(self.sns_nLD2PUdef * NflxPUdef),
            np.exp(self.sns_nSR2PUdef * NflxPUdef),
        )
