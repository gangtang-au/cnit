"""
Land Use Flux Calculator for the CNit Model.

This module defines the `LanduseCalculator` class, which is responsible
for processing driver time series of net or gross land-use emissions
(for both Carbon and Nitrogen) and translating them into the specific
fluxes required by the model, such as deforestation, afforestation,
and regrowth fluxes.
"""

from typing import Dict
from attrs import define, field
import numpy as np
import pint

from ..utils.units import check_units, Q


@define
class LanduseCalculator:
    """
    Calculator for land use related fluxes.
    """

    frc_rgr: pint.Quantity = field(
        default=Q(0.9, "1"),
        validator=check_units("1"),
    )
    """
    Fraction of deforestation that can regrow [dimensionless].
    """

    tau_rgr: pint.Quantity = field(
        default=Q(100, "yr"),
        validator=check_units("yr"),
    )
    """
    Time for regrowth (years).
    """

    method_cLU: pint.Quantity = field(
        default=Q(0, "1"),
        validator=check_units("1"),
    )
    """
    Carbon flux calculation method [dimensionless]:

        - 0: Use net emission (CemsLUnet) as input (e.g., AFOLU emissions).
        - 1: Use gross emission (CemsLUgrs) as input (e.g., gross deforestation).
    """

    method_nLU: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """
    Nitrogen flux calculation method [dimensionless]:

        - 0: Use net emission (NemsLUnet) as input (not commonly used).
        - 1: Use gross emission (NemsLUgrs) as input (e.g., gross deforestation).
    """

    def calculate_CflxLU_series(
            self,
            CemsLUnet_s: pint.Quantity,
            CemsLUgrs_s: pint.Quantity,
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate carbon flux series from land use changes.

        This method computes various time series of carbon fluxes associated with land
        use change, including gross flux, regrowth, deforestation, afforestation,
        and decay after afforestation, using the specified calculation method.

        Parameters
        ----------
        CemsLUnet_s : pint.Quantity
            Net land use carbon emission time series.
        CemsLUgrs_s : pint.Quantity
            Gross land use carbon emission time series.

        Returns
        -------
            Dictionary with time series for each flux:
                - CflxLUgrs: Gross land use carbon flux
                - CflxLUrgr: Regrowth carbon flux
                - CflxLUdfst: Deforestation carbon flux
                - CflxLUafst: Afforestation carbon flux
                - CflxLUafst_decay: Afforestation decay carbon flux
                - cumsum_CflxLUrgr: Cumulative regrowth carbon
                - cumsum_CflxLUdfst: Cumulative deforestation carbon
                - cumsum_CflxLUafst: Cumulative afforestation carbon
                - cumsum_CflxLUafst_decay: Cumulative afforestation decay carbon

        Notes
        -----
        For the land use emissions (for both carbon and nitrogen), the net flux (LUnet)
        is the difference between the gross emission (LUgrs) and the regrowth (
        LUrgr, i.e., recovery after deforestation).
            .. math::
                LUnet(t) = LUgrs(t) - LUrgr(t)

        The gross flux (LUgrs) is separated into deforestation (LUdfst) and
        afforestation (LUafst), both of which are always positive values:
            .. math::
                LUnet(t) = LUdfst(t) - LUrgr(t) - LUafst(t)

        The separation is especially important when using LUnet as input and when the
        LUnet is negative, the simulated LUgrs = LUnet + LUrgr might be negative, which
        essentially means afforestation or reforestation. In this case,
        we set LUdfst = 0 and get the LUafst = -LUgrs.

        We acknowledge the LUgrs might be confusing and misleading because in the
        common sense the "gross emission" is always positive, but here we
        keep the LUgrs as

        For a single deforestation event at time t, the regrowth flux is a continued
        constant flux over the regrowth time period (tau_rgr), calculated as:
            .. math::
                LUrgr(T) = \\frac{\\varphi \\times LUdfst(t)}{\\tau_{rgr}}, \\forall
                T\\in [t+1, t+\\tau_{rgr}]

        where :math:`\\varphi` is the fraction of deforestation that can regrow.

        The total regrowth flux at time t is the sum of contributions from all past
        deforestation events that are still regrowing:
            .. math::
                LUrgr(t) = \\begin{cases}
                    \\sum_{T = \\max(0, t - \\tau_{rgr})}^{t-1} \\frac{\\varphi
                    \\times LUdfst(T)}{\\tau_{rgr}} & t \\geq 1 \\\\
                    0 & t < 1 \\\\
                \\end{cases}


        The carbon fluxes are calculated based on the chosen method:

        - Method 0: Uses net emissions (CemsLUnet) as input, suitable for AFOLU
        emissions.

        - Method 1: Uses gross emissions (CemsLUgrs) as input, suitable for gross
        deforestation data.
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
        Calculate nitrogen flux series from land use changes.

        This method computes various time series of nitrogen fluxes associated with land
        use change, including gross flux, regrowth, deforestation, afforestation,
        afforestation decay, using the specified calculation method. Nitrogen emitted
        directly from the mineral pool (e.g., most AFOLU nitrogen emissions of NH₃, NOₓ,
        and N₂O from agricultural land use change) is included as a separate flux.


        Parameters
        ----------
        NemsLUnet_s : pint.Quantity
            Net land use nitrogen emission time series.
        NemsLUgrs_s : pint.Quantity
            Gross land use nitrogen emission time series.
        NemsLUmin_s : pint.Quantity
            Land use nitrogen emission directly from mineral pool time series.

        Returns
        -------
            Dictionary with time series for each flux:
                - NflxLUgrs: Gross land use flux
                - NflxLUrgr: Regrowth flux
                - NflxLUdfst: Deforestation flux
                - NflxLUafst: Afforestation flux
                - NflxLUafst_decay: Afforestation decay flux
                - cumsum_NflxLUrgr: Cumulative regrowth
                - cumsum_NflxLUdfst: Cumulative deforestation
                - cumsum_NflxLUafst: Cumulative afforestation
                - cumsum_NflxLUafst_decay: Cumulative afforestation decay
                - NflxLUmin: Land use nitrogen emission directly from mineral pool
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
        Calculates and returns flux series with appropriate naming.

        Parameters
        ----------
        net_flux : array-like
            Net flux time series.
        gross_flux : array-like
            Gross flux time series.
        method_calc : int
            Method for calculation (0 or 1).
        unit : str
            Unit string for the quantities.
        prefix : str
            Prefix for the returned dictionary keys (e.g., "Cflx" or "Nflx").

        Returns
        -------
        dict
            Dictionary with named flux series results.

        Raises
        ------
        ValueError
            If method_calc is not 0 or 1.
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
