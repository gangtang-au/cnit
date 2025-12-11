"""
MAGICC Coupled Carbon-Nitrogen Cycle Model (CNit) Implementation
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

            CflxNPP = (NPP0 + CflxLUrgr) \times eff\_dT\_NPP \times eff\_CO2\_NPP \times eff\_CN\_NPP \times eff\_LU\_NPP

            CflxLPR = LPR0 \times eff\_dT\_LPR \times eff\_CO2\_NPP \times eff\_CN\_NPP \times eff\_LU\_NPP
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
    """
    Calculator of carbon turnover for plant, litter, and soil carbon pools
    """

    tau_CplsP: pint.Quantity = field(
        default=Q(10, "yr"), validator=check_units("yr")
    )
    """
        Plant carbon pool turnover time [yr]
    """

    tau_CplsL: pint.Quantity = field(
        default=Q(1, "yr"), validator=check_units("yr")
    )
    """
        Litter carbon pool turnover time [yr]
    """

    tau_CplsS: pint.Quantity = field(
        default=Q(100, "yr"), validator=check_units("yr")
    )
    """
        Soil carbon pool turnover time [yr]
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
        Calculate carbon turnover fluxes from each carbon pool based on turnover times and effect factors.

        Parameters
        ----------
        CplsP : pint.Quantity
            Plant carbon pool size [GtC].
        CplsL : pint.Quantity
            Litter carbon pool size [GtC].
        CplsS : pint.Quantity
            Soil carbon pool size [GtC].
        eff_dT_cLP : pint.Quantity
            Temperature effect on plant-to-litter turnover [1].
        eff_dT_cLD : pint.Quantity
            Temperature effect on litter decomposition [1].
        eff_dT_cSR : pint.Quantity
            Temperature effect on soil respiration [1].
        eff_N_cLP : pint.Quantity
            CN effect on plant-to-litter turnover [1].
        eff_N_cLD : pint.Quantity
            CN effect on litter decomposition [1].
        eff_N_cSR : pint.Quantity
            CN effect on soil respiration [1].
        eff_LU_cLP : pint.Quantity
            Land use effect on plant-to-litter turnover [1].
        eff_LU_cLD : pint.Quantity
            Land use effect on litter decomposition [1].
        eff_LU_cSR : pint.Quantity
            Land use effect on soil respiration [1].

        Returns
        -------
            Dictionary with keys:

            - CflxLP: carbon flux from plant to litter [GtC/yr]
            - CflxLD: carbon flux from litter to soil [GtC/yr]
            - CflxSR: carbon flux from soil through respiration [GtC/yr]

        Notes
        -----
        The following formulas are used:

        .. math::

            CflxLP = \\frac{\\mathrm{CplsP}}{\\tau_{\\mathrm{CplsP}}} \\times \\mathrm{eff\\_dT\\_cLP} \\times \\mathrm{eff\\_N\\_cLP} \\times \\mathrm{eff\\_LU\\_cLP}

            CflxLD = \\frac{\\mathrm{CplsL}}{\\tau_{\\mathrm{CplsL}}} \\times \\mathrm{eff\\_dT\\_cLD} \\times \\mathrm{eff\\_N\\_cLD} \\times \\mathrm{eff\\_LU\\_cLD}

            CflxSR = \\frac{\\mathrm{CplsS}}{\\tau_{\\mathrm{CplsS}}} \\times \\mathrm{eff\\_dT\\_cSR} \\times \\mathrm{eff\\_N\\_cSR} \\times \\mathrm{eff\\_LU\\_cSR}
        """
        return {
            "CflxLP": CplsP / self.tau_CplsP * eff_dT_cLP * eff_N_cLP * eff_LU_cLP,
            "CflxLD": CplsL / self.tau_CplsL * eff_dT_cLD * eff_N_cLD * eff_LU_cLD,
            "CflxSR": CplsS / self.tau_CplsS * eff_dT_cSR * eff_N_cSR * eff_LU_cSR,
        }


@define
class CarbonCycleCalculator:
    """
    Calculator for carbon cycle fluxes and pool changes
    """

    frc_NPP2P: pint.Quantity = field(
        default=Q(0.55, "1"),
        validator=check_units("1"),
    )
    """
        Fraction of NPP which goes into the plant pool [1]
    """

    frc_NPP2L: pint.Quantity = field(
        default=Q(0.35, "1"),
        validator=check_units("1"),
    )
    """
        Fraction of NPP which goes into the litter pool [1]
    """

    frc_cLP2L: pint.Quantity = field(
        default=Q(0.7, "1"),
        validator=check_units("1"),
    )
    """
        Fraction of litter production which goes into the litter pool [1]
        The rest goes into the soil pool
    """

    frc_cLD2S: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """
        Fraction of litter decomposition which goes into the soil pool [1]
        The rest is released directly into the atmosphere
    """

    frc_cLUgrs2P: pint.Quantity = field(
        default=Q(0.6, "1"),
        validator=check_units("1"),
    )
    """
        Fraction of land use emissions which is assumed to occur in the plant pool [1]
    """

    frc_cLUgrs2L: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """
        Fraction of land use emissions which is assumed to occur in the litter pool [1]
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
        Calculate carbon pool size changes over time

        Parameters
        ----------
        CflxNPP
            Net Primary Production carbon flux [GtC/yr]
        CflxLPR
            Litter Production Respiration carbon flux [GtC/yr]
        CflxLP
            Litter Production carbon flux [GtC/yr]
        CflxLD
            Litter Decomposition carbon flux [GtC/yr]
        CflxSR
            Soil Respiration carbon flux [GtC/yr]
        CflxLUgrs
            Gross Land Use carbon flux [GtC/yr]

        Returns
        -------
            List containing:

            - dCplsP_dt: Rate of change of plant carbon pool [GtC/yr]
            - dCplsL_dt: Rate of change of litter carbon pool [GtC/yr]
            - dCplsS_dt: Rate of change of soil carbon pool [GtC/yr]

        Notes
        -----
        The mass balance equations used are:

        .. math::

            \\frac{dC_{plsP}}{dt} &= C_{flxNPP2P} - C_{flxLPR} - C_{flxLP} - C_{flxLUgrs2P} \\\\
            \\frac{dC_{plsL}}{dt} &= C_{flxNPP2L} + C_{flxLP2L} - C_{flxLD} - C_{flxLUgrs2L} \\\\
            \\frac{dC_{plsS}}{dt} &= C_{flxNPP2S} + C_{flxLP2S} + C_{flxLD2S} - C_{flxSR} - C_{flxLUgrs2S}
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
        Calculate all carbon fluxes in the system

        Parameters
        ----------
        CflxNPP
            Net Primary Production carbon flux [GtC/yr]
        CflxLPR
            Litter Production Respiration carbon flux [GtC/yr]
        CflxLP
            Litter Production carbon flux [GtC/yr]
        CflxLD
            Litter Decomposition carbon flux [GtC/yr]
        CflxSR
            Soil Respiration carbon flux [GtC/yr]
        CflxLUgrs
            Gross Land Use carbon flux [GtC/yr]

        Returns
        -------
            Dictionary containing all carbon fluxes in the system
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
    """Results of running the carbon cycle model.

    Attributes:
        ts_collection: Collection of time series results from the run
        calc_c_turnover: Calculator used for carbon turnover
        calc_c_cycle: Calculator used for carbon cycle
        time_axis: Time points for simulation [yr]
        CflxNPP_t: Time-interpolated net primary production flux [GtC/yr]
        CflxLPR_t: Time-interpolated litter production flux [GtC/yr]
        CflxLUgrs_t: Time-interpolated land use gross flux [GtC/yr]
        eff_dT_cLP_t: Time-interpolated temperature effect on litter production [1]
        eff_dT_cLD_t: Time-interpolated temperature effect on litter decomposition [1]
        eff_dT_cSR_t: Time-interpolated temperature effect on soil respiration [1]
        eff_N_cLP_t: Time-interpolated CN effect on litter production [1]
        eff_N_cLD_t: Time-interpolated CN effect on litter decomposition [1]
        eff_N_cSR_t: Time-interpolated CN effect on soil respiration [1]
        eff_LU_cLP_t: Time-interpolated land use effect on litter production [1]
        eff_LU_cLD_t: Time-interpolated land use effect on litter decomposition [1]
        eff_LU_cSR_t: Time-interpolated land use effect on soil respiration [1]

    Note:
        All *_t functions are time-continuous series created through linear interpolation
        of discrete time points. These functions map from time [yr] to their respective units.
    """

    c_state: Dict[str, pint.Quantity]
    calc_c_turnover: CarbonTurnoverCalculator
    calc_c_cycle: CarbonCycleCalculator
    time_axis: pint.Quantity

    # Time-interpolated flux functions
    CflxNPP_t: Callable[[pint.Quantity], pint.Quantity]
    CflxLPR_t: Callable[[pint.Quantity], pint.Quantity]
    CflxLUgrs_t: Callable[[pint.Quantity], pint.Quantity]

    # Time-interpolated effect functions
    eff_dT_cLP_t: Callable[[pint.Quantity], pint.Quantity]
    eff_dT_cLD_t: Callable[[pint.Quantity], pint.Quantity]
    eff_dT_cSR_t: Callable[[pint.Quantity], pint.Quantity]
    eff_N_cLP_t: Callable[[pint.Quantity], pint.Quantity]
    eff_N_cLD_t: Callable[[pint.Quantity], pint.Quantity]
    eff_N_cSR_t: Callable[[pint.Quantity], pint.Quantity]
    eff_LU_cLP_t: Callable[[pint.Quantity], pint.Quantity]
    eff_LU_cLD_t: Callable[[pint.Quantity], pint.Quantity]
    eff_LU_cSR_t: Callable[[pint.Quantity], pint.Quantity]

    def add_non_state_variables(self) -> xr.Dataset:
        """
        Add non-state variables to the results.

        Returns
        -------
        xr.Dataset
            Results including derived variables.
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
    """MAGICC's terrestrial carbon cycle model.

    This model simulates the dynamics of carbon pools in the terrestrial biosphere.
    Time-continuous functions (suffixed with _t) are created through linear interpolation
    of discrete time points and provide continuous values over the simulation period.

    Attributes:
        calc_c_turnover: Calculator for carbon turnover processes
        calc_c_cycle: Calculator for carbon cycle processes
        CplsP0: Initial plant carbon pool size [GtC]
        CplsL0: Initial litter carbon pool size [GtC]
        CplsS0: Initial soil carbon pool size [GtC]
        time0: Initial time point [yr]
        switch_Cpls: Switches for enabling/disabling carbon pools [1]
    """

    calc_c_turnover: CarbonTurnoverCalculator
    calc_c_cycle: CarbonCycleCalculator

    CplsP0: pint.Quantity = field(
        default=Q(400, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial plant carbon pool size [GtC]"""

    CplsL0: pint.Quantity = field(
        default=Q(50, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial litter carbon pool size [GtC]"""

    CplsS0: pint.Quantity = field(
        default=Q(1550, "GtC"),
        validator=check_units("GtC"),
    )
    """Initial soil carbon pool size [GtC]"""

    time0: pint.Quantity = field(
        default=Q(1850, "yr"),
        validator=check_units("yr"),
    )
    """Time at which the initialisation values apply [yr]"""

    switch_Cpls: list[int] = field(
        default=[1, 1, 1],
    )
    """Carbon pool switches (1=enabled, 0=disabled) for [plant, litter, soil]"""

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

        Parameters
        ----------
        time_axis : pint.Quantity
            Time points for simulation [yr].
        CflxNPP_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated net primary production flux [GtC/yr].
        CflxLPR_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated litter production flux [GtC/yr].
        CflxLUgrs_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use gross flux [GtC/yr].
        eff_dT_cLP_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated temperature effect on litter production [1].
        eff_dT_cLD_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated temperature effect on litter decomposition [1].
        eff_dT_cSR_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated temperature effect on soil respiration [1].
        eff_N_cLP_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated CN effect on litter production [1].
        eff_N_cLD_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated CN effect on litter decomposition [1].
        eff_N_cSR_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated CN effect on soil respiration [1].
        eff_LU_cLP_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use effect on litter production [1].
        eff_LU_cLD_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use effect on litter decomposition [1].
        eff_LU_cSR_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use effect on soil respiration [1].

        Returns
        -------
        CarbonCycleModelResult
            Simulation results.

        Raises
        ------
        ValueError
            If `time_axis` starts before `time0`.
        SolveError
            If the numerical solution fails to converge.

        Notes
        -----
        All *_t arguments are time-continuous series created through linear interpolation
        of discrete time points. These arguments are functions mapping from time [yr]
        to their respective units.
        """

        def func_to_solve(t: float, y: np.ndarray) -> list:
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
            dydt = [v.to("GtC/yr").m * switch for v, switch in zip(dydt, self.switch_Cpls)]
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
