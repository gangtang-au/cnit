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
from ..utils.maths import _get_interp, SolveError
from ..utils.data_utils import make_dataset_from_var_dict


@define
class NitrogenPUBNFCalculator:
    """Calculator for nitrogen plant uptake and biological nitrogen fixation.

    This calculator determines nitrogen fluxes related to plant uptake (PU) and
    biological nitrogen fixation (BNF) in the terrestrial biosphere.

    Attributes grouped by function:

    Plant Uptake parameters:
        PU0: Initial plant uptake rate [GtN/yr]
        eff_C_PU_max: Maximum NPP effect on plant uptake [1]
        sns_PU2NPPrd: Sensitivity of PU to NPP [1]
        sns_PU2PUdef: Sensitivity of PU to N deficit [yr/GtN]

    Biological Nitrogen Fixation parameters:
        BNF0: Initial BNF rate [GtN/yr]
        eff_C_BNF_max: Maximum NPP effect on BNF [1]
        sns_BNF2NPPrd: Sensitivity of BNF to NPP [1]
        sns_BNF2PUdef: Sensitivity of BNF to N deficit [yr/GtN]
        constant_BNF: Switch for using constant BNF (True/False) [1]

    Net Mineralization parameters:
        NetMIN0: Initial net mineralization rate [GtN/yr]
        NetMINlt0: Long-term net mineralization rate [GtN/yr]
        sns_NetMIN2dT: Temperature sensitivity of net mineralization [1/K]
        sns_NetMINlt2AD: Sensitivity of net mineralization to N deposition [yr/GtN]

    Stoichiometric Ratios:
        NCratio_PUreq: Required N:C ratio for plant uptake [GtN/GtC]

    Reference Values:
        NPP0: Reference net primary production [GtC/yr]

    Last modified: 2025-04-27 11:22:52 UTC
    Author: gangtang-au
    """

    PU0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Initial plant uptake rate.

    Base rate of nitrogen uptake by plants under reference conditions [GtN/yr].
    Used to scale the actual uptake rate based on environmental conditions.
    """

    BNF0: pint.Quantity = field(
        default=Q(0.1, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Initial biological nitrogen fixation rate.

    Base rate of nitrogen fixation under reference conditions [GtN/yr].
    Used to scale actual fixation rate based on environmental conditions.
    """

    # Stoichiometric Ratios
    NCratio_PUreq: pint.Quantity = field(
        default=Q(0.02, "GtN/GtC"),
        validator=check_units("GtN/GtC"),
    )
    """Required nitrogen to carbon ratio for plant uptake.

    Stoichiometric ratio determining N demand based on C uptake [GtN/GtC].
    Used to calculate required N uptake from NPP.
    """

    NetMIN0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Initial net mineralization rate.

    Base rate of nitrogen mineralization under reference conditions [GtN/yr].
    Part of the temperature-sensitive component of mineralization.
    """

    sns_NetMIN2dT: pint.Quantity = field(
        default=Q(0, "1/K"),
        validator=check_units("1/K"),
    )
    """Temperature sensitivity of net mineralization.

    Controls how strongly mineralization responds to temperature changes [1/K].
    Higher values indicate stronger temperature effect on mineralization.
    """

    NetMINlt0: pint.Quantity = field(
        default=Q(0.5, "GtN/yr"),
        validator=check_units("GtN/yr"),
    )
    """Long-term net mineralization rate.

    Baseline mineralization rate under steady state [GtN/yr].
    Part of the deposition and NPP-sensitive component of mineralization.
    """

    sns_NetMINlt2AD: pint.Quantity = field(
        default=Q(0, "yr/GtN"),
        validator=check_units("yr/GtN"),
    )
    """Sensitivity of net mineralization to N deposition.

    Controls how strongly mineralization responds to N deposition [yr/GtN].
    Higher values indicate stronger deposition effect on mineralization.
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
        eff_dT_PU : pint.Quantity
            Temperature effect on plant uptake [1].
        eff_C_PU : pint.Quantity
            Carbon cycle effect on plant uptake [1].
        eff_N_PU : pint.Quantity
            Nitrogen cycle effect on plant uptake [1].
        eff_dT_BNF : pint.Quantity
            Temperature effect on biological nitrogen fixation [1].
        eff_C_BNF : pint.Quantity
            Carbon cycle effect on biological nitrogen fixation [1].
        eff_N_BNF : pint.Quantity
            Nitrogen cycle effect on biological nitrogen fixation [1].

        Returns
        -------
        dict of str to pint.Quantity
            Dictionary containing:
            - 'NflxPU': Plant uptake nitrogen flux [GtN/yr]
            - 'NflxBNF': Biological nitrogen fixation flux [GtN/yr]
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
        Calculate required and available nitrogen for plant uptake.

        Parameters
        ----------
        CflxNPP : pint.Quantity
            Net primary production carbon flux [GtC/yr].
        NflxAD : pint.Quantity
            Nitrogen deposition flux [GtN/yr].
        NflxFT : pint.Quantity
            Nitrogen fertilizer flux [GtN/yr].
        dT : pint.Quantity
            Temperature change [K].

        Returns
        -------
        dict of str to pint.Quantity
            Dictionary containing:
            - 'NflxPUreq': Required nitrogen for plant uptake [GtN/yr]
            - 'NflxPUavail': Available nitrogen for plant uptake [GtN/yr]
            - 'NflxPUavail_netMIN': Available nitrogen from net mineralization [GtN/yr]
            - 'NflxPUdef': Nitrogen deficit [GtN/yr]
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
    environmental effects.

    Attributes:
        tau_NplsP: Plant nitrogen pool turnover time [yr]
            Controls the rate of nitrogen transfer from plants to litter

        tau_NplsL: Litter nitrogen pool turnover time [yr]
            Controls the rate of nitrogen transfer from litter to soil

        tau_NplsS: Soil nitrogen pool turnover time [yr]
            Controls the rate of nitrogen release from soil through respiration

        tau_NplsM: Mineral nitrogen pool turnover time [yr]
            Controls the rate of gaseous nitrogen loss from mineral pool

        frc_nLSgas: Fraction of mineral nitrogen loss that is gaseous [1]
            Partitions mineral nitrogen loss between gaseous and non-gaseous pathways

    Methods:
        calculate: Computes nitrogen fluxes between pools based on current states
                  and environmental effects
    """

    # Pool turnover times
    tau_NplsP: pint.Quantity = field(
        default=Q(1, "yr"),
        validator=check_units("yr"),
    )
    """Plant nitrogen pool turnover time [yr]"""

    tau_NplsL: pint.Quantity = field(
        default=Q(1, "yr"),
        validator=check_units("yr"),
    )
    """Litter nitrogen pool turnover time [yr]"""

    tau_NplsS: pint.Quantity = field(
        default=Q(1, "yr"),
        validator=check_units("yr"),
    )
    """Soil nitrogen pool turnover time [yr]"""

    tau_NplsM: pint.Quantity = field(
        default=Q(1, "yr"),
        validator=check_units("yr"),
    )
    """Mineral nitrogen pool turnover time (for gaseous loss) [yr]"""

    # Fractionation parameters
    frc_nLSgas: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """Fraction of mineral nitrogen loss that is gaseous [1]"""

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
        Calculate nitrogen turnover fluxes from each nitrogen pool based on
        turnover times and effect factors.

        Parameters
        ----------
        NplsP : pint.Quantity
            Plant nitrogen pool size [GtN].
        NplsL : pint.Quantity
            Litter nitrogen pool size [GtN].
        NplsS : pint.Quantity
            Soil nitrogen pool size [GtN].
        NplsM : pint.Quantity
            Mineral nitrogen pool size [GtN].
        eff_dT_nLP : pint.Quantity
            Temperature effect on plant-to-litter turnover [1].
        eff_dT_nLD : pint.Quantity
            Temperature effect on litter decomposition [1].
        eff_dT_nSR : pint.Quantity
            Temperature effect on soil turnover [1].
        eff_N_nLP : pint.Quantity
            CN effect on plant-to-litter turnover [1].
        eff_N_nLD : pint.Quantity
            CN effect on litter decomposition [1].
        eff_N_nSR : pint.Quantity
            CN effect on soil turnover [1].
        eff_dT_nLSgas : pint.Quantity
            Temperature effect on gaseous nitrogen loss from mineral pool [1].
        eff_LU_nLP : pint.Quantity
            Land use effect on plant-to-litter turnover [1].
        eff_LU_nLD : pint.Quantity
            Land use effect on litter decomposition [1].
        eff_LU_nSR : pint.Quantity
            Land use effect on soil turnover [1].

        Returns
        -------
        dict of str to pint.Quantity
            Dictionary containing:
            - 'NflxLP': Nitrogen flux from plant to litter [GtN/yr]
            - 'NflxLD': Nitrogen flux from litter to soil [GtN/yr]
            - 'NflxSR': Nitrogen flux from soil through respiration [GtN/yr]
            - 'NflxLS': Nitrogen loss from mineral pool [GtN/yr]
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
        Calculate nitrogen turnover fluxes from each nitrogen pool based on
        turnover times and effect factors.

        Parameters
        ----------
        NplsP : pint.Quantity
            Plant nitrogen pool size [GtN].
        NplsL : pint.Quantity
            Litter nitrogen pool size [GtN].
        NplsS : pint.Quantity
            Soil nitrogen pool size [GtN].
        NplsM : pint.Quantity
            Mineral nitrogen pool size [GtN].
        eff_dT_nLP : pint.Quantity
            Temperature effect on plant-to-litter turnover [1].
        eff_dT_nLD : pint.Quantity
            Temperature effect on litter decomposition [1].
        eff_dT_nSR : pint.Quantity
            Temperature effect on soil turnover [1].
        eff_N_nLP : pint.Quantity
            CN effect on plant-to-litter turnover [1].
        eff_N_nLD : pint.Quantity
            CN effect on litter decomposition [1].
        eff_N_nSR : pint.Quantity
            CN effect on soil turnover [1].
        eff_dT_nLSgas : pint.Quantity
            Temperature effect on gaseous nitrogen loss from mineral pool [1].
        eff_LU_nLP : pint.Quantity
            Land use effect on plant-to-litter turnover [1].
        eff_LU_nLD : pint.Quantity
            Land use effect on litter decomposition [1].
        eff_LU_nSR : pint.Quantity
            Land use effect on soil turnover [1].

        Returns
        -------
        dict of str to pint.Quantity
            Dictionary containing:
            - 'NflxLP': Nitrogen flux from plant to litter [GtN/yr]
            - 'NflxLD': Nitrogen flux from litter to soil [GtN/yr]
            - 'NflxSR': Nitrogen flux from soil through respiration [GtN/yr]
            - 'NflxLS': Nitrogen loss from mineral pool [GtN/yr]
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
    """
    Calculator of nitrogen cycle.

    This class handles calculations related to nitrogen cycle processes including:
    - Biological nitrogen fixation (BNF)
    - Plant uptake (PU)
    - Litter production and decomposition
    - Soil processes
    - Land use changes:
        * Gross deforestation (NflxLUgrs)
        * Direct mineral nitrogen loss (NflxLUmin)

    Created by: gangtang-au
    Last updated: 2025-04-27 11:43:32 UTC
    """

    # Biological Nitrogen Fixation (BNF) fractions
    frc_BNF2P: pint.Quantity = field(
        default=Q(0.55, "1"),
        validator=check_units("1"),
    )
    """Fraction of nitrogen from BNF that goes to the plant pool [1]"""

    frc_BNF2L: pint.Quantity = field(
        default=Q(0.35, "1"),
        validator=check_units("1"),
    )
    """Fraction of nitrogen from BNF that goes to the litter pool [1]"""

    # Plant Uptake (PU) fractions
    frc_PU2P: pint.Quantity = field(
        default=Q(0.55, "1"),
        validator=check_units("1"),
    )
    """Fraction of nitrogen uptake flux which goes to the plant pool [1]"""

    frc_PU2L: pint.Quantity = field(
        default=Q(0.35, "1"),
        validator=check_units("1"),
    )
    """Fraction of nitrogen uptake flux which goes to the litter pool [1]"""

    # Litter fractions
    frc_nLP2L: pint.Quantity = field(
        default=Q(0.7, "1"),
        validator=check_units("1"),
    )
    """Fraction of litter production nitrogen which goes to the litter pool [1]"""

    frc_nLD2S: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """Fraction of litter decomposition nitrogen which goes to the soil pool [1]"""

    frc_nLSgas: pint.Quantity = field(
        default=Q(1, "1"),
        validator=check_units("1"),
    )
    """Fraction of mineral nitrogen loss that is gaseous [1]"""

    # Grazing fractions
    frc_nLUgrs2P: pint.Quantity = field(
        default=Q(0.6, "1"),
        validator=check_units("1"),
    )
    """Fraction of nitrogen uptake flux from landuse gross deforestation which goes to
    the plant pool [1]"""

    frc_nLUgrs2L: pint.Quantity = field(
        default=Q(0.3, "1"),
        validator=check_units("1"),
    )
    """Fraction of nitrogen uptake flux from landuse gross deforestation which goes to the litter pool [1]"""

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
        Calculate the rate of change of nitrogen pools.

        Parameters
        ----------
        NflxPU : pint.Quantity
            Plant uptake flux
        NflxBNF : pint.Quantity
            Biological nitrogen fixation flux
        NflxLP : pint.Quantity
            Litter production flux
        NflxLD : pint.Quantity
            Litter decomposition flux
        NflxSR : pint.Quantity
            Soil respiration flux
        NflxLS : pint.Quantity
            Litter to soil flux
        NflxAD : pint.Quantity
            Atmospheric deposition flux
        NflxFT : pint.Quantity
            Fertilizer transfer flux
        NflxLUgrs : pint.Quantity
            Nitrogen land use flux resulting from gross deforestation
        NflxLUmin : pint.Quantity
            Nitrogen land use flux directly contributing to mineral nitrogen loss

        Returns
        -------
        List[pint.Quantity]
            List containing changes in nitrogen pools [dNplsP_dt, dNplsL_dt, dNplsS_dt, dNplsM_dt]
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
        Calculate the rate of change of nitrogen pools.

        Parameters
        ----------
        NflxPU : pint.Quantity
            Plant uptake flux
        NflxBNF : pint.Quantity
            Biological nitrogen fixation flux
        NflxLP : pint.Quantity
            Litter production flux
        NflxLD : pint.Quantity
            Litter decomposition flux
        NflxSR : pint.Quantity
            Soil respiration flux
        NflxLS : pint.Quantity
            Litter to soil flux
        NflxAD : pint.Quantity
            Atmospheric deposition flux
        NflxFT : pint.Quantity
            Fertilizer transfer flux
        NflxLUgrs : pint.Quantity
            Nitrogen land use flux resulting from gross deforestation
        NflxLUmin : pint.Quantity
            Nitrogen land use flux directly contributing to mineral nitrogen loss

        Returns
        -------
        List[pint.Quantity]
            List containing changes in nitrogen pools [dNplsP_dt, dNplsL_dt, dNplsS_dt, dNplsM_dt]
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
        Calculate all nitrogen fluxes in the system.

        Returns
        -------
        Dict[str, pint.Quantity]
            Dictionary containing all calculated nitrogen fluxes
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
    """Results of running the nitrogen cycle model.

    Attributes:
        ts_collection: Collection of time series results from the run
        calc_n_turnover: Calculator used for nitrogen turnover
        calc_n_cycle: Calculator used for nitrogen cycle
        time_axis: Time points for simulation [yr]
        NflxPU_t: Time-interpolated plant uptake flux [GtN/yr]
        NflxBNF_t: Time-interpolated biological nitrogen fixation flux [GtN/yr]
        NflxAD_t: Time-interpolated atmospheric deposition flux [GtN/yr]
        NflxFT_t: Time-interpolated fertilizer flux [GtN/yr]
        NflxLUgrs_t: Time-interpolated land use gross flux [GtN/yr]
        NflxLUmin_t: Time-interpolated land use mineralization flux [GtN/yr]
        eff_dT_nLP_t: Time-interpolated temperature effect on litter production [1]
        eff_dT_nLD_t: Time-interpolated temperature effect on litter decomposition [1]
        eff_dT_nSR_t: Time-interpolated temperature effect on soil respiration [1]
        eff_dT_nLSgas_t: Time-interpolated temperature effect on gaseous losses [1]
        eff_N_nLP_t: Time-interpolated CN effect on litter production [1]
        eff_N_nLD_t: Time-interpolated CN effect on litter decomposition [1]
        eff_N_nSR_t: Time-interpolated CN effect on soil respiration [1]
    """

    n_state: Dict[str, pint.Quantity]
    calc_n_turnover: NitrogenTurnoverCalculator
    calc_n_cycle: NitrogenCycleCalculator
    time_axis: pint.Quantity

    # Time-interpolated flux functions
    NflxPU_t: Callable[[pint.Quantity], pint.Quantity]
    NflxBNF_t: Callable[[pint.Quantity], pint.Quantity]
    NflxAD_t: Callable[[pint.Quantity], pint.Quantity]
    NflxFT_t: Callable[[pint.Quantity], pint.Quantity]
    NflxLUgrs_t: Callable[[pint.Quantity], pint.Quantity]
    NflxLUmin_t: Callable[[pint.Quantity], pint.Quantity]

    # Time-interpolated effect functions
    eff_dT_nLP_t: Callable[[pint.Quantity], pint.Quantity]
    eff_dT_nLD_t: Callable[[pint.Quantity], pint.Quantity]
    eff_dT_nSR_t: Callable[[pint.Quantity], pint.Quantity]
    eff_dT_nLSgas_t: Callable[[pint.Quantity], pint.Quantity]
    eff_N_nLP_t: Callable[[pint.Quantity], pint.Quantity]
    eff_N_nLD_t: Callable[[pint.Quantity], pint.Quantity]
    eff_N_nSR_t: Callable[[pint.Quantity], pint.Quantity]
    eff_LU_nLP_t: Callable[[pint.Quantity], pint.Quantity]
    eff_LU_nLD_t: Callable[[pint.Quantity], pint.Quantity]
    eff_LU_nSR_t: Callable[[pint.Quantity], pint.Quantity]

    def add_non_state_variables(self) -> xr.Dataset:
        """Add non state variables to the results.

        Returns
        _______
            xr.Dataset: Results including derived variables
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
    """MAGICC's terrestrial nitrogen cycle model.

    This model simulates the dynamics of nitrogen pools in the terrestrial biosphere.
    Time-continuous functions (suffixed with _t) are created through linear interpolation
    of discrete time points and provide continuous values over the simulation period.

    Attributes:
        calc_n_turnover: Calculator for nitrogen turnover processes
        calc_n_cycle: Calculator for nitrogen cycle processes
        NplsP0: Initial plant nitrogen pool size [GtN]
        NplsL0: Initial litter nitrogen pool size [GtN]
        NplsS0: Initial soil nitrogen pool size [GtN]
        NplsM0: Initial mineral nitrogen pool size [GtN]
        time0: Initial time point [yr]
        switch_Npls: Switches for enabling/disabling nitrogen pools [1]
    """

    calc_n_turnover: NitrogenTurnoverCalculator
    calc_n_cycle: NitrogenCycleCalculator

    NplsP0: pint.Quantity = field(
        default=Q(1, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial plant nitrogen pool size [GtN]"""

    NplsL0: pint.Quantity = field(
        default=Q(1, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial litter nitrogen pool size [GtN]"""

    NplsS0: pint.Quantity = field(
        default=Q(1, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial soil nitrogen pool size [GtN]"""

    NplsM0: pint.Quantity = field(
        default=Q(1, "GtN"),
        validator=check_units("GtN"),
    )
    """Initial mineral nitrogen pool size [GtN]"""

    time0: pint.Quantity = field(
        default=Q(1850, "yr"),
        validator=check_units("yr"),
    )
    """Time at which the initialisation values apply [yr]"""

    switch_Npls: list[int] = field(
        default=[1, 1, 1, 1],
    )
    """Nitrogen pool switches (1=enabled, 0=disabled) for [plant, litter, soil, mineral]"""

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

        Parameters
        ----------
        time_axis : pint.Quantity
            Time points for simulation [yr].
        NflxPU_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated plant uptake flux [GtN/yr].
        NflxBNF_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated biological nitrogen fixation flux [GtN/yr].
        NflxAD_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated atmospheric nitrogen deposition flux [GtN/yr].
        NflxFT_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated nitrogen fertilizer flux [GtN/yr].
        NflxLUgrs_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use gross flux [GtN/yr].
        NflxLUmin_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use mineralization flux [GtN/yr].
        eff_dT_nLP_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated temperature effect on litter production [1].
        eff_dT_nLD_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated temperature effect on litter decomposition [1].
        eff_dT_nSR_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated temperature effect on soil respiration [1].
        eff_dT_nLSgas_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated temperature effect on gaseous nitrogen losses [1].
        eff_N_nLP_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated CN effect on litter production [1].
        eff_N_nLD_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated CN effect on litter decomposition [1].
        eff_N_nSR_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated CN effect on soil respiration [1].
        eff_LU_nLP_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use effect on litter production [1].
        eff_LU_nLD_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use effect on litter decomposition [1].
        eff_LU_nSR_t : Callable[[pint.Quantity], pint.Quantity]
            Time-interpolated land use effect on soil respiration [1].

        Returns
        -------
        NitrogenCycleModelResult
            Object containing simulation results.

        Raises
        ------
        ValueError
            If `time_axis` starts before `time0`.
        SolveError
            If the numerical solution fails to converge.

        Notes
        -----
        All *_t arguments are time-continuous functions created through linear interpolation
        of discrete time points. These functions map from time [yr] to their respective units.
        """

        def func_to_solve_NplsPLS(t: float, y: np.ndarray) -> list:
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