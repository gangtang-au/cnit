"""
MAGICC Coupled Carbon-Nitrogen Cycle Model (CNit) Implementation
"""
from dataclasses import asdict, fields
from typing import Dict, List
from attrs import define, field
import numpy as np
from scipy.interpolate import interp1d
import xarray as xr
import pint

from .landuse import LanduseCalculator
from .effects import (
    EffectCO2Calculator,
    EffectTemperatureCalculator,
    EffectCarbonNitrogenCouplingCalculator,
    EffectLanduseCalculator
)
from .carbon_cycle import (
    CarbonNPPLPRCalculator,
    CarbonTurnoverCalculator,
    CarbonCycleCalculator,
    CarbonCycleModel
)
from .nitrogen_cycle import (
    NitrogenPUBNFCalculator,
    NitrogenTurnoverCalculator,
    NitrogenCycleCalculator,
    NitrogenCycleModel
)

from ..configs.model_config import CarbonNitrogenCycleModelConfig
from ..configs.experiment_config import CarbonNitrogenCycleExperimentConfig
from ..utils.units import Q
from ..utils.maths import _get_interp
from ..utils.data_utils import make_dataset_from_var_dict

one = Q(1, "1")
tau_fallback = Q(99.99, "yr")  # if the computed turnover time is 0, NaN, or invalid


@define
class CarbonNitrogenCycleModel:
    """MAGICC's terrestrial carbon nitrogen cycle model.

    This model combines carbon and nitrogen cycle processes to simulate
    terrestrial biogeochemistry dynamics.

    Attributes:
        calc_lu: Calculator for land use processes
        calc_eff_CO2: Calculator for CO2 effects
        calc_eff_dT: Calculator for temperature effects
        calc_eff_CN: Calculator for carbon-nitrogen coupling effects
        calc_eff_LU: Calculator for land use effects
        calc_npp_lpr: Calculator for NPP and litter production
        mdl_c_cycle: Carbon cycle model component
        calc_pu_bnf: Calculator for plant uptake and BNF
        mdl_n_cycle: Nitrogen cycle model component
        switch_N: Switch to enable/disable nitrogen cycle [dimensionless]
    """

    # Calculators
    calc_lu: LanduseCalculator
    calc_eff_CO2: EffectCO2Calculator
    calc_eff_dT: EffectTemperatureCalculator
    calc_eff_CN: EffectCarbonNitrogenCouplingCalculator
    calc_eff_LU: EffectLanduseCalculator
    calc_npp_lpr: CarbonNPPLPRCalculator
    calc_pu_bnf: NitrogenPUBNFCalculator

    # Model components
    mdl_c_cycle: CarbonCycleModel
    mdl_n_cycle: NitrogenCycleModel

    # Model switches
    switch_N: pint.Quantity = field(
        default=Q(1, "1"),
        # validator=check_units("1")
    )

    @classmethod
    def from_dict(cls, params: dict) -> "CarbonNitrogenCycleModel":
        """
        Build a fully configured CarbonNitrogenCycleModel from a parameter dictionary.

        Parameters
        ----------
        params : dict
            Parameter dictionary for configuring calculators and sub-models.

        Returns
        -------
        CarbonNitrogenCycleModel
            Configured model instance.
        """

        # Optionally: check params keys against config fields for typo catching
        valid_keys = {f.name for f in fields(CarbonNitrogenCycleModelConfig)}
        invalid_keys = set(params) - valid_keys
        if invalid_keys:
            raise ValueError(f"Unknown parameter(s) in config: {invalid_keys}")

        # Instantiation (explicit and clear)
        calc_lu = LanduseCalculator()
        calc_eff_CO2 = EffectCO2Calculator()
        calc_eff_dT = EffectTemperatureCalculator()
        calc_eff_CN = EffectCarbonNitrogenCouplingCalculator()
        calc_eff_LU = EffectLanduseCalculator()
        calc_npp_lpr = CarbonNPPLPRCalculator()
        calc_pu_bnf = NitrogenPUBNFCalculator()
        calc_c_turnover = CarbonTurnoverCalculator()
        calc_c_cycle = CarbonCycleCalculator()
        mdl_c_cycle = CarbonCycleModel(
            calc_c_turnover=calc_c_turnover,
            calc_c_cycle=calc_c_cycle,
        )
        calc_n_turnover = NitrogenTurnoverCalculator()
        calc_n_cycle = NitrogenCycleCalculator()
        mdl_n_cycle = NitrogenCycleModel(
            calc_n_turnover=calc_n_turnover,
            calc_n_cycle=calc_n_cycle,
        )
        mdl_cn_cycle = cls(
            calc_lu=calc_lu,
            calc_eff_CO2=calc_eff_CO2,
            calc_eff_dT=calc_eff_dT,
            calc_eff_CN=calc_eff_CN,
            calc_eff_LU=calc_eff_LU,
            calc_npp_lpr=calc_npp_lpr,
            mdl_c_cycle=mdl_c_cycle,
            calc_pu_bnf=calc_pu_bnf,
            mdl_n_cycle=mdl_n_cycle,
        )
        # Set parameters on all components
        for obj in [
            mdl_cn_cycle,
            calc_lu, calc_eff_CO2, calc_eff_dT, calc_eff_CN, calc_eff_LU,
            calc_npp_lpr, calc_pu_bnf, calc_c_turnover, calc_c_cycle,
            calc_n_turnover, calc_n_cycle, mdl_c_cycle, mdl_n_cycle
        ]:
            for attr, value in params.items():
                if hasattr(obj, attr):
                    setattr(obj, attr, value)

        return mdl_cn_cycle

    @classmethod
    def from_config(
            cls,
            config: CarbonNitrogenCycleModelConfig
    ) -> "CarbonNitrogenCycleModel":
        """
        Build a fully configured CarbonNitrogenCycleModel from a config dataclass.

        Parameters
        ----------
        config : CarbonNitrogenCycleModelConfig
            Configuration dataclass for calculators and sub-models.

        Returns
        -------
        CarbonNitrogenCycleModel
            Configured model instance.
        """
        cfg_dict = asdict(config)
        return cls.from_dict(cfg_dict)

    def run_experiments(self, experiments: List["CarbonNitrogenCycleExperimentConfig"]):
        """
        Run multiple experiments with the model.

        Parameters
        ----------
        experiments : List[CarbonNitrogenCycleExperimentConfig]
            List of experiment configurations.

        Returns
        -------
        xr.Dataset
            dataset containing results for all experiments.
        """

        input_dicts = [exp.to_input_dict() for exp in experiments]

        results = []
        for exp, input_dict in zip(experiments, input_dicts):
            self.switch_N = exp.switch_N
            res_exp = self.run(**input_dict)
            res_exp = res_exp.assign_coords(experiment=exp.name).expand_dims(
                "experiment")
            results.append(res_exp)
        return xr.concat(results, dim="experiment", join="outer")

    def run(
            self,
            time_axis: pint.Quantity,
            dT_s: pint.Quantity,
            CO2_s: pint.Quantity,
            CemsLUnet_s: pint.Quantity,
            CemsLUgrs_s: pint.Quantity,
            NflxAD_s: pint.Quantity,
            NflxFT_s: pint.Quantity,
            NemsLUnet_s: pint.Quantity,
            NemsLUgrs_s: pint.Quantity,
            NemsLUmin_s: pint.Quantity,
    ):
        r"""
        Run the carbon-nitrogen cycle model

        Parameters
        ----------
        time_axis : pint.Quantity
            The time points for the simulation [yr].
        dT_s : pint.Quantity
            Time series of temperature anomaly [K].
        CO2_s : pint.Quantity
            Time series of atmospheric CO2 concentration [ppm].
        CemsLUnet_s : pint.Quantity
            Time series of net land-use carbon emissions [GtC/yr].
        CemsLUgrs_s : pint.Quantity
            Time series of gross land-use carbon emissions [GtC/yr].
        NflxAD_s : pint.Quantity
            Time series of atmospheric nitrogen deposition flux [GtN/yr].
        NflxFT_s : pint.Quantity
            Time series of nitrogen fertilizer flux [GtN/yr].
        NemsLUnet_s : pint.Quantity
            Time series of net land-use nitrogen emissions [GtN/yr].
        NemsLUgrs_s : pint.Quantity
            Time series of gross land-use nitrogen emissions [GtN/yr].
        NemsLUmin_s : pint.Quantity
            Time series of land-use nitrogen emissions from the mineral pool [GtN/yr].

        Returns
        -------
        xarray.Dataset
            A dataset containing all carbon-nitrogen fluxes and state variables (e.g.,
            CplsP, CplsL, CplsS, NplsP, NplsL, NplsS, NplsM) for the simulation period.
        """
        CflxLU_s = self.calc_lu.calculate_CflxLU_series(
            CemsLUnet_s=CemsLUnet_s,
            CemsLUgrs_s=CemsLUgrs_s,
        )
        CflxLUrgr_s = CflxLU_s["CflxLUrgr"]
        CflxLUgrs_s = CflxLU_s["CflxLUgrs"]

        CflxNPPLPR_NflxPUBNF_s = self.calculate_CflxNPPLPR_NflxPUBNF(
            dT=dT_s,
            CO2=CO2_s,
            CflxLUrgr=CflxLUrgr_s,
            cumsum_CflxLUrgr=CflxLU_s["cumsum_CflxLUrgr"],
            cumsum_CflxLUdfst=CflxLU_s["cumsum_CflxLUdfst"],
            cumsum_CflxLUafst=CflxLU_s["cumsum_CflxLUafst"],
            cumsum_CflxLUafst_decay=CflxLU_s["cumsum_CflxLUafst_decay"],
            NflxAD=NflxAD_s,
            NflxFT=NflxFT_s,
        )
        CflxNPP_s = CflxNPPLPR_NflxPUBNF_s["CflxNPP"]
        CflxLPR_s = CflxNPPLPR_NflxPUBNF_s["CflxLPR"]
        NflxPUdef_s = CflxNPPLPR_NflxPUBNF_s["NflxPUdef"]
        NflxPU_s = CflxNPPLPR_NflxPUBNF_s["NflxPU"]
        NflxBNF_s = CflxNPPLPR_NflxPUBNF_s["NflxBNF"]

        res_dict_extra = {
            k: v for k, v in (CflxLU_s | CflxNPPLPR_NflxPUBNF_s).items()
            if k not in ["CflxLUgrs", "CflxNPP", "CflxLPR", "NflxPU", "NflxBNF"]
        }
        res_extra = make_dataset_from_var_dict(res_dict_extra, time_axis)

        ones = np.ones_like(time_axis)
        if sum(self.mdl_c_cycle.switch_Cpls) != 0:
            eff_dT_cLP, eff_dT_cLD, eff_dT_cSR = (
                self.calc_eff_dT.calculate_eff_dT_cLPLDSR(dT=dT_s)
            )
            eff_N_cLP, eff_N_cLD, eff_N_cSR = self.calc_eff_CN.calculate_eff_N_cLPLDSR(
                NflxPUdef=NflxPUdef_s,
                switch_N=self.switch_N,
            )
            eff_LU_cLP, eff_LU_cLD, eff_LU_cSR = (
                self.calc_eff_LU.calculate_eff_LU_cLPLDSR(CflxLUrgr=CflxLUrgr_s)
            )

            res_c = (
                self.mdl_c_cycle.run(
                    time_axis=time_axis,
                    CflxNPP_t=_get_interp(CflxNPP_s, time_axis),
                    CflxLPR_t=_get_interp(CflxLPR_s, time_axis),
                    CflxLUgrs_t=_get_interp(CflxLUgrs_s, time_axis),
                    eff_dT_cLP_t=interp1d(time_axis, eff_dT_cLP),
                    eff_dT_cLD_t=interp1d(time_axis, eff_dT_cLD),
                    eff_dT_cSR_t=interp1d(time_axis, eff_dT_cSR),
                    eff_N_cLP_t=interp1d(time_axis, eff_N_cLP),
                    eff_N_cLD_t=interp1d(time_axis, eff_N_cLD),
                    eff_N_cSR_t=interp1d(time_axis, eff_N_cSR),
                    eff_LU_cLP_t=interp1d(time_axis, eff_LU_cLP),
                    eff_LU_cLD_t=interp1d(time_axis, eff_LU_cLD),
                    eff_LU_cSR_t=interp1d(time_axis, eff_LU_cSR),
                )
                .add_non_state_variables()
            )
        else:
            res_dict_c = {
                "CplsP": self.mdl_c_cycle.CplsP0 * ones,
                "CplsL": self.mdl_c_cycle.CplsL0 * ones,
                "CplsS": self.mdl_c_cycle.CplsS0 * ones,
            }
            res_c = make_dataset_from_var_dict(res_dict_c, time_axis)

        if sum(self.mdl_n_cycle.switch_Npls) != 0 and self.switch_N != 0:
            NflxLU_s = self.calc_lu.calculate_NflxLU_series(
                NemsLUnet_s=NemsLUnet_s,
                NemsLUgrs_s=NemsLUgrs_s,
                NemsLUmin_s=NemsLUmin_s,
            )
            NflxLUgrs_s = NflxLU_s["NflxLUgrs"]
            NflxLUmin_s = NflxLU_s["NflxLUmin"]

            eff_dT_nLP, eff_dT_nLD, eff_dT_nSR, eff_dT_nLSgas = (
                self.calc_eff_dT.calculate_eff_dT_nLPLDSRLS(dT=dT_s)
            )
            eff_N_nLP, eff_N_nLD, eff_N_nSR = self.calc_eff_CN.calculate_eff_N_nLPLDSR(
                NflxPUdef=NflxPUdef_s,
            )
            eff_LU_nLP, eff_LU_nLD, eff_LU_nSR = (
                self.calc_eff_LU.calculate_eff_LU_nLPLDSR(CflxLUrgr=CflxLUrgr_s)
            )
            res_n = (
                self.mdl_n_cycle.run(
                    time_axis=time_axis,
                    NflxPU_t=_get_interp(NflxPU_s, time_axis),
                    NflxBNF_t=_get_interp(NflxBNF_s, time_axis),
                    NflxAD_t=_get_interp(NflxAD_s, time_axis),
                    NflxFT_t=_get_interp(NflxFT_s, time_axis),
                    NflxLUgrs_t=_get_interp(NflxLUgrs_s, time_axis),
                    NflxLUmin_t=_get_interp(NflxLUmin_s, time_axis),
                    eff_dT_nLP_t=interp1d(time_axis, eff_dT_nLP),
                    eff_dT_nLD_t=interp1d(time_axis, eff_dT_nLD),
                    eff_dT_nSR_t=interp1d(time_axis, eff_dT_nSR),
                    eff_dT_nLSgas_t=interp1d(time_axis, eff_dT_nLSgas),
                    eff_N_nLP_t=interp1d(time_axis, eff_N_nLP),
                    eff_N_nLD_t=interp1d(time_axis, eff_N_nLD),
                    eff_N_nSR_t=interp1d(time_axis, eff_N_nSR),
                    eff_LU_nLP_t=interp1d(time_axis, eff_LU_nLP),
                    eff_LU_nLD_t=interp1d(time_axis, eff_LU_nLD),
                    eff_LU_nSR_t=interp1d(time_axis, eff_LU_nSR),
                )
                .add_non_state_variables()
            )
        else:
            res_dict_n = {
                "NplsP": self.mdl_n_cycle.NplsP0 * ones,
                "NplsL": self.mdl_n_cycle.NplsL0 * ones,
                "NplsS": self.mdl_n_cycle.NplsS0 * ones,
                "NplsM": self.mdl_n_cycle.NplsM0 * ones,
            }
            res_n = make_dataset_from_var_dict(res_dict_n, time_axis)

        return xr.merge([res_extra, res_c, res_n])

    def calculate_CflxNPPLPR_NflxPUBNF(
            self,
            dT: pint.Quantity,
            CO2: pint.Quantity,
            CflxLUrgr: pint.Quantity,
            cumsum_CflxLUrgr: pint.Quantity,
            cumsum_CflxLUdfst: pint.Quantity,
            cumsum_CflxLUafst: pint.Quantity,
            cumsum_CflxLUafst_decay: pint.Quantity,
            NflxAD: pint.Quantity,
            NflxFT: pint.Quantity,
    ):
        """
        Calculates key C and N fluxes based on environmental drivers.

        This method first calculates the "potential" NPP and LPR (assuming
        no nitrogen limitation). It then calculates the N-limitation effect
        (NflxPUdef) and re-calculates the "actual" (N-limited) NPP and LPR.
        Finally, it calculates the resulting nitrogen uptake and fixation fluxes.

        Parameters
        ----------
        dT : pint.Quantity
            Time series of temperature anomaly [K].
        CO2 : pint.Quantity
            Time series of atmospheric CO2 concentration [ppm].
        CflxLUrgr : pint.Quantity
            Time series of land-use regrowth carbon flux [GtC/yr].
        cumsum_CflxLUrgr : pint.Quantity
            Time series of cumulative land-use regrowth [GtC].
        cumsum_CflxLUdfst : pint.Quantity
            Time series of cumulative land-use deforestation [GtC].
        cumsum_CflxLUafst : pint.Quantity
            Time series of cumulative land-use afforestation [GtC].
        cumsum_CflxLUafst_decay : pint.Quantity
            Time series of cumulative land-use afforestation decay [GtC].
        NflxAD : pint.Quantity
            Time series of atmospheric nitrogen deposition flux [GtN/yr].
        NflxFT : pint.Quantity
            Time series of nitrogen fertilizer flux [GtN/yr].

        Returns
        -------
        Dict[str, pint.Quantity]
            A dictionary containing the following time series:
            - CflxNPP: Actual (N-limited) Net Primary Production [GtC/yr].
            - CflxLPR: Actual (N-limited) Litter Production Respiration [GtC/yr].
            - CflxNPPpot: Potential (non-N-limited) NPP [GtC/yr].
            - CflxLPRpot: Potential (non-N-limited) LPR [GtC/yr].
            - NflxPU: Plant Uptake nitrogen flux [GtN/yr].
            - NflxBNF: Biological Nitrogen Fixation flux [GtN/yr].
            - NflxPUreq: Required N for plant uptake [GtN/yr].
            - NflxPUavail: Available N for plant uptake [GtN/yr].
            - NflxPUavail_netMIN: N from net mineralization [GtN/yr].
            - NflxPUdef: N deficit for plant uptake [GtN/yr].
        """

        # get the potential CflxNPP
        eff_CO2_NPP = self.calc_eff_CO2.calculate_eff_CO2_NPP(CO2=CO2)
        eff_dT_NPP, eff_dT_LPR = self.calc_eff_dT.calculate_eff_dT_NPPLPR(dT=dT)
        eff_LU_NPP = self.calc_eff_LU.calculate_eff_LU_NPP(
            cumsum_CflxLUrgr=cumsum_CflxLUrgr,
            cumsum_CflxLUdfst=cumsum_CflxLUdfst,
            cumsum_CflxLUafst=cumsum_CflxLUafst,
            cumsum_CflxLUafst_decay=cumsum_CflxLUafst_decay,
            CplsPLS0=self.mdl_c_cycle.CplsP0
                     + self.mdl_c_cycle.CplsL0
                     + self.mdl_c_cycle.CplsS0,
            frc_rgr=self.calc_lu.frc_rgr,
            frc_afst_decay=1 - self.calc_lu.frc_rgr,
        )
        CflxNPPLPRpot = self.calc_npp_lpr.calculate(
            eff_CO2_NPP=eff_CO2_NPP,
            eff_dT_NPP=eff_dT_NPP,
            eff_dT_LPR=eff_dT_LPR,
            eff_N_NPP=one,
            eff_LU_NPP=eff_LU_NPP,
            CflxLUrgr=CflxLUrgr,
        )
        CflxNPPpot = CflxNPPLPRpot["CflxNPP"]
        CflxLPRpot = CflxNPPLPRpot["CflxLPR"]

        # check if there is a nitrogen cycle and it is on
        if self.mdl_n_cycle.NplsP0 != 0 and self.switch_N != 0:
            NflxPU_req_avail_def = self.calc_pu_bnf.calculate_NflxPU_req_avail_def(
                CflxNPP=CflxNPPpot,
                CflxNPP0=self.calc_npp_lpr.NPP0,
                NflxAD=NflxAD,
                NflxFT=NflxFT,
                dT=dT,
            )
            NflxPUdef = NflxPU_req_avail_def["NflxPUdef"]
            eff_N_NPP = self.calc_eff_CN.calculate_eff_N_NPP(
                NflxPUdef=NflxPUdef,
                switch_N=self.switch_N,
            )
            # N effect on CflxNPP and CflxLPR based on the N requirement
            CflxNPPLPR = self.calc_npp_lpr.calculate(
                CflxLUrgr=CflxLUrgr,
                eff_CO2_NPP=eff_CO2_NPP,
                eff_dT_NPP=eff_dT_NPP,
                eff_N_NPP=eff_N_NPP,
                eff_LU_NPP=eff_LU_NPP,
                eff_dT_LPR=eff_dT_LPR,
            )
            CflxNPP = CflxNPPLPR["CflxNPP"]
            CflxLPR = CflxNPPLPR["CflxLPR"]

            # use the N limited CflxNPP to get real NflxPU
            eff_dT_PU, eff_dT_BNF = self.calc_eff_dT.calculate_eff_dT_PUBNF(dT=dT)
            eff_C_PU, eff_C_BNF = self.calc_eff_CN.calculate_eff_C_PUBNF(
                CflxNPP=CflxNPP,
                CflxNPP0=self.calc_npp_lpr.NPP0,
            )
            eff_N_PU, eff_N_BNF = self.calc_eff_CN.calculate_eff_N_PUBNF(
                NflxPUdef=NflxPUdef
            )

            NflxPUBNF_s = self.calc_pu_bnf.calculate(
                eff_dT_PU=eff_dT_PU,
                eff_C_PU=eff_C_PU,
                eff_N_PU=eff_N_PU,
                eff_dT_BNF=eff_dT_BNF,
                eff_C_BNF=eff_C_BNF,
                eff_N_BNF=eff_N_BNF,
            )
            NflxPUBNF_s = NflxPUBNF_s | NflxPU_req_avail_def

        else:
            CflxNPP = CflxNPPpot
            CflxLPR = CflxLPRpot
            NflxPUBNF_s = {
                key: NflxAD * np.nan
                for key in [
                    "NflxPU",
                    "NflxBNF",
                    "NflxPUreq",
                    "NflxPUavail",
                    "NflxPUavail_netMIN",
                    "NflxPUdef",
                ]
            }

        CflxNPPLPR_s = {
            "CflxNPP": CflxNPP,
            "CflxLPR": CflxLPR,
            "CflxNPPpot": CflxNPPpot,
            "CflxLPRpot": CflxLPRpot,
        }

        return CflxNPPLPR_s | NflxPUBNF_s

    def calculate_mdl_eqm0(
            self,
            dT_s: pint.Quantity,
            CO2_s: pint.Quantity,
            CemsLUnet_s: pint.Quantity,
            CemsLUgrs_s: pint.Quantity,
            NflxAD_s: pint.Quantity,
            NflxFT_s: pint.Quantity,
            NemsLUnet_s: pint.Quantity,
            NemsLUgrs_s: pint.Quantity,
            NemsLUmin_s: pint.Quantity,
    ):
        CflxLU_s = self.calc_lu.calculate_CflxLU_series(
            CemsLUnet_s=CemsLUnet_s[:1],
            CemsLUgrs_s=CemsLUgrs_s[:1],
        )
        CflxLUrgr_s = CflxLU_s["CflxLUrgr"]

        CflxNPPLPR_NflxPUBNF = self.calculate_CflxNPPLPR_NflxPUBNF(
            dT=dT_s[0],
            CO2=CO2_s[0],
            CflxLUrgr=CflxLUrgr_s[0],
            cumsum_CflxLUrgr=CflxLU_s["cumsum_CflxLUrgr"][0],
            cumsum_CflxLUdfst=CflxLU_s["cumsum_CflxLUdfst"][0],
            cumsum_CflxLUafst=CflxLU_s["cumsum_CflxLUafst"][0],
            cumsum_CflxLUafst_decay=CflxLU_s["cumsum_CflxLUafst_decay"][0],
            NflxAD=NflxAD_s[0],
            NflxFT=NflxFT_s[0],
        )

        NflxLU_s = self.calc_lu.calculate_NflxLU_series(
            NemsLUnet_s=NemsLUnet_s[:1],
            NemsLUgrs_s=NemsLUgrs_s[:1],
            NemsLUmin_s=NemsLUmin_s[:1],
        )

        return self.calculate_eqm0(
            dT=dT_s[0],
            CflxNPP=CflxNPPLPR_NflxPUBNF["CflxNPP"],
            CflxLPR=CflxNPPLPR_NflxPUBNF["CflxLPR"],
            CflxLUgrs=CflxLU_s["CflxLUgrs"][0],
            CflxLUrgr=CflxLUrgr_s[0],
            NflxPUdef=CflxNPPLPR_NflxPUBNF["NflxPUdef"],
            NflxPU=CflxNPPLPR_NflxPUBNF["NflxPU"],
            NflxBNF=CflxNPPLPR_NflxPUBNF["NflxBNF"],
            NflxAD=NflxAD_s[0],
            NflxFT=NflxFT_s[0],
            NflxLUgrs=NflxLU_s["NflxLUgrs"][0],
            NflxLUmin=NflxLU_s["NflxLUmin"][0],
        )

    def calculate_eqm0(
            self,
            dT: pint.Quantity,
            CflxNPP: pint.Quantity,
            CflxLPR: pint.Quantity,
            CflxLUgrs: pint.Quantity,
            CflxLUrgr: pint.Quantity,
            NflxPUdef: pint.Quantity,
            NflxPU: pint.Quantity,
            NflxBNF: pint.Quantity,
            NflxAD: pint.Quantity,
            NflxFT: pint.Quantity,
            NflxLUgrs: pint.Quantity,
            NflxLUmin: pint.Quantity,
    ):
        eff_dT_cLP, eff_dT_cLD, eff_dT_cSR = self.calc_eff_dT.calculate_eff_dT_cLPLDSR(
            dT=dT
        )
        eff_N_cLP, eff_N_cLD, eff_N_cSR = self.calc_eff_CN.calculate_eff_N_cLPLDSR(
            NflxPUdef=NflxPUdef,
            switch_N=self.switch_N,
        )
        eff_LU_cLP, eff_LU_cLD, eff_LU_cSR = self.calc_eff_LU.calculate_eff_LU_cLPLDSR(
            CflxLUrgr=CflxLUrgr,
        )

        carbon_cycle_eqm = self.inverse_carbon_turnover(
            CflxNPP=CflxNPP,
            CflxLPR=CflxLPR,
            CflxLUgrs=CflxLUgrs,
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

        # check if there is a nitrogen cycle, and it is on
        if self.mdl_n_cycle.NplsP0 != 0 and self.switch_N != 0:
            eff_dT_nLP, eff_dT_nLD, eff_dT_nSR, eff_dT_nLSgas = (
                self.calc_eff_dT.calculate_eff_dT_nLPLDSRLS(dT=dT)
            )
            eff_N_nLP, eff_N_nLD, eff_N_nSR = self.calc_eff_CN.calculate_eff_N_nLPLDSR(
                NflxPUdef=NflxPUdef,
            )
            eff_LU_nLP, eff_LU_nLD, eff_LU_nSR = (
                self.calc_eff_LU.calculate_eff_LU_nLPLDSR(
                    CflxLUrgr=CflxLUrgr,
                )
            )
            nitrogen_cycle_eqm = self.inverse_nitrogen_turnover(
                NflxPU=NflxPU,
                NflxBNF=NflxBNF,
                NflxAD=NflxAD,
                NflxFT=NflxFT,
                NflxLUgrs=NflxLUgrs,
                NflxLUmin=NflxLUmin,
                eff_dT_nLP=eff_dT_nLP,
                eff_dT_nLD=eff_dT_nLD,
                eff_dT_nSR=eff_dT_nSR,
                eff_dT_nLSgas=eff_dT_nLSgas,
                eff_N_nLP=eff_N_nLP,
                eff_N_nLD=eff_N_nLD,
                eff_N_nSR=eff_N_nSR,
                eff_LU_nLP=eff_LU_nLP,
                eff_LU_nLD=eff_LU_nLD,
                eff_LU_nSR=eff_LU_nSR,
            )
            return carbon_cycle_eqm | nitrogen_cycle_eqm

        else:
            return carbon_cycle_eqm

    def inverse_turnover_flux(
            self,
            influx: pint.Quantity,
            outflux: pint.Quantity,
            pool_size_change: pint.Quantity,
    ):
        return influx - outflux - pool_size_change

    def inverse_carbon_turnover(
            self,
            CflxNPP: pint.Quantity,
            CflxLPR: pint.Quantity,
            CflxLUgrs: pint.Quantity,
            eff_dT_cLP: pint.Quantity,
            eff_dT_cLD: pint.Quantity,
            eff_dT_cSR: pint.Quantity,
            eff_N_cLP: pint.Quantity,
            eff_N_cLD: pint.Quantity,
            eff_N_cSR: pint.Quantity,
            eff_LU_cLP: pint.Quantity,
            eff_LU_cLD: pint.Quantity,
            eff_LU_cSR: pint.Quantity,
            CplsP: pint.Quantity = None,
            CplsL: pint.Quantity = None,
            CplsS: pint.Quantity = None,
            dCplsP_dt: pint.Quantity = Q(0, "GtC/yr"),
            dCplsL_dt: pint.Quantity = Q(0, "GtC/yr"),
            dCplsS_dt: pint.Quantity = Q(0, "GtC/yr"),
    ) -> Dict[str, pint.Quantity]:
        """
        Calculate equilibrium values for carbon fluxes and turnover times

        Parameters
        ----------
        CplsP
            Plant carbon pool [GtC]
        CplsL
            Litter carbon pool [GtC]
        CplsS
            Soil carbon pool [GtC]
        eff_dT_cLP
            Temperature effect on carbon litter production [dimensionless]
        eff_dT_cLD
            Temperature effect on carbon litter decomposition [dimensionless]
        eff_dT_cSR
            Temperature effect on carbon soil respiration [dimensionless]
        eff_N_cLP
            Carbon-nitrogen coupling effect on carbon litter production [dimensionless]
        eff_N_cLD
            Carbon-nitrogen coupling effect on carbon litter decomposition [dimensionless]
        eff_N_cSR
            Carbon-nitrogen coupling effect on carbon soil respiration [dimensionless]
        eff_LU_cLP
            Land use effect on carbon litter production [dimensionless]
        eff_LU_cLD
            Land use effect on carbon litter decomposition [dimensionless]
        eff_LU_cSR
            Land use effect on carbon soil respiration [dimensionless]
        CflxNPP
            Net Primary Production carbon flux [GtC/yr]
        CflxLPR
            Litter Production Respiration carbon flux [GtC/yr]
        CflxLUgrs
            Gross Land Use carbon flux [GtC/yr]

        Returns
        -------
            Dictionary with equilibrium values for carbon fluxes and turnover times
        """
        calc_c_cycle = self.mdl_c_cycle.calc_c_cycle

        # Calculate NPP fractions to each pool
        CflxNPP2P = CflxNPP * calc_c_cycle.frc_NPP2P
        CflxNPP2L = CflxNPP * calc_c_cycle.frc_NPP2L
        CflxNPP2S = CflxNPP - CflxNPP2P - CflxNPP2L

        # Calculate land use emissions fractions
        CflxLUgrs2P = CflxLUgrs * calc_c_cycle.frc_cLUgrs2P
        CflxLUgrs2L = CflxLUgrs * calc_c_cycle.frc_cLUgrs2L
        CflxLUgrs2S = CflxLUgrs - CflxLUgrs2P - CflxLUgrs2L

        # Calculate equilibrium fluxes
        CflxLP_inverse = self.inverse_turnover_flux(
            influx=CflxNPP2P,
            outflux=CflxLPR + CflxLUgrs2P,
            pool_size_change=dCplsP_dt,
        )
        CflxLP2L_inverse = CflxLP_inverse * calc_c_cycle.frc_cLP2L
        CflxLP2S_inverse = CflxLP_inverse - CflxLP2L_inverse

        CflxLD_inverse = self.inverse_turnover_flux(
            influx=CflxNPP2L + CflxLP2L_inverse,
            outflux=CflxLUgrs2L,
            pool_size_change=dCplsL_dt,
        )
        CflxLD2S_inverse = CflxLD_inverse * calc_c_cycle.frc_cLD2S
        CflxLD2A_inverse = CflxLD_inverse - CflxLD2S_inverse

        CflxSR_inverse = self.inverse_turnover_flux(
            influx=CflxNPP2S + CflxLP2S_inverse + CflxLD2S_inverse,
            outflux=CflxLUgrs2S,
            pool_size_change=dCplsS_dt,
        )

        # Calculate equilibrium turnover times
        CplsP = CplsP if CplsP is not None else self.mdl_c_cycle.CplsP0
        CplsL = CplsL if CplsL is not None else self.mdl_c_cycle.CplsL0
        CplsS = CplsS if CplsS is not None else self.mdl_c_cycle.CplsS0

        tau_CplsP_inverse = CplsP / CflxLP_inverse * eff_dT_cLP * eff_N_cLP * eff_LU_cLP
        tau_CplsL_inverse = CplsL / CflxLD_inverse * eff_dT_cLD * eff_N_cLD * eff_LU_cLD
        tau_CplsS_inverse = CplsS / CflxSR_inverse * eff_dT_cSR * eff_N_cSR * eff_LU_cSR

        # Handle edge case for litter pool
        tau_CplsL_inverse = (
            tau_fallback
            if not tau_CplsL_inverse or np.isnan(tau_CplsL_inverse)
            else tau_CplsL_inverse
        )

        # Return dictionary of equilibrium values
        return {
            "CflxLP_inverse": CflxLP_inverse,
            "CflxLP2L_inverse": CflxLP2L_inverse,
            "CflxLP2S_inverse": CflxLP2S_inverse,
            "CflxLD_inverse": CflxLD_inverse,
            "CflxLD2S_inverse": CflxLD2S_inverse,
            "CflxLD2A_inverse": CflxLD2A_inverse,
            "CflxSR_inverse": CflxSR_inverse,
            "tau_CplsP_inverse": tau_CplsP_inverse,
            "tau_CplsL_inverse": tau_CplsL_inverse,
            "tau_CplsS_inverse": tau_CplsS_inverse,
        }

    def inverse_nitrogen_turnover(
            self,
            NflxPU: pint.Quantity,
            NflxBNF: pint.Quantity,
            NflxAD: pint.Quantity,
            NflxFT: pint.Quantity,
            NflxLUgrs: pint.Quantity,
            NflxLUmin: pint.Quantity,
            eff_dT_nLP: pint.Quantity,
            eff_dT_nLD: pint.Quantity,
            eff_dT_nSR: pint.Quantity,
            eff_dT_nLSgas: pint.Quantity,
            eff_N_nLP: pint.Quantity,
            eff_N_nLD: pint.Quantity,
            eff_N_nSR: pint.Quantity,
            eff_LU_nLP: pint.Quantity,
            eff_LU_nLD: pint.Quantity,
            eff_LU_nSR: pint.Quantity,
            NplsP: pint.Quantity = None,
            NplsL: pint.Quantity = None,
            NplsS: pint.Quantity = None,
            NplsM: pint.Quantity = None,
            dNplsP_dt: pint.Quantity = Q(0, "GtN/yr"),
            dNplsL_dt: pint.Quantity = Q(0, "GtN/yr"),
            dNplsS_dt: pint.Quantity = Q(0, "GtN/yr"),
            dNplsM_dt: pint.Quantity = Q(0, "GtN/yr"),
    ):
        calc_n_cycle = self.mdl_n_cycle.calc_n_cycle

        NflxPU2P = NflxPU * calc_n_cycle.frc_PU2P
        NflxPU2L = NflxPU * calc_n_cycle.frc_PU2L
        NflxPU2S = NflxPU - NflxPU2P - NflxPU2L

        # Calculate BNF partitioning
        NflxBNF2P = NflxBNF * calc_n_cycle.frc_BNF2P
        NflxBNF2L = NflxBNF * calc_n_cycle.frc_BNF2L
        NflxBNF2S = NflxBNF - NflxBNF2P - NflxBNF2L

        # Calculate gross land use emission partitioning
        NflxLUgrs2P = NflxLUgrs * calc_n_cycle.frc_nLUgrs2P
        NflxLUgrs2L = NflxLUgrs * calc_n_cycle.frc_nLUgrs2L
        NflxLUgrs2S = NflxLUgrs - NflxLUgrs2P - NflxLUgrs2L

        # Calculate equilibrium fluxes
        NflxLP_inverse = self.inverse_turnover_flux(
            influx=NflxPU2P + NflxBNF2P,
            outflux=NflxLUgrs2P,
            pool_size_change=dNplsP_dt,
        )
        NflxLP2L_inverse = NflxLP_inverse * calc_n_cycle.frc_nLP2L
        NflxLP2S_inverse = NflxLP_inverse - NflxLP2L_inverse

        NflxLD_inverse = self.inverse_turnover_flux(
            influx=NflxPU2L + NflxBNF2L + NflxLP2L_inverse,
            outflux=NflxLUgrs2L,
            pool_size_change=dNplsL_dt,
        )
        NflxLD2S_inverse = NflxLD_inverse * calc_n_cycle.frc_nLD2S
        NflxLD2M_inverse = NflxLD_inverse - NflxLD2S_inverse

        NflxSR_inverse = self.inverse_turnover_flux(
            influx=NflxPU2S + NflxBNF2S + NflxLP2S_inverse + NflxLD2S_inverse,
            outflux=NflxLUgrs2S,
            pool_size_change=dNplsS_dt,
        )

        NflxLS_inverse = self.inverse_turnover_flux(
            influx=NflxAD + NflxFT + NflxSR_inverse + NflxLD2M_inverse,
            outflux=NflxPU + NflxLUmin,
            pool_size_change=dNplsM_dt,
        )

        # Calculate equilibrium turnover times
        NplsP = NplsP if NplsP is not None else self.mdl_n_cycle.NplsP0
        NplsL = NplsL if NplsL is not None else self.mdl_n_cycle.NplsL0
        NplsS = NplsS if NplsS is not None else self.mdl_n_cycle.NplsS0
        NplsM = NplsM if NplsM is not None else self.mdl_n_cycle.NplsM0

        tau_NplsP_inverse = NplsP / NflxLP_inverse * eff_dT_nLP * eff_N_nLP * eff_LU_nLP
        tau_NplsL_inverse = NplsL / NflxLD_inverse * eff_dT_nLD * eff_N_nLD * eff_LU_nLD
        tau_NplsS_inverse = NplsS / NflxSR_inverse * eff_dT_nSR * eff_N_nSR * eff_LU_nSR

        frc_nLSgas = self.mdl_n_cycle.calc_n_turnover.frc_nLSgas
        tau_NplsM_inverse = (
                NplsM / NflxLS_inverse * (frc_nLSgas * eff_dT_nLSgas + (1 - frc_nLSgas))
        )

        # Handle edge case for litter pool
        tau_NplsL_inverse = (
            tau_fallback
            if not tau_NplsL_inverse or np.isnan(tau_NplsL_inverse)
            else tau_NplsL_inverse
        )

        return {
            "NflxLP_inverse": NflxLP_inverse,
            "NflxLP2L_inverse": NflxLP2L_inverse,
            "NflxLP2S_inverse": NflxLP2S_inverse,
            "NflxLD_inverse": NflxLD_inverse,
            "NflxLD2S_inverse": NflxLD2S_inverse,
            "NflxLD2M_inverse": NflxLD2M_inverse,
            "NflxSR_inverse": NflxSR_inverse,
            "NflxLS_inverse": NflxLS_inverse,
            "tau_NplsP_inverse": tau_NplsP_inverse,
            "tau_NplsL_inverse": tau_NplsL_inverse,
            "tau_NplsS_inverse": tau_NplsS_inverse,
            "tau_NplsM_inverse": tau_NplsM_inverse,
        }
