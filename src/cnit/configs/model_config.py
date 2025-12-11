from dataclasses import dataclass, field
import pint

from ..physics.landuse import LanduseCalculator
from ..physics.effects import (
    EffectCO2Calculator,
    EffectTemperatureCalculator,
    EffectCarbonNitrogenCouplingCalculator,
    EffectLanduseCalculator
)
from ..physics.carbon_cycle import (
    CarbonNPPLPRCalculator,
    CarbonTurnoverCalculator,
    CarbonCycleCalculator,
    CarbonCycleModel
)
from ..physics.nitrogen_cycle import (
    NitrogenPUBNFCalculator,
    NitrogenTurnoverCalculator,
    NitrogenCycleCalculator,
    NitrogenCycleModel
)

from ..utils.units import Q

calc_lu = LanduseCalculator()
calc_eff_CO2 = EffectCO2Calculator()
calc_eff_dT = EffectTemperatureCalculator()
calc_eff_CN = EffectCarbonNitrogenCouplingCalculator()
calc_eff_LU = EffectLanduseCalculator()
calc_npp_lpr = CarbonNPPLPRCalculator()
calc_c_turnover = CarbonTurnoverCalculator()
calc_c_cycle = CarbonCycleCalculator()
mdl_c_cycle = CarbonCycleModel(
    calc_c_turnover=calc_c_turnover,
    calc_c_cycle=calc_c_cycle,
)
calc_pu_bnf = NitrogenPUBNFCalculator()
calc_n_turnover = NitrogenTurnoverCalculator()
calc_n_cycle = NitrogenCycleCalculator()
mdl_n_cycle = NitrogenCycleModel(
    calc_n_turnover=calc_n_turnover,
    calc_n_cycle=calc_n_cycle,
)


@dataclass(slots=True)
class CarbonNitrogenCycleModelConfig:
    frc_rgr: pint.Quantity = calc_lu.frc_rgr
    tau_rgr: pint.Quantity = calc_lu.tau_rgr
    method_cLU: pint.Quantity = calc_lu.method_cLU
    method_nLU: pint.Quantity = calc_lu.method_nLU

    CO2ref: pint.Quantity = calc_eff_CO2.CO2ref
    CO2b: pint.Quantity = calc_eff_CO2.CO2b
    sns_CO2_log: pint.Quantity = calc_eff_CO2.sns_CO2_log
    sns_CO2_sig: pint.Quantity = calc_eff_CO2.sns_CO2_sig
    eff_CO2_sig_max: pint.Quantity = calc_eff_CO2.eff_CO2_sig_max
    method_CO2_NPP: pint.Quantity = calc_eff_CO2.method_CO2_NPP

    sns_NPP2dT: pint.Quantity = calc_eff_dT.sns_NPP2dT
    sns_NPP2dT_sig: pint.Quantity = calc_eff_dT.sns_NPP2dT_sig
    method_dT_NPP: pint.Quantity = calc_eff_dT.method_dT_NPP
    sns_LPR2dT: pint.Quantity = calc_eff_dT.sns_LPR2dT
    sns_cLP2dT: pint.Quantity = calc_eff_dT.sns_cLP2dT
    sns_cLD2dT: pint.Quantity = calc_eff_dT.sns_cLD2dT
    sns_cSR2dT: pint.Quantity = calc_eff_dT.sns_cSR2dT
    sns_PU2dT: pint.Quantity = calc_eff_dT.sns_PU2dT
    sns_BNF2dT: pint.Quantity = calc_eff_dT.sns_BNF2dT
    sns_nLP2dT: pint.Quantity = calc_eff_dT.sns_nLP2dT
    sns_nLD2dT: pint.Quantity = calc_eff_dT.sns_nLD2dT
    sns_nSR2dT: pint.Quantity = calc_eff_dT.sns_nSR2dT
    sns_nLSgas2dT: pint.Quantity = calc_eff_dT.sns_nLSgas2dT

    sns_NPP2PUdef: pint.Quantity = calc_eff_CN.sns_NPP2PUdef
    sns_cLP2PUdef: pint.Quantity = calc_eff_CN.sns_cLP2PUdef
    sns_cLD2PUdef: pint.Quantity = calc_eff_CN.sns_cLD2PUdef
    sns_cSR2PUdef: pint.Quantity = calc_eff_CN.sns_cSR2PUdef
    sns_PU2PUdef: pint.Quantity = calc_eff_CN.sns_PU2PUdef
    sns_BNF2PUdef: pint.Quantity = calc_eff_CN.sns_BNF2PUdef
    sns_nLP2PUdef: pint.Quantity = calc_eff_CN.sns_nLP2PUdef
    sns_nLD2PUdef: pint.Quantity = calc_eff_CN.sns_nLD2PUdef
    sns_nSR2PUdef: pint.Quantity = calc_eff_CN.sns_nSR2PUdef
    eff_C_PU_max: pint.Quantity = calc_eff_CN.eff_C_PU_max
    sns_PU2NPPrd: pint.Quantity = calc_eff_CN.sns_PU2NPPrd
    eff_C_BNF_max: pint.Quantity = calc_eff_CN.eff_C_BNF_max
    sns_BNF2NPPrd: pint.Quantity = calc_eff_CN.sns_BNF2NPPrd

    sns_cLP2LUrgr: pint.Quantity = calc_eff_LU.sns_cLP2LUrgr
    sns_cLD2LUrgr: pint.Quantity = calc_eff_LU.sns_cLD2LUrgr
    sns_cSR2LUrgr: pint.Quantity = calc_eff_LU.sns_cSR2LUrgr
    sns_nLP2LUrgr: pint.Quantity = calc_eff_LU.sns_nLP2LUrgr
    sns_nLD2LUrgr: pint.Quantity = calc_eff_LU.sns_nLD2LUrgr
    sns_nSR2LUrgr: pint.Quantity = calc_eff_LU.sns_nSR2LUrgr

    NPP0: pint.Quantity = calc_npp_lpr.NPP0
    LPR0: pint.Quantity = calc_npp_lpr.LPR0

    frc_NPP2P: pint.Quantity = calc_c_cycle.frc_NPP2P
    frc_NPP2L: pint.Quantity = calc_c_cycle.frc_NPP2L
    frc_cLP2L: pint.Quantity = calc_c_cycle.frc_cLP2L
    frc_cLD2S: pint.Quantity = calc_c_cycle.frc_cLD2S
    frc_cLUgrs2P: pint.Quantity = calc_c_cycle.frc_cLUgrs2P
    frc_cLUgrs2L: pint.Quantity = calc_c_cycle.frc_cLUgrs2L

    tau_CplsP: pint.Quantity = calc_c_turnover.tau_CplsP
    tau_CplsL: pint.Quantity = calc_c_turnover.tau_CplsL
    tau_CplsS: pint.Quantity = calc_c_turnover.tau_CplsS

    CplsP0: pint.Quantity = mdl_c_cycle.CplsP0
    CplsL0: pint.Quantity = mdl_c_cycle.CplsL0
    CplsS0: pint.Quantity = mdl_c_cycle.CplsS0

    PU0: pint.Quantity = calc_pu_bnf.PU0
    BNF0: pint.Quantity = calc_pu_bnf.BNF0
    NCratio_PUreq: pint.Quantity = calc_pu_bnf.NCratio_PUreq
    NetMIN0: pint.Quantity = calc_pu_bnf.NetMIN0
    sns_NetMIN2dT: pint.Quantity = calc_pu_bnf.sns_NetMIN2dT
    NetMINlt0: pint.Quantity = calc_pu_bnf.NetMINlt0
    sns_NetMINlt2AD: pint.Quantity = calc_pu_bnf.sns_NetMINlt2AD

    tau_NplsP: pint.Quantity = calc_n_turnover.tau_NplsP
    tau_NplsL: pint.Quantity = calc_n_turnover.tau_NplsL
    tau_NplsS: pint.Quantity = calc_n_turnover.tau_NplsS
    tau_NplsM: pint.Quantity = calc_n_turnover.tau_NplsM

    frc_BNF2P: pint.Quantity = calc_n_cycle.frc_BNF2P
    frc_BNF2L: pint.Quantity = calc_n_cycle.frc_BNF2L
    frc_PU2P: pint.Quantity = calc_n_cycle.frc_PU2P
    frc_PU2L: pint.Quantity = calc_n_cycle.frc_PU2L
    frc_nLP2L: pint.Quantity = calc_n_cycle.frc_nLP2L
    frc_nLD2S: pint.Quantity = calc_n_cycle.frc_nLD2S
    frc_nLSgas: pint.Quantity = calc_n_cycle.frc_nLSgas
    frc_nLUgrs2P: pint.Quantity = calc_n_cycle.frc_nLUgrs2P
    frc_nLUgrs2L: pint.Quantity = calc_n_cycle.frc_nLUgrs2L

    NplsP0: pint.Quantity = mdl_n_cycle.NplsP0
    NplsL0: pint.Quantity = mdl_n_cycle.NplsL0
    NplsS0: pint.Quantity = mdl_n_cycle.NplsS0
    NplsM0: pint.Quantity = mdl_n_cycle.NplsM0

    time0: pint.Quantity = mdl_c_cycle.time0
    switch_N: pint.Quantity = Q(1, "1")
    switch_Cpls: list = field(default_factory=lambda: [1, 1, 1])
    switch_Npls: list = field(default_factory=lambda: [1, 1, 1, 1])
