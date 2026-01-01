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
    """MAGICC's coupled terrestrial carbon-nitrogen cycle model (CNit).

    This is the top-level model class that integrates the carbon and nitrogen cycles
    with their environmental drivers and feedback mechanisms.  The model simulates the
    full terrestrial biogeochemistry system including:

    - Carbon cycle dynamics (3 pools: plant, litter, soil)
    - Nitrogen cycle dynamics (4 pools: plant, litter, soil, mineral)
    - Carbon-nitrogen coupling through:

        * Nitrogen limitation effects on carbon processes (NPP, decomposition)
        * Carbon-driven nitrogen demand (stoichiometric requirements)
        * Nitrogen deficit feedbacks
    - Environmental responses:

        * CO2 fertilization effects on NPP
        * Temperature effects on all processes
        * Land use change impacts (deforestation, afforestation, regrowth)

    The model operates on annual timesteps and solves the carbon and nitrogen cycles
    iteratively to capture their bidirectional coupling. A key feature is the explicit
    calculation of nitrogen deficit (difference between nitrogen required and available
    for plant uptake), which creates feedback loops regulating both carbon and nitrogen
    cycling.

    **Model Architecture:**

    The model consists of several interconnected calculator components:

    - **Land Use Calculator** (:class:`LanduseCalculator`): Processes land use change
      emissions and regrowth fluxes for both carbon and nitrogen
    - **Effect Calculators**:  Compute environmental modifiers
        * :class:`EffectCO2Calculator`: CO2 fertilization effects
        * :class:`EffectTemperatureCalculator`: Temperature effects on all processes
        * :class:`EffectCarbonNitrogenCouplingCalculator`: C-N coupling effects
        * :class:`EffectLanduseCalculator`: Land use effects on turnover rates
    - **Flux Calculators**:
        * :class:`CarbonNPPLPRCalculator`: NPP and LPR with environmental modifiers
        * :class:`NitrogenPUBNFCalculator`: PU and BNF with nitrogen deficit feedbacks
    - **Cycle Models**:
        * :class:`CarbonCycleModel`: Solves carbon pool dynamics (3 pools)
        * :class:`NitrogenCycleModel`: Solves nitrogen pool dynamics (4 pools)

    **Coupling Strategy:**

    The carbon and nitrogen cycles are coupled through an iterative approach:

    1. Calculate potential NPP (without nitrogen limitation)
    2. Calculate nitrogen requirement based on NPP and stoichiometry
    3. Calculate nitrogen availability from mineralization and external inputs
    4. Calculate nitrogen deficit (required - available)
    5. Apply nitrogen limitation effect to get actual NPP
    6. Calculate nitrogen uptake and fixation based on actual NPP
    7. Solve carbon cycle with nitrogen effects
    8. Solve nitrogen cycle with carbon-driven demand

    This approach ensures that carbon and nitrogen dynamics are mutually consistent
    while avoiding numerical instability from simultaneous solution.

    **Configuration:**

    The model can be initialized in two ways:

    - :meth:`from_dict`: Build from a parameter dictionary
    - :meth:`from_config`: Build from a :class:`CarbonNitrogenCycleModelConfig` dataclass

    Both methods automatically distribute parameters to the appropriate sub-components.
    """

    # Calculators
    calc_lu:  LanduseCalculator
    """Calculator for land use fluxes. 
    
    Processes land use change emissions (gross deforestation, net emissions) and
    regrowth for both carbon and nitrogen cycles.  Tracks cumulative land use changes
    for calculating land use effects on productivity. 
    """

    calc_eff_CO2: EffectCO2Calculator
    """Calculator for CO2 fertilization effects. 
    
    Computes the enhancement of NPP due to elevated atmospheric CO2 concentration
    through increased photosynthetic efficiency and water use efficiency. 
    """

    calc_eff_dT: EffectTemperatureCalculator
    """Calculator for temperature effects. 
    
    Computes temperature modifiers for all carbon and nitrogen processes including
    NPP, LPR, decomposition, mineralization, plant uptake, biological fixation, and
    gaseous losses.
    """

    calc_eff_CN: EffectCarbonNitrogenCouplingCalculator
    """Calculator for carbon-nitrogen coupling effects.
    
    Computes nitrogen limitation effects on carbon processes (NPP, decomposition)
    and carbon availability effects on nitrogen processes (plant uptake, fixation)
    based on nitrogen deficit and NPP changes.
    """

    calc_eff_LU: EffectLanduseCalculator
    """Calculator for land use effects. 
    
    Computes land use change effects on process rates (turnover, decomposition)
    based on land use history (regrowth, deforestation, afforestation). 
    """

    calc_npp_lpr: CarbonNPPLPRCalculator
    """Calculator for net primary production (NPP) and litter production respiration (LPR).
    
    Computes carbon fixation by plants (NPP) and associated autotrophic respiration
    (LPR) with environmental modifiers (CO2, temperature, nitrogen, land use).
    """

    calc_pu_bnf: NitrogenPUBNFCalculator
    """Calculator for plant nitrogen uptake (PU) and biological nitrogen fixation (BNF).
    
    Computes nitrogen inputs to organic pools through plant uptake from mineral pool
    and biological fixation from atmosphere, with environmental modifiers and nitrogen
    deficit feedbacks. 
    """

    # Model components
    mdl_c_cycle: CarbonCycleModel
    """Carbon cycle model component.
    
    Solves the carbon pool dynamics (plant, litter, soil) using ODE integration.
    Receives nitrogen limitation effects from nitrogen cycle and provides NPP for
    nitrogen demand calculation.
    """

    mdl_n_cycle: NitrogenCycleModel
    """Nitrogen cycle model component.
    
    Solves the nitrogen pool dynamics (plant, litter, soil, mineral) using two-stage
    ODE integration. Receives carbon-driven nitrogen demand and provides nitrogen
    limitation effects to carbon cycle.
    """

    # Model switches
    switch_N:  pint. Quantity = field(
        default=Q(1, "1"),
    )
    """Switch to enable/disable nitrogen cycle [dimensionless]. 
    
    Controls whether nitrogen limitation effects are applied to carbon cycle: 
    - 1: Full carbon-nitrogen coupling (default)
    - 0: Carbon cycle only (no nitrogen limitation)
    
    When set to 0, the model runs as a carbon-only model with no nitrogen feedbacks.
    """

    @classmethod
    def from_dict(cls, params: dict) -> "CarbonNitrogenCycleModel":
        """
        Build a fully configured CarbonNitrogenCycleModel from a parameter dictionary.

        This factory method creates all calculator and model components, then distributes
        parameters from the dictionary to the appropriate components based on attribute
        names.

        Parameters
        ----------
        params
            Parameter dictionary for configuring calculators and sub-models.  Keys should
            match attribute names of calculator/model components.  Unknown keys raise an
            error.

        Returns
        -------
        CarbonNitrogenCycleModel
            Fully configured model instance ready to run.

        Raises
        ------
        ValueError
            If params contains keys that don't match any component attribute names.

        Examples
        --------
        >>> params = {
        ...     'NPP0': Q(50, 'GtC/yr'),
        ...     'tau_CplsP':  Q(10, 'yr'),
        ...     'NCratio_PUreq': Q(0.02, 'GtN/GtC'),
        ...     'switch_N': Q(1, '1'),
        ... }
        >>> model = CarbonNitrogenCycleModel.from_dict(params)
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

        This is a convenience wrapper around : meth:`from_dict` that accepts a
        configuration dataclass instead of a dictionary.

        Parameters
        ----------
        config
            Configuration dataclass containing all model parameters.

        Returns
        -------
        CarbonNitrogenCycleModel
            Fully configured model instance ready to run.

        See Also
        --------
        from_dict : Build model from parameter dictionary
        CarbonNitrogenCycleModelConfig : Configuration dataclass definition
        """
        cfg_dict = asdict(config)
        return cls.from_dict(cfg_dict)

    def run_experiments(self, experiments: List["CarbonNitrogenCycleExperimentConfig"]):
        """"
        Run multiple experiments with the model.

        This method facilitates running multiple model scenarios with different
        configurations and/or forcing data. Results from all experiments are
        concatenated into a single dataset with an 'experiment' dimension.

        Parameters
        ----------
        experiments
            List of experiment configurations.  Each experiment specifies forcing
            data, nitrogen switch setting, and a name for identification.

        Returns
        -------
        xr.Dataset
            Dataset containing results for all experiments, with an 'experiment'
            dimension for comparison.  Variables are aligned across experiments
            using outer join (missing values filled with NaN).

        Examples
        --------
        >>> exp1 = CarbonNitrogenCycleExperimentConfig(
        ...     name="control",
        ...     switch_N=Q(1, "1"),
        ...     # ...  forcing data ...
        ... )
        >>> exp2 = CarbonNitrogenCycleExperimentConfig(
        ...     name="no_nitrogen",
        ...     switch_N=Q(0, "1"),
        ...     # ... forcing data ...
        ... )
        >>> results = model.run_experiments([exp1, exp2])
        >>> control_npp = results.sel(experiment="control")["CflxNPP"]
        >>> no_n_npp = results.sel(experiment="no_nitrogen")["CflxNPP"]
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
        """
        Run the coupled carbon-nitrogen cycle model.

        This is the main simulation method that:

        1. Processes land use change forcing
        2. Calculates potential and actual NPP with nitrogen limitation
        3. Calculates nitrogen deficit and feedbacks
        4. Solves carbon cycle with nitrogen effects
        5. Solves nitrogen cycle with carbon-driven demand (if nitrogen cycle exists)
        6. Returns combined results

        Parameters
        ----------
        time_axis
            Time points for simulation [yr]. Must start at or after initialization
            time (time0) specified in model components.
        dT_s
            Time series of temperature anomaly relative to preindustrial [K].
        CO2_s
            Time series of atmospheric CO2 concentration [ppm].
        CemsLUnet_s
            Time series of net land use carbon emissions [GtC/yr].  Includes both
            gross emissions and regrowth.
        CemsLUgrs_s
            Time series of gross land use carbon emissions [GtC/yr].  From
            deforestation and land conversion only.
        NflxAD_s
            Time series of atmospheric nitrogen deposition flux [GtN/yr].  Includes
            both natural and anthropogenic sources.
        NflxFT_s
            Time series of nitrogen fertilizer application flux [GtN/yr].
            Anthropogenic nitrogen input from agriculture.
        NemsLUnet_s
            Time series of net land use nitrogen emissions [GtN/yr].  Includes both
            gross emissions and regrowth.
        NemsLUgrs_s
            Time series of gross land use nitrogen emissions from organic pools [GtN/yr].
            From deforestation and land conversion.
        NemsLUmin_s
            Time series of land use nitrogen emissions from mineral pool [GtN/yr].
            Direct mineral nitrogen losses during land conversion.

        Returns
        -------
        xr.Dataset
            Complete simulation results containing:

            **Carbon Cycle Variables:**

            - Pool sizes:  CplsP, CplsL, CplsS, CplsPLS
            - Primary fluxes: CflxNPP, CflxNPPpot, CflxLPR, CflxLPRpot
            - Turnover fluxes: CflxLP, CflxLD, CflxSR
            - Partitioned fluxes: CflxNPP2P/L/S, CflxLP2L/S, CflxLD2S/A, CflxLUgrs2P/L/S
            - Derived fluxes: CflxRH, CflxNetPLS
            - Land use:  CflxLUgrs, CflxLUrgr, cumulative land use changes

            **Nitrogen Cycle Variables** (if switch_N=1):

            - Pool sizes: NplsP, NplsL, NplsS, NplsM, NplsPLS, NplsPLSM
            - Primary fluxes: NflxPU, NflxBNF, NflxAD, NflxFT
            - Turnover fluxes: NflxLP, NflxLD, NflxSR, NflxLS
            - Partitioned fluxes: NflxPU2P/L/S, NflxBNF2P/L/S, NflxLP2L/S,
              NflxLD2S/M, NflxLUgrs2P/L/S
            - Derived fluxes:  NflxNetMIN, NflxNetPLSM
            - Nitrogen budget: NflxPUreq, NflxPUavail, NflxPUavail_netMIN, NflxPUdef
            - Land use: NflxLUgrs, NflxLUmin

        Notes
        -----
        **Solution Sequence:**

        The model solves carbon and nitrogen cycles in sequence to maintain consistency:

        1. **Calculate potential NPP:**  NPP without nitrogen limitation (ε_N = 1)
        2. **Calculate nitrogen budget:**  Required vs available nitrogen for plant uptake
        3. **Calculate nitrogen deficit:**  PUdef = PUreq - PUavail
        4. **Apply nitrogen limitation:**  Calculate actual NPP with nitrogen effect
        5. **Solve carbon cycle:**  With nitrogen effects on all processes
        6. **Solve nitrogen cycle:**  With carbon-driven demand

        **Nitrogen Limitation:**

        When switch_N = 1, nitrogen limitation affects carbon processes through:

        ..  math::

            NPP_{actual} = NPP_{potential} \\times \\epsilon_{N(NPP)}

        where the nitrogen effect depends on nitrogen deficit:

        .. math::

            \\epsilon_{N(NPP)} = e^{s_{NPP2PUdef} \\times PUdef}

        Negative deficit (nitrogen surplus) can enhance NPP, while positive deficit
        (nitrogen limitation) reduces NPP.

        **Model Switches:**

        The behavior can be controlled through switches:
        - switch_N = 0: Carbon-only mode (no nitrogen limitation)
        - switch_Cpls = [0,0,0]:  Freeze all carbon pools
        - switch_Npls = [0,0,0,0]:  Freeze all nitrogen pools

        See Also
        --------
        calculate_CflxNPPLPR_NflxPUBNF: Core coupling calculation
        :py:meth:`cnit.physics.carbon_cycle.CarbonCycleModel.run`: Carbon cycle solution
        :py:meth:`cnit.physics.nitrogen_cycle.NitrogenCycleModel.run`: Nitrogen cycle solution
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
        Calculate coupled carbon and nitrogen fluxes with nitrogen limitation feedback.

        This is the core coupling method that implements the nitrogen limitation feedback
        on NPP.   It follows a two-step approach:

        1. **Calculate potential fluxes** (without nitrogen limitation):

           - Potential NPP and LPR assuming no nitrogen constraint
           - Nitrogen requirement based on potential NPP and stoichiometry
           - Nitrogen availability from mineralization and external inputs
           - Nitrogen deficit (required - available)

        2. **Calculate actual fluxes** (with nitrogen limitation):

           - Apply nitrogen limitation effect to get actual NPP and LPR
           - Calculate nitrogen uptake and fixation based on actual conditions

        This approach ensures that nitrogen limitation is calculated based on the system's
        potential carbon uptake, while the actual fluxes reflect nitrogen constraints.

        Parameters
        ----------
        dT
            Temperature anomaly [K]. Used for temperature effects on all processes.
        CO2
            Atmospheric CO2 concentration [ppm]. Used for CO2 fertilization effect.
        CflxLUrgr
            Land use regrowth carbon flux [GtC/yr]. Adds to baseline NPP.
        cumsum_CflxLUrgr
            Cumulative land use regrowth [GtC]. Used for land use effects.
        cumsum_CflxLUdfst
            Cumulative deforestation [GtC]. Used for land use effects.
        cumsum_CflxLUafst
            Cumulative afforestation [GtC]. Used for land use effects.
        cumsum_CflxLUafst_decay
            Cumulative afforestation with decay [GtC]. Used for land use effects.
        NflxAD
            Atmospheric nitrogen deposition flux [GtN/yr]. External nitrogen input.
        NflxFT
            Nitrogen fertilizer flux [GtN/yr]. Anthropogenic nitrogen input.

        Returns
        -------
        Dictionary containing the following time series:

            **Carbon Fluxes:**

            - CflxNPP:  Actual (nitrogen-limited) net primary production [GtC/yr]
            - CflxLPR: Actual (nitrogen-limited) litter production respiration [GtC/yr]
            - CflxNPPpot: Potential (non-nitrogen-limited) NPP [GtC/yr]
            - CflxLPRpot: Potential (non-nitrogen-limited) LPR [GtC/yr]

            **Nitrogen Fluxes:**

            - NflxPU: Plant nitrogen uptake [GtN/yr]
            - NflxBNF: Biological nitrogen fixation [GtN/yr]

            **Nitrogen Budget:**

            - NflxPUreq:  Nitrogen required for plant uptake [GtN/yr]
            - NflxPUavail: Nitrogen available for plant uptake [GtN/yr]
            - NflxPUavail_netMIN: N available from net mineralization [GtN/yr]
            - NflxPUdef: Nitrogen deficit for plant uptake [GtN/yr]

        Notes
        -----
        **Nitrogen Limitation Feedback:**

        The nitrogen limitation on NPP follows:

        .. math::

            NPP_{actual} = NPP_{potential} \\times \\epsilon_{N(NPP)}

        where the nitrogen effect depends on nitrogen deficit:

        .. math::

            \\epsilon_{N(NPP)} = e^{s_{NPP2PUdef} \\times PUdef}

        **Nitrogen Budget:**

        Nitrogen requirement is based on stoichiometry:

        .. math::

            PUreq = NPP_{potential} \\times NCratio_{PUreq}

        Nitrogen availability includes:

        .. math::

            PUavail = NetMIN + AD

        where NetMIN has temperature-sensitive and long-term components.

        **Switch Behavior:**

        When nitrogen cycle is disabled (switch_N = 0 or NplsP0 = 0):
        - NPP = NPPpot (no nitrogen limitation)
        - All nitrogen fluxes set to NaN
        - Model runs in carbon-only mode

        See Also
        --------
        :py:meth:`cnit.physics.carbon_cycle.CarbonNPPLPRCalculator.calculate`: NPP and LPR calculation
        :py:meth:`cnit.physics.nitrogen_cycle.NitrogenPUBNFCalculator.calculate_NflxPU_req_avail_def`: Nitrogen plant uptake budget
        :py:meth:`cnit.physics.effects.EffectCarbonNitrogenCouplingCalculator.calculate_eff_N_NPP`: Nitrogen effect
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
        """
                Calculate equilibrium initial state for the model.

                This method calculates the equilibrium pool sizes and turnover times that
                would result from steady forcing conditions. It is typically used to initialize
                the model in equilibrium with preindustrial conditions before running transient
                simulations.

                The method uses the first timestep values of all forcing to calculate steady-state
                fluxes, then uses inverse turnover calculations to determine equilibrium pool
                sizes and turnover times.

                Parameters
                ----------
                dT_s
                    Time series of temperature anomaly [K].  Only first value used.
                CO2_s
                    Time series of atmospheric CO2 concentration [ppm]. Only first value used.
                CemsLUnet_s
                    Time series of net land use carbon emissions [GtC/yr]. Only first value used.
                CemsLUgrs_s
                    Time series of gross land use carbon emissions [GtC/yr].  Only first value used.
                NflxAD_s
                    Time series of atmospheric nitrogen deposition [GtN/yr]. Only first value used.
                NflxFT_s
                    Time series of nitrogen fertilizer [GtN/yr]. Only first value used.
                NemsLUnet_s
                    Time series of net land use nitrogen emissions [GtN/yr]. Only first value used.
                NemsLUgrs_s
                    Time series of gross land use nitrogen emissions [GtN/yr]. Only first value used.
                NemsLUmin_s
                    Time series of land use mineral nitrogen emissions [GtN/yr]. Only first value used.

                Returns
                -------
                Dict[str, pint.Quantity]
                    Dictionary containing equilibrium values for:

                    **Carbon Equilibrium:**

                    - CflxLP_inverse, CflxLD_inverse, CflxSR_inverse: Turnover fluxes
                    - CflxLP2L_inverse, CflxLP2S_inverse:  Litter production partitioning
                    - CflxLD2S_inverse, CflxLD2A_inverse: Litter decomposition partitioning
                    - tau_CplsP_inverse, tau_CplsL_inverse, tau_CplsS_inverse: Turnover times

                    **Nitrogen Equilibrium** (if nitrogen cycle active):

                    - NflxLP_inverse, NflxLD_inverse, NflxSR_inverse, NflxLS_inverse: Turnover fluxes
                    - NflxLP2L_inverse, NflxLP2S_inverse: Litter production partitioning
                    - NflxLD2S_inverse, NflxLD2M_inverse: Litter decomposition partitioning
                    - tau_NplsP_inverse, tau_NplsL_inverse, tau_NplsS_inverse, tau_NplsM_inverse:
                      Turnover times

                Notes
                -----
                **Equilibrium Concept:**

                At equilibrium, pool sizes are constant (dC/dt = 0, dN/dt = 0), which means:

                - Inputs = Outputs for each pool
                - Fluxes are determined by forcing conditions
                - Pool sizes and turnover times are mutually consistent

                **Usage:**

                Typical workflow for equilibrium initialization:

                ..  code-block:: python

                    # Calculate equilibrium
                    eqm = model.calculate_mdl_eqm0(
                        dT_s=preindustrial_temp,
                        CO2_s=preindustrial_CO2,
                        # ... other preindustrial forcing ...
                    )

                    # Set initial pool sizes (if desired)
                    model.mdl_c_cycle.CplsP0 = eqm_pool_sizes["CplsP"]

                    # Set equilibrium turnover times
                    model.mdl_c_cycle.calc_c_turnover. tau_CplsP = eqm["tau_CplsP_inverse"]

                See Also
                --------
                calculate_eqm0 : Core equilibrium calculation
                inverse_carbon_turnover : Carbon equilibrium equations
                inverse_nitrogen_turnover : Nitrogen equilibrium equations
        """


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
        """
                Calculate equilibrium state from specified fluxes.

                This method calculates equilibrium pool sizes and turnover times given a set
                of fluxes and environmental conditions. It assumes steady state (dC/dt = 0,
                dN/dt = 0) and solves for the pool sizes and turnover times that satisfy
                mass balance.

                Parameters
                ----------
                dT
                    Temperature anomaly [K]. Used for calculating environmental effects.
                CflxNPP
                    Net primary production [GtC/yr]. Carbon input to system.
                CflxLPR
                    Litter production respiration [GtC/yr].  Autotrophic respiration.
                CflxLUgrs
                    Gross land use carbon emissions [GtC/yr]. Carbon output from system.
                CflxLUrgr
                    Land use regrowth [GtC/yr]. Used for land use effects.
                NflxPUdef
                    Nitrogen deficit [GtN/yr]. Used for nitrogen effects.
                NflxPU
                    Plant nitrogen uptake [GtN/yr]. N transfer from mineral to organic pools.
                NflxBNF
                    Biological nitrogen fixation [GtN/yr]. External N input.
                NflxAD
                    Atmospheric nitrogen deposition [GtN/yr].  External N input.
                NflxFT
                    Nitrogen fertilizer [GtN/yr]. External N input.
                NflxLUgrs
                    Gross land use nitrogen emissions [GtN/yr]. N output from organic pools.
                NflxLUmin
                    Land use mineral nitrogen losses [GtN/yr]. N output from mineral pool.

                Returns
                -------
                Dict[str, pint.Quantity]
                    Dictionary containing equilibrium values (see : meth:`calculate_mdl_eqm0`
                    for complete list).

                Notes
                -----
                **Equilibrium Mass Balance:**

                For each pool, at equilibrium (dX/dt = 0):

                .. math::

                    \\sum Inputs = \\sum Outputs

                This allows solving for either:
                - Pool sizes (given fluxes and turnover times), or
                - Turnover times (given fluxes and pool sizes)

                The inverse methods solve for turnover times from the equilibrium condition.

                **Environmental Effects:**

                All environmental effects are calculated at the specified conditions (dT, CO2, etc.)
                and applied to determine equilibrium turnover rates.

                See Also
                --------
                inverse_carbon_turnover :  Solve carbon equilibrium equations
                inverse_nitrogen_turnover : Solve nitrogen equilibrium equations
        """
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
        """
                Calculate turnover flux from mass balance.

                This helper method calculates the turnover flux from a pool given the mass
                balance equation:

                .. math::

                    \\frac{dPool}{dt} = Influx - Turnover - Outflux

                Rearranging for turnover:

                .. math::

                    Turnover = Influx - Outflux - \\frac{dPool}{dt}

                At equilibrium (dPool/dt = 0):

                .. math::

                    Turnover = Influx - Outflux

                Parameters
                ----------
                influx
                    Total input to pool [GtC/yr or GtN/yr].
                outflux
                    Total output from pool (excluding turnover) [GtC/yr or GtN/yr].
                pool_size_change
                    Rate of change of pool size [GtC/yr or GtN/yr].  Zero at equilibrium.

                Returns
                -------
                    Turnover flux from pool [GtC/yr or GtN/yr].

                Notes
                -----
                This method is used by both carbon and nitrogen equilibrium calculations.
                For equilibrium calculations, pool_size_change = 0.
        """
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
        Calculate equilibrium carbon fluxes and turnover times from mass balance.

        This method solves the carbon cycle equilibrium equations to determine:

        1. Turnover fluxes (cLP, cLD, cSR) from mass balance
        2. Partitioned fluxes using allocation fractions
        3. Equilibrium turnover times from pool sizes and fluxes

        The method works backwards from pool mass balance equations to infer the
        turnover fluxes that satisfy equilibrium conditions.

        Parameters
        ----------
        CflxNPP
            Net primary production [GtC/yr]. Carbon input flux.
        CflxLPR
            Litter production respiration [GtC/yr].  Autotrophic respiration.
        CflxLUgrs
            Gross land use carbon emissions [GtC/yr].  Carbon output flux.
        eff_dT_cLP, eff_dT_cLD, eff_dT_cSR
            Temperature effects on carbon processes [dimensionless].
        eff_N_cLP, eff_N_cLD, eff_N_cSR
            Nitrogen effects on carbon processes [dimensionless].
        eff_LU_cLP, eff_LU_cLD, eff_LU_cSR
            Land use effects on carbon processes [dimensionless].
        CplsP, CplsL, CplsS
            Current carbon pool sizes [GtC].  If None, uses initial values from model.
        dCplsP_dt, dCplsL_dt, dCplsS_dt
            Pool size changes [GtC/yr]. Zero for equilibrium calculation.

        Returns
        -------
            Dictionary containing:

            **Turnover Fluxes:**

            - CflxLP_inverse: Litter production [GtC/yr]
            - CflxLD_inverse:  Litter decomposition [GtC/yr]
            - CflxSR_inverse:  Soil respiration [GtC/yr]

            **Partitioned Fluxes:**

            - CflxLP2L_inverse, CflxLP2S_inverse: LP allocation
            - CflxLD2S_inverse, CflxLD2A_inverse: LD allocation

            **Equilibrium Turnover Times:**

            - tau_CplsP_inverse: Plant pool turnover time [yr]
            - tau_CplsL_inverse:  Litter pool turnover time [yr]
            - tau_CplsS_inverse: Soil pool turnover time [yr]

        Notes
        -----
        **Mass Balance Equations:**

        For each pool at equilibrium:

        .. math::

            \\frac{dC_P}{dt} = NPP2P - LPR - cLP - cLUgrs_P = 0

            \\frac{dC_L}{dt} = NPP2L + cLP2L - cLD - cLUgrs_L = 0

            \\frac{dC_S}{dt} = NPP2S + cLP2S + cLD2S - cSR - cLUgrs_S = 0

        These are solved sequentially to find cLP, cLD, and cSR.

        **Turnover Time Calculation:**

        From first-order kinetics:

        .. math::

            Flux = \\frac{Pool}{\\tau} \\times \\prod Effects

        Rearranging:

        .. math::

            \\tau = \\frac{Pool}{Flux} \\times \\prod Effects

        See Also
        --------
        inverse_turnover_flux: Mass balance calculation
        :py:meth:`cnit.physics.carbon_cycle.CarbonTurnoverCalculator.calculate`: Forward turnover calculation
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
        """
                Calculate equilibrium nitrogen fluxes and turnover times from mass balance.

                This method solves the nitrogen cycle equilibrium equations to determine:
                1. Turnover fluxes (nLP, nLD, nSR, nLS) from mass balance
                2. Partitioned fluxes using allocation fractions
                3. Equilibrium turnover times from pool sizes and fluxes

                Similar to : meth:`inverse_carbon_turnover` but includes the mineral nitrogen
                pool and its associated fluxes.

                Parameters
                ----------
                NflxPU
                    Plant nitrogen uptake [GtN/yr]. Mineral to organic transfer.
                NflxBNF
                    Biological nitrogen fixation [GtN/yr]. External input.
                NflxAD
                    Atmospheric deposition [GtN/yr].  External input.
                NflxFT
                    Fertilizer application [GtN/yr]. External input.
                NflxLUgrs
                    Gross land use nitrogen emissions [GtN/yr]. Output from organic pools.
                NflxLUmin
                    Land use mineral nitrogen losses [GtN/yr]. Output from mineral pool.
                eff_dT_nLP, eff_dT_nLD, eff_dT_nSR, eff_dT_nLSgas
                    Temperature effects on nitrogen processes [dimensionless].
                eff_N_nLP, eff_N_nLD, eff_N_nSR
                    Carbon-nitrogen coupling effects [dimensionless].
                eff_LU_nLP, eff_LU_nLD, eff_LU_nSR
                    Land use effects on nitrogen processes [dimensionless].
                NplsP, NplsL, NplsS, NplsM
                    Current nitrogen pool sizes [GtN].  If None, uses initial values from model.
                dNplsP_dt, dNplsL_dt, dNplsS_dt, dNplsM_dt
                    Pool size changes [GtN/yr]. Zero for equilibrium calculation.

                Returns
                -------
                Dict[str, pint.Quantity]
                    Dictionary containing:

                    **Turnover Fluxes:**
                    - NflxLP_inverse, NflxLD_inverse, NflxSR_inverse, NflxLS_inverse

                    **Partitioned Fluxes:**
                    - NflxLP2L_inverse, NflxLP2S_inverse:  LP allocation
                    - NflxLD2S_inverse, NflxLD2M_inverse: LD allocation

                    **Equilibrium Turnover Times:**
                    - tau_NplsP_inverse, tau_NplsL_inverse, tau_NplsS_inverse, tau_NplsM_inverse

                Notes
                -----
                **Mass Balance Equations:**

                For equilibrium in each pool:

                .. math::

                    \\frac{dN_P}{dt} = BNF2P + PU2P - nLP - nLUgrs_P = 0

                    \\frac{dN_L}{dt} = BNF2L + PU2L + nLP2L - nLD - nLUgrs_L = 0

                    \\frac{dN_S}{dt} = BNF2S + PU2S + nLP2S + nLD2S - nSR - nLUgrs_S = 0

                    \\frac{dN_M}{dt} = AD + FT + nLD2M + nSR - PU - nLS - nLUmin = 0

                These are solved sequentially to find nLP, nLD, nSR, and nLS.

                See Also
                --------
                inverse_carbon_turnover :  Analogous carbon calculation
                :py:meth:`cnit.physics.nitrogen_cycle.NitrogenTurnoverCalculator.calculate_nLPLDSR`: Forward nitrogen turnover
                :py:meth:`cnit.physics.nitrogen_cycle.NitrogenTurnoverCalculator.calculate_nLS`: Forward mineral N losses
        """
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
