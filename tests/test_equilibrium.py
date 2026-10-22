import numpy as np

from cnit.configs.model_config import CarbonNitrogenCycleModelConfig
from cnit.physics.carbon_nitrogen_cycle import CarbonNitrogenCycleModel
from cnit.utils.units import Q


def test_equilibrium_initialization_keeps_carbon_and_nitrogen_pools_stable():
    config = CarbonNitrogenCycleModelConfig()
    model = CarbonNitrogenCycleModel.from_config(config)

    equilibrium_forcing = {
        "dT_s": Q([0.0], "K"),
        "CO2_s": Q([config.CO2ref.to("ppm").m], "ppm"),
        "CemsLUnet_s": Q([0.0], "GtC/yr"),
        "CemsLUgrs_s": Q([0.0], "GtC/yr"),
        "NflxAD_s": Q([0.2], "GtN/yr"),
        "NflxFT_s": Q([0.1], "GtN/yr"),
        "NemsLUnet_s": Q([0.0], "GtN/yr"),
        "NemsLUgrs_s": Q([0.0], "GtN/yr"),
        "NemsLUmin_s": Q([0.0], "GtN/yr"),
    }
    equilibrium = model.calculate_mdl_eqm0(**equilibrium_forcing)

    for pool in ("P", "L", "S"):
        tau = equilibrium[f"tau_Cpls{pool}_inverse"]
        assert np.isfinite(tau.m) and tau.m > 0
        setattr(
            model.mdl_c_cycle.calc_c_turnover,
            f"tau_Cpls{pool}",
            tau,
        )
    for pool in ("P", "L", "S", "M"):
        tau = equilibrium[f"tau_Npls{pool}_inverse"]
        assert np.isfinite(tau.m) and tau.m > 0
        setattr(
            model.mdl_n_cycle.calc_n_turnover,
            f"tau_Npls{pool}",
            tau,
        )

    years = np.arange(1850, 1861)
    result = model.run(
        time_axis=Q(years, "yr"),
        dT_s=Q(np.zeros(years.size), "K"),
        CO2_s=Q(np.full(years.size, config.CO2ref.to("ppm").m), "ppm"),
        CemsLUnet_s=Q(np.zeros(years.size), "GtC/yr"),
        CemsLUgrs_s=Q(np.zeros(years.size), "GtC/yr"),
        NflxAD_s=Q(np.full(years.size, 0.2), "GtN/yr"),
        NflxFT_s=Q(np.full(years.size, 0.1), "GtN/yr"),
        NemsLUnet_s=Q(np.zeros(years.size), "GtN/yr"),
        NemsLUgrs_s=Q(np.zeros(years.size), "GtN/yr"),
        NemsLUmin_s=Q(np.zeros(years.size), "GtN/yr"),
    )

    for pool in ("CplsP", "CplsL", "CplsS", "NplsP", "NplsL", "NplsS", "NplsM"):
        values = result[pool].values
        np.testing.assert_allclose(values, values[0], rtol=1e-5, atol=1e-7)
