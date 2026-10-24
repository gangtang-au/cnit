import numpy as np
import matplotlib
from pathlib import Path

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cnit.configs.model_config import CarbonNitrogenCycleModelConfig
from cnit.physics.carbon_nitrogen_cycle import CarbonNitrogenCycleModel
from cnit.physics.landuse import LanduseCalculator
from cnit.utils.units import Q


def test_single_deforestation_pulse_recovers_configured_land_carbon_fraction():
    calculator = LanduseCalculator()
    # Gross-input mode represents one directly observed deforestation event.
    calculator.method_cLU = Q(1, "1")
    years = int(calculator.tau_rgr.to("yr").m) + 1
    gross_emissions = np.zeros(years)
    gross_emissions[0] = 1.0

    fluxes = calculator.calculate_CflxLU_series(
        CemsLUnet_s=Q(np.zeros(years), "GtC/yr"),
        CemsLUgrs_s=Q(gross_emissions, "GtC/yr"),
    )

    regrowth = fluxes["CflxLUrgr"].to("GtC/yr").m
    np.testing.assert_allclose(regrowth[0], 0.0)
    np.testing.assert_allclose(
        regrowth[1:],
        calculator.frc_rgr.to("1").m / calculator.tau_rgr.to("yr").m,
    )

    recovered_carbon = fluxes["cumsum_CflxLUrgr"].to("GtC").m[-1]
    expected_recovery = calculator.frc_rgr.to("1").m * gross_emissions[0]
    np.testing.assert_allclose(recovered_carbon, expected_recovery)

    # Gross emissions retain only the one input pulse while regrowth proceeds.
    calculated_gross = fluxes["CflxLUgrs"].to("GtC/yr").m
    np.testing.assert_allclose(calculated_gross, gross_emissions)


def test_subannual_gross_flux_rounds_regrowth_and_decay_to_nearest_step_count():
    calculator = LanduseCalculator(
        frc_rgr=Q(0.8, "1"),
        tau_rgr=Q(1.4, "yr"),
        method_cLU=Q(1, "1"),
    )
    time_axis = Q(np.arange(5) * 0.5, "yr")
    gross_emissions = Q([2.0, -2.0, 0.0, 0.0, 0.0], "GtC/yr")

    fluxes = calculator.calculate_CflxLU_series(
        CemsLUnet_s=Q(np.zeros(5), "GtC/yr"),
        CemsLUgrs_s=gross_emissions,
        time_axis=time_axis,
    )

    # Each rate applies for half a year. A 1.4-year timescale becomes three steps
    # (1.5 years), and each pulse begins one interval after its input event.
    np.testing.assert_allclose(
        fluxes["CflxLUrgr"].to("GtC/yr").m,
        [0.0, 0.8 / 1.5, 0.8 / 1.5, 0.8 / 1.5, 0.0],
    )
    np.testing.assert_allclose(
        fluxes["CflxLUafst_decay"].to("GtC/yr").m,
        [0.0, 0.0, 0.2 / 1.5, 0.2 / 1.5, 0.2 / 1.5],
    )
    np.testing.assert_allclose(
        fluxes["cumsum_CflxLUrgr"].to("GtC").m[-1],
        0.8,
    )
    np.testing.assert_allclose(
        fluxes["cumsum_CflxLUafst_decay"].to("GtC").m[-1],
        0.2,
    )


def test_subannual_net_flux_reconstructs_gross_flux_with_regrowth():
    calculator = LanduseCalculator(
        frc_rgr=Q(0.8, "1"),
        tau_rgr=Q(1.25, "yr"),
        method_cLU=Q(0, "1"),
    )
    time_axis = Q(np.arange(4) * 0.5, "yr")
    net_emissions = Q([2.0, 0.0, 0.0, 0.0], "GtC/yr")

    fluxes = calculator.calculate_CflxLU_series(
        CemsLUnet_s=net_emissions,
        CemsLUgrs_s=Q(np.zeros(4), "GtC/yr"),
        time_axis=time_axis,
    )

    gross = fluxes["CflxLUgrs"].to("GtC/yr").m
    regrowth = fluxes["CflxLUrgr"].to("GtC/yr").m
    np.testing.assert_allclose(gross, net_emissions.m + regrowth)
    assert np.all(np.isfinite(gross))


def test_model_land_carbon_stock_approaches_gross_land_use_adjusted_equilibrium():
    cfg_dict = dict(
        frc_rgr=Q(0.8, "1"),
        tau_rgr=Q(100, "yr"),
        method_cLU=Q(1, "1"),
        CO2ref=Q(284.875, "ppm"),
        NPP0=Q(50, "GtC/yr"),
        CplsP0=Q(500, "GtC"),
        CplsL0=Q(50, "GtC"),
        CplsS0=Q(1450, "GtC"),
        time0=Q(1700, "yr"),
        frc_NPP2P=Q(1, "1"),
        frc_NPP2L=Q(0, "1"),
        frc_cLP2L=Q(1, "1"),
        frc_cLD2S=Q(1, "1"),
        frc_cLUgrs2P=Q(1, "1"),
        frc_cLUgrs2L=Q(0, "1"),
        switch_N=Q(0, "1"),
        sns_cLP2LUrgr=Q(0, "yr/GtC"),
        # These turnover times make the specified pool sizes an equilibrium:
        # 500 / 50 = 10 yr, 50 / 50 = 1 yr, and 1450 / 50 = 29 yr.
        tau_CplsP=Q(10, "yr"),
        tau_CplsL=Q(1, "yr"),
        tau_CplsS=Q(29, "yr"),
    )
    config = CarbonNitrogenCycleModelConfig.from_dict(cfg_dict)
    pulse_start = 10
    pulse_duration = 10
    years = np.arange(int(config.time0.m), int(config.time0.m) + 500)
    zero_carbon_flux = Q(np.zeros(years.size), "GtC/yr")
    zero_nitrogen_emissions = Q(np.zeros(years.size), "GtN/yr")
    nitrogen_deposition = Q(np.full(years.size, 0.2), "GtN/yr")
    nitrogen_fertilizer = Q(np.full(years.size, 0.1), "GtN/yr")
    temperature = Q(np.zeros(years.size), "K")
    co2 = Q(np.full(years.size, config.CO2ref.to("ppm").m), "ppm")

    common_forcing = {
        "time_axis": Q(years, "yr"),
        "dT_s": temperature,
        "CO2_s": co2,
        "CemsLUnet_s": zero_carbon_flux,
        "NflxAD_s": nitrogen_deposition,
        "NflxFT_s": nitrogen_fertilizer,
        "NemsLUnet_s": zero_nitrogen_emissions,
        "NemsLUgrs_s": zero_nitrogen_emissions,
        "NemsLUmin_s": zero_nitrogen_emissions,
    }

    control = CarbonNitrogenCycleModel.from_config(config).run(
        CemsLUgrs_s=zero_carbon_flux,
        **common_forcing,
    )
    gross_land_use = np.zeros(years.size)
    gross_land_use[pulse_start : pulse_start + pulse_duration] = 1.0
    disturbed = CarbonNitrogenCycleModel.from_config(config).run(
        CemsLUgrs_s=Q(gross_land_use, "GtC/yr"),
        **common_forcing,
    )

    control_land_carbon = sum(control[pool] for pool in ("CplsP", "CplsL", "CplsS"))
    disturbed_land_carbon = sum(disturbed[pool] for pool in ("CplsP", "CplsL", "CplsS"))
    pulse_effect = (disturbed_land_carbon - control_land_carbon).values

    np.testing.assert_allclose(
        control_land_carbon.values,
        control_land_carbon.values[0],
        rtol=1e-5,
        atol=1e-7,
    )

    figure, (stock_axis, difference_axis, flux_axis, npp_axis) = plt.subplots(
        4, 1, figsize=(9, 13), sharex=True, constrained_layout=True
    )
    stock_axis.plot(years, control_land_carbon.values, label="No land-use emissions")
    stock_axis.plot(
        years,
        disturbed_land_carbon.values,
        label="Gross emissions: 1 GtC/yr for 10 years",
    )
    stock_axis.axvline(
        years[pulse_start + pulse_duration],
        color="tab:gray",
        linestyle="--",
        label="Gross emissions forcing ends",
    )
    stock_axis.set_ylabel("Total land carbon (GtC)")
    stock_axis.legend()
    stock_axis.grid(True, alpha=0.3)

    smoothing_window = 100
    difference_axis.plot(
        years, pulse_effect, color="tab:red", label="Annual difference"
    )
    difference_axis.axhline(0, color="black", linewidth=0.8)
    equilibrium_effect = np.mean(pulse_effect[-smoothing_window:])
    stock_axis.set_xlabel("Year")
    difference_axis.set_xlabel("Year")
    difference_axis.set_ylabel("Pulse minus control (GtC)")
    difference_axis.legend()
    difference_axis.grid(True, alpha=0.3)

    flux_axis.plot(
        years,
        disturbed["CflxLUgrs"].values,
        label="Gross land-use emissions",
    )
    flux_axis.plot(
        years,
        disturbed["CflxLUrgr"].values,
        label="Land regrowth",
    )
    flux_axis.axhline(0, color="black", linewidth=0.8)
    flux_axis.set_xlabel("Year")
    flux_axis.set_ylabel("Land-use flux (GtC/yr)")
    flux_axis.legend()
    flux_axis.grid(True, alpha=0.3)

    npp_axis.plot(
        years,
        control["CflxNPP"].values,
        label="No land-use emissions",
    )
    npp_axis.plot(
        years,
        disturbed["CflxNPP"].values,
        label="Gross emissions: 1 GtC/yr for 10 years",
    )
    npp_axis.set_xlabel("Year")
    npp_axis.set_ylabel("NPP (GtC/yr)")
    npp_axis.legend()
    npp_axis.grid(True, alpha=0.3)
    npp_axis.set_xlim(years[0], years[-1])
    figure.suptitle("Land carbon response to 10 years of gross land-use emissions")
    figure.savefig(
        Path(__file__).with_name("landuse_recovery.png"),
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(figure)

    np.testing.assert_allclose(pulse_effect[0], 0.0, atol=1e-8)
    assert pulse_effect[pulse_start + pulse_duration] < 0
    assert abs(pulse_effect[50]) < abs(pulse_effect[pulse_start + pulse_duration])
    # Over 500 years, the response is still approaching the new equilibrium.
    previous_window_mean = np.mean(
        pulse_effect[-2 * smoothing_window : -smoothing_window]
    )
    assert equilibrium_effect < previous_window_mean
    assert abs(equilibrium_effect - previous_window_mean) < 0.75
    assert equilibrium_effect < 0
    # Gross emissions remove 10 GtC; 80% regrows, leaving a 2 GtC permanent loss.
    expected_permanent_loss = -(1 - config.frc_rgr.to("1").m) * pulse_duration
    np.testing.assert_allclose(equilibrium_effect, expected_permanent_loss, atol=0.1)
    np.testing.assert_allclose(
        disturbed["cumsum_CflxLUrgr"].values[-1],
        pulse_duration * config.frc_rgr.to("1").m,
        rtol=1e-6,
    )
