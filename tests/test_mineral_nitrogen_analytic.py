"""Accuracy checks for the analytical mineral nitrogen pool update."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cnit.physics.nitrogen_cycle import (
    NitrogenCycleCalculator,
    NitrogenCycleModel,
    NitrogenTurnoverCalculator,
)
from cnit.utils.units import Q


@pytest.mark.parametrize(
    "tau_npls_m, gas_effect, relative_tolerance",
    [
        (0.02, [1.0] * 6, 1e-7),
        (100_000.0, [1.0] * 6, 1e-7),
        (0.02, np.linspace(0.99, 1.01, 6), 1e-2),
    ],
    ids=("fast-turnover", "small-decay-limit", "slowly-varying-turnover"),
)
def test_analytic_mineral_nitrogen_matches_high_accuracy_ode_reference(
    tau_npls_m, gas_effect, relative_tolerance
):
    """Compare the interval update with a tight ODE solution of the same balance.

    The cases cover fast turnover, the small-decay limit, and modestly varying
    temperature effects.
    """
    time = np.arange(1850.0, 1856.0)
    gas_effect = np.asarray(gas_effect, dtype=float)
    uptake = np.array([0.04, 0.06, 0.05, 0.08, 0.07, 0.09])
    fixation = np.array([0.02, 0.03, 0.035, 0.025, 0.04, 0.03])
    deposition = np.array([0.10, 0.11, 0.09, 0.12, 0.10, 0.13])
    fertilizer = np.array([0.01, 0.02, 0.01, 0.03, 0.02, 0.01])
    zeros = np.zeros(time.size)
    ones = np.ones(time.size)

    turnover = NitrogenTurnoverCalculator(
        tau_NplsM=Q(tau_npls_m, "yr"),
        frc_nLSgas=Q(0.5, "1"),
    )
    model = NitrogenCycleModel(
        calc_n_turnover=turnover,
        calc_n_cycle=NitrogenCycleCalculator(),
    )
    result = model.run(
        time_axis=Q(time, "yr"),
        NflxPU_s=Q(uptake, "GtN/yr"),
        NflxBNF_s=Q(fixation, "GtN/yr"),
        NflxAD_s=Q(deposition, "GtN/yr"),
        NflxFT_s=Q(fertilizer, "GtN/yr"),
        NflxLUgrs_s=Q(zeros, "GtN/yr"),
        NflxLUmin_s=Q(zeros, "GtN/yr"),
        eff_dT_nLP_s=Q(ones, "1"),
        eff_dT_nLD_s=Q(ones, "1"),
        eff_dT_nSR_s=Q(ones, "1"),
        eff_dT_nLSgas_s=Q(gas_effect, "1"),
        eff_N_nLP_s=Q(ones, "1"),
        eff_N_nLD_s=Q(ones, "1"),
        eff_N_nSR_s=Q(ones, "1"),
        eff_LU_nLP_s=Q(ones, "1"),
        eff_LU_nLD_s=Q(ones, "1"),
        eff_LU_nSR_s=Q(ones, "1"),
    )

    # Reconstruct the mineral-pool source term from the organic-pool trajectory.
    litter = result.n_state["NplsL"].to("GtN").m
    soil = result.n_state["NplsS"].to("GtN").m
    net_mineralization = (
        litter
        / turnover.tau_NplsL.to("yr").m
        * (1 - model.calc_n_cycle.frc_nLD2S.to("1").m)
        + soil / turnover.tau_NplsS.to("yr").m
    )
    net_inflow = deposition + fertilizer + net_mineralization - uptake

    def mineral_pool_rhs(t, state):
        inflow = np.interp(t, time, net_inflow)
        temperature_effect = np.interp(t, time, gas_effect)
        loss_rate = (
            turnover.frc_nLSgas.to("1").m * temperature_effect
            + (1 - turnover.frc_nLSgas.to("1").m)
        ) / tau_npls_m
        return [inflow - loss_rate * state[0]]

    reference = solve_ivp(
        mineral_pool_rhs,
        (time[0], time[-1]),
        [model.NplsM0.to("GtN").m],
        t_eval=time,
        method="DOP853",
        rtol=1e-11,
        atol=1e-13,
        max_step=min(np.diff(time).min() / 10, tau_npls_m / 5),
    )
    assert reference.success

    np.testing.assert_allclose(
        result.n_state["NplsM"].to("GtN").m,
        reference.y[0],
        rtol=relative_tolerance,
        atol=1e-9,
    )
