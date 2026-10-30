from unittest.mock import patch

import numpy as np
import pytest
import xarray as xr

from cnit.configs.experiment_config import CarbonNitrogenCycleExperimentConfig
from cnit.configs.model_config import CarbonNitrogenCycleModelConfig
from cnit.physics.carbon_nitrogen_cycle import CarbonNitrogenCycleModel
from cnit.utils.units import Q


def _experiment(name, switch_N):
    time_axis = Q([1850, 1851], "yr")
    return CarbonNitrogenCycleExperimentConfig(
        name=name,
        switch_N=Q(switch_N, "1"),
        time_axis=time_axis,
        dT_s=Q(np.zeros(2), "K"),
        CO2_s=Q(np.full(2, 280.0), "ppm"),
        CemsLUnet_s=Q(np.zeros(2), "GtC/yr"),
        CemsLUgrs_s=Q(np.zeros(2), "GtC/yr"),
        NflxAD_s=Q(np.zeros(2), "GtN/yr"),
        NflxFT_s=Q(np.zeros(2), "GtN/yr"),
        NemsLUnet_s=Q(np.zeros(2), "GtN/yr"),
        NemsLUgrs_s=Q(np.zeros(2), "GtN/yr"),
        NemsLUmin_s=Q(np.zeros(2), "GtN/yr"),
    )


@pytest.mark.parametrize(
    "model_switch",
    [0, 1],
)
def test_run_experiments_uses_experiment_switch_without_mutating_model(model_switch):
    model = CarbonNitrogenCycleModel.from_config(CarbonNitrogenCycleModelConfig())
    model.switch_N = Q(model_switch, "1")
    model_switch_before = model.switch_N
    observed_switches = []
    observed_models = []

    def run_stub(experiment_model, **_):
        observed_switches.append(experiment_model.switch_N.to("1").m)
        observed_models.append(experiment_model)
        return xr.Dataset(
            {"value": ("time", np.array([len(observed_switches)], dtype=float))},
            coords={"time": [1850]},
        )

    experiments = [_experiment("disabled", 0), _experiment("enabled", 1)]
    with patch.object(CarbonNitrogenCycleModel, "run", run_stub):
        result = model.run_experiments(experiments)

    assert observed_switches == [0, 1]
    assert model.switch_N == model_switch_before
    assert all(experiment_model is not model for experiment_model in observed_models)
    assert all(
        experiment_model.calc_lu is model.calc_lu
        for experiment_model in observed_models
    )
    assert result.experiment.values.tolist() == ["disabled", "enabled"]
