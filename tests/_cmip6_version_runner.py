"""Run default and calibrated CMIP6 cases against a CNit source tree."""

import json
from pathlib import Path
import sys

import numpy as np
import xarray as xr


source_path = Path(sys.argv[1])
forcing_path = Path(sys.argv[2])
config_path = Path(sys.argv[3])
output_path = Path(sys.argv[4])
truncate_tau_rgr = len(sys.argv) > 5 and sys.argv[5] == "truncate_tau_rgr"
sys.path.insert(0, str(source_path))

from cnit.physics.carbon_nitrogen_cycle import CarbonNitrogenCycleModel  # noqa: E402
from cnit.utils.units import Q  # noqa: E402


FORCING_UNITS = {
    "dT": "K",
    "CO2": "ppm",
    "CemsLUnet": "GtC/yr",
    "CemsLUgrs": "GtC/yr",
    "NflxAD": "GtN/yr",
    "NflxFT": "GtN/yr",
    "NemsLUnet": "GtN/yr",
    "NemsLUgrs": "GtN/yr",
    "NemsLUmin": "GtN/yr",
}
EXPERIMENTS = ("1pctCO2", "hist_ssp585")


def _load_model(values):
    params = {
        name: (
            Q(value["magnitude"], value["unit"])
            if isinstance(value, dict) and "unit" in value
            else value
        )
        for name, value in values.items()
    }
    if truncate_tau_rgr:
        params["tau_rgr"] = Q(int(params["tau_rgr"].to("yr").m), "yr")
    return CarbonNitrogenCycleModel.from_dict(params)


def _experiment_forcing(forcing, model_name, experiment_name):
    if experiment_name == "hist_ssp585":
        historical = forcing.sel(experiment="historical", model=model_name)
        ssp585 = forcing.sel(experiment="ssp585", model=model_name)
        # Historical data end in 2014; use SSP585 starting in 2015.
        return xr.where(historical.time < 2015, historical, ssp585).fillna(0)
    return forcing.sel(experiment=experiment_name, model=model_name)


model_configs = json.loads(config_path.read_text())
forcing = xr.open_dataset(forcing_path)
outputs = {}
try:
    for config_name, config_values in model_configs.items():
        model = _load_model(config_values)
        model_name = "ACCESS-ESM1-5" if config_name == "default" else config_name
        for experiment_name in EXPERIMENTS:
            experiment_forcing = _experiment_forcing(
                forcing,
                model_name,
                experiment_name,
            )
            valid_time = np.isfinite(experiment_forcing["CO2"].values)
            experiment_forcing = experiment_forcing.isel(time=valid_time).fillna(0)
            time = experiment_forcing.time.values
            run_inputs = {
                f"{name}_s": Q(experiment_forcing[name].values, unit)
                for name, unit in FORCING_UNITS.items()
            }
            result = model.run(time_axis=Q(time, "yr"), **run_inputs)
            prefix = f"{config_name}::{experiment_name}::"
            for variable, data_array in result.data_vars.items():
                outputs[f"{prefix}{variable}"] = data_array.values
            outputs[f"{prefix}__time"] = time
finally:
    forcing.close()

np.savez_compressed(output_path, **outputs)
