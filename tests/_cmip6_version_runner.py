"""Run calibrated CMIP6 cases against a selected CNit source tree."""

import csv
import json
from pathlib import Path
import sys

import numpy as np
import xarray as xr


source_path = Path(sys.argv[1])
parameter_path = Path(sys.argv[2])
forcing_path = Path(sys.argv[3])
config_dir = Path(sys.argv[4])
output_path = Path(sys.argv[5])
revision = sys.argv[6]
sys.path.insert(0, str(source_path))

from cnit.configs.model_config import CarbonNitrogenCycleModelConfig  # noqa: E402
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


def _load_model(model_name):
    if revision == "current":
        config = CarbonNitrogenCycleModelConfig.from_parameter_csv(
            parameter_path,
            model_name=model_name,
        )
        return CarbonNitrogenCycleModel.from_config(config)

    values = json.loads((config_dir / f"{model_name}.json").read_text())
    params = {
        name: (
            Q(value["magnitude"], value["unit"])
            if isinstance(value, dict) and "unit" in value
            else value
        )
        for name, value in values.items()
    }
    return CarbonNitrogenCycleModel.from_dict(params)


def _experiment_forcing(forcing, model_name, experiment_name):
    if experiment_name == "hist_ssp585":
        historical = forcing.sel(experiment="historical", model=model_name)
        ssp585 = forcing.sel(experiment="ssp585", model=model_name)
        # Historical data end in 2014; use SSP585 starting in 2015.
        return xr.where(historical.time < 2015, historical, ssp585).fillna(0)
    return forcing.sel(experiment=experiment_name, model=model_name)


with parameter_path.open(newline="", encoding="utf-8") as parameter_file:
    model_names = next(csv.reader(parameter_file))[2:]

forcing = xr.open_dataset(forcing_path)
outputs = {}
try:
    for model_name in model_names:
        model = _load_model(model_name)
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
            prefix = f"{model_name}::{experiment_name}::"
            for variable, data_array in result.data_vars.items():
                outputs[f"{prefix}{variable}"] = data_array.values
            outputs[f"{prefix}__time"] = time
finally:
    forcing.close()

np.savez_compressed(output_path, **outputs)
