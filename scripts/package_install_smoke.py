"""Smoke-test a CNit installation by loading parameters and running the model."""

from pathlib import Path
import sys

import numpy as np

import cnit
from cnit import CNitModel, CNitModelConfig, Q
from cnit.configs.model_config import CarbonNitrogenCycleModelConfig


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_PATH = Path(cnit.__file__).resolve()
try:
    PACKAGE_PATH.relative_to(REPO_ROOT)
except ValueError:
    pass
else:
    raise RuntimeError(
        f"CNit was imported from the source checkout, not the installed package: "
        f"{PACKAGE_PATH}"
    )

parameter_file = Path(sys.argv[1])
calibrated_config = CarbonNitrogenCycleModelConfig.from_parameter_csv(
    parameter_file,
    model_name="ACCESS-ESM1-5",
)
calibrated_model = CNitModel.from_config(calibrated_config)
assert calibrated_model.calc_lu.frc_rgr == calibrated_config.frc_rgr

time_axis = Q(np.arange(2000, 2005), "yr")
n_steps = len(time_axis)
zero_carbon_flux = Q(np.zeros(n_steps), "GtC/yr")
zero_nitrogen_flux = Q(np.zeros(n_steps), "GtN/yr")

simulation_config = CNitModelConfig.from_dict({"time0": Q(2000, "yr")})
result = CNitModel.from_config(simulation_config).run(
    time_axis=time_axis,
    dT_s=Q(np.linspace(0, 0.2, n_steps), "K"),
    CO2_s=Q(np.linspace(280, 285, n_steps), "ppm"),
    CemsLUnet_s=zero_carbon_flux,
    CemsLUgrs_s=zero_carbon_flux,
    NflxAD_s=Q(np.full(n_steps, 0.05), "GtN/yr"),
    NflxFT_s=zero_nitrogen_flux,
    NemsLUnet_s=zero_nitrogen_flux,
    NemsLUgrs_s=zero_nitrogen_flux,
    NemsLUmin_s=zero_nitrogen_flux,
)

for variable in ("CplsP", "CplsL", "CplsS", "NplsP", "NplsL", "NplsS", "NplsM"):
    values = result[variable].data
    assert values.shape == (
        n_steps,
    ), f"Unexpected shape for {variable}: {values.shape}"
    assert np.isfinite(values).all(), f"Non-finite values in {variable}"

print(f"Installed CNit simulation OK ({PACKAGE_PATH})")
